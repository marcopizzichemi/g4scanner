//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
// #include "g4matrixPrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "CreateTree.hh"
#include "G4SPSRandomGenerator.hh"

#include <iostream>
#include <fstream>
#include <assert.h>     /* assert */


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixPrimaryGeneratorAction::g4matrixPrimaryGeneratorAction(ConfigFile& config)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
//    ,fConfig(config)
{

  //READ FROM CONFIG FILE

  G4double sphereRadius = config.read<double>("sphereRadius");
  radius = sphereRadius; //mm

  //read position and orientation of plates in space
  plates = config.read<int>("plates",2);
  // so we just put 0 to all y
  for(int i = 0 ; i < plates; i++)
    plate_y.push_back(0);
  modules = config.read<int>("modules",1);
  //   G4cout << "Plates in detector: " << plates << G4endl;
  //   G4cout << "Modules per plate: " << modules << G4endl;

  plate_x_s    = config.read<std::string>("plateCenterX","0");
  //   plate_y_s    = config.read<std::string>("plateCenterY",0);
  plate_z_s    = config.read<std::string>("plateCenterZ","0");
  rotation_s   = config.read<std::string>("plateRotation","0");
  config.split( plate_x_f, plate_x_s, "," );
  config.split( plate_z_f, plate_z_s, "," );
  config.split( rotation_f, rotation_s, "," );
  for(unsigned int i = 0 ; i < plate_x_f.size() ; i++)
  {
    config.trim(plate_x_f[i]);
    plate_x.push_back(atof(plate_x_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < plate_z_f.size() ; i++)
  {
    config.trim(plate_z_f[i]);
    plate_z.push_back(atof(plate_z_f[i].c_str()));
  }
  for(unsigned int i = 0 ; i < rotation_f.size() ; i++)
  {
    config.trim(rotation_f[i]);
    rotation.push_back(atof(rotation_f[i].c_str()));
  }
  assert( plates == plate_x.size() );
  assert( plate_x.size() == plate_z.size() );
  assert( plate_z.size() == rotation.size() );

  //retrieve dimensions of the other elements and calculate the position of the source
  crystalz = config.read<double>("crystalz");
  crystalz   = crystalz*mm;

  //energy of gammas
  energy = config.read<double>("energy");
  G4cout << "Energy of gamma source [keV]: " << energy << G4endl;

  //find the x and y of matrix
  crystalx = config.read<double>("crystalx");
  crystaly = config.read<double>("crystaly");
  ncrystalx = config.read<int>("ncrystalx");
  ncrystaly = config.read<int>("ncrystaly");
  esrThickness = config.read<double>("esrThickness");

  //read source positions
  sourcex_s    = config.read<std::string>("sourcex","0");
  sourcey_s    = config.read<std::string>("sourcey","0");
  sourcez_s    = config.read<std::string>("sourcez","0");
  config.split( sourcex_f, sourcex_s, "," );
  config.split( sourcey_f, sourcey_s, "," );
  config.split( sourcez_f, sourcez_s, "," );
  for(unsigned int i = 0 ; i < sourcex_f.size() ; i++)
  {
    config.trim(sourcex_f[i]);
    sourcex.push_back(atof(sourcex_f[i].c_str()));
    G4cout << "Gamma source x position [mm]: " << sourcex[0] << G4endl;
  }
  for(unsigned int i = 0 ; i < sourcey_f.size() ; i++)
  {
    config.trim(sourcey_f[i]);
    sourcey.push_back(atof(sourcey_f[i].c_str()));
    G4cout << "Gamma source y position [mm]: " << sourcey[0] << G4endl;
  }

  for(unsigned int i = 0 ; i < sourcez_f.size() ; i++)
  {
    config.trim(sourcez_f[i]);
    sourcez.push_back(atof(sourcez_f[i].c_str()));
  }

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  //definition of particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  fParticleGun->SetParticleTime(0.0*ns);
  fParticleGun->SetParticleDefinition(particle);


  //LOOP ON SOURCES: DEFAULT KINEMATICS
  for(unsigned int i=0; i<sourcex.size(); i++)
  {
    //set default position
    fParticleGun->SetParticlePosition(G4ThreeVector(sourcex[i],sourcey[i],sourcez[i]));

    //set default momentum direction
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    //set energy
    fParticleGun->SetParticleEnergy(energy*keV);

    theta.push_back(0);
    phi.push_back(0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

g4matrixPrimaryGeneratorAction::~g4matrixPrimaryGeneratorAction()
{
  delete fParticleGun;
  //   delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void g4matrixPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //random number of the source that will produce gammas
  //!! ONLY ONE SOURCE AT A TIME.. FIXME?
  int nSources = sourcex.size();
  int i = round(G4UniformRand()*nSources-0.5);

  //position of the source: randomize around chosen center in sphere volume


  randx = radius*2.;
  randy = radius*2.;
  randz = radius*2.;
  G4SPSRandomGenerator* posRndm = new G4SPSRandomGenerator();

  while(((randx*randx)+(randy*randy)+(randz*randz)) > (radius*radius))
  {
    randx = ((posRndm->GenRandX())*2.*radius) - radius;
    randy = ((posRndm->GenRandY())*2.*radius) - radius;
    randz = ((posRndm->GenRandZ())*2.*radius) - radius;
  }
  fParticleGun->SetParticlePosition(G4ThreeVector((sourcex[i]+randx)*mm,(sourcey[i]+randy)*mm,(sourcez[i]+randz)*mm));

  //emission of the source

  G4double Diagonal =  sqrt(pow((crystalx + esrThickness) * ncrystalx,2.0) + pow( (crystaly + esrThickness) * ncrystaly ,2.0));
  //we take into account that the source might not be in the center
  distance_plate0.push_back(sqrt((plate_x[0]-sourcex[i])*(plate_x[0]-sourcex[i]) + (plate_y[0]-sourcey[i])*(plate_y[0]-sourcey[i]) + (plate_z[0]-sourcez[i])*(plate_z[0]-sourcez[i])) - crystalz/2.0);
  distance_plate1.push_back(sqrt((plate_x[1]-sourcex[i])*(plate_x[1]-sourcex[i]) + (plate_y[1]-sourcey[i])*(plate_y[1]-sourcey[i]) + (plate_z[1]-sourcez[i])*(plate_z[1]-sourcez[i])) - crystalz/2.0);

  //choose the minimum distance, so we are sure the gamma will hit all the crystals
  G4double angleLimit = atan(Diagonal / std::max(distance_plate0[i],distance_plate1[i]));
  //find the limit for acos
  double acosMin = cos(angleLimit);
  //so acos will have to be generated uniformely between acosMin and +1
  double randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
  theta[i] = acos(randomNum);


  phi[i] = G4UniformRand() * 2.0 * CLHEP::pi;


  //number of plate (0 or 1) furthest to the source (we're only interested in coincidences)
  int Nmaxplate = 0;
  if (std::max(distance_plate0[i],distance_plate1[i]) == distance_plate1[i]) { Nmaxplate = 1; }

  //random vector to cover all the plate
  G4ThreeVector randomDirection = G4ThreeVector(sin(theta[i])*sin(phi[i]),sin(theta[i])*cos(phi[i]),cos(theta[i]));

  //generate first gamma
  fParticleGun->SetParticleMomentumDirection(pow(-1,Nmaxplate)*randomDirection);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //generate second gamma (back to back)
  fParticleGun->SetParticleMomentumDirection(-(pow(-1,Nmaxplate)*randomDirection));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // std::ofstream myfile;
  // myfile.open ("SourcePosition.txt",std::ios::app);
  // myfile << EventSourcex << " " << EventSourcey << " " << EventSourcez << " " << EventTargetx << " " << EventTargety << " " << EventTargetz << std::endl;
  // myfile.close();
  // float ScintillatorsLengthX = (crystalx + esrThickness) * ncrystalx;
  // float ScintillatorsLengthY = (crystaly + esrThickness) * ncrystaly;
  // float ScintillatorsLengthZ = crystalz;

  // if(backgroudSimulation)
  // {
  //   // G4cout << "----------------- O -----------------" << G4endl;
  //   G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  //   if (particle == G4ChargedGeantino::ChargedGeantino()) {
  //     //lu176
  //     G4int Z = 71, A = 176;
  //     G4double ionCharge   = 0.*eplus;
  //     G4double excitEnergy = 0.*keV;
  //
  //     G4ParticleDefinition* ion
  //     = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
  //     fParticleGun->SetParticleDefinition(ion);
  //     fParticleGun->SetParticleCharge(ionCharge);
  //   }
  //   EventTargetx = EventTargety = EventTargetz = 0;  // target does not makes sense in this configuration
  //   // randomized position
  //   // inside the matrix volume, but only in LYSO
  //   // generate random position inside the matrix volume,
  //   //generate position
  //   // bool isInLYSO = false;
  //   // while(!isInLYSO)
  //   // {
  //     EventSourcex = (G4UniformRand()*ScintillatorsLengthX) - (ScintillatorsLengthX/2.0);
  //     EventSourcey = (G4UniformRand()*ScintillatorsLengthY) - (ScintillatorsLengthY/2.0);
  //     EventSourcez = (G4UniformRand()*ScintillatorsLengthZ) - (ScintillatorsLengthZ/2.0);
  //     // check if it is in the space occupied by LYSO or not
  //     // EventSourcez is not a problem here
  //
  //   // }
  //   //check that it's inside the LYSO material?
  //   G4cout << EventSourcex << " " << EventSourcex << " " << EventSourcez << G4endl;
  //   fParticleGun->SetParticlePosition(G4ThreeVector(EventSourcex,EventSourcey,EventSourcez));
  //
  //   //create vertex
  //   //
  //   fParticleGun->GeneratePrimaryVertex(anEvent);
  // }
  // else
  // {
  //   // source as a sphere
  //   // if sphereRadius > 0 is specified in the config file, the source is no more a point but a sphere centered in sourcex,sourcey,sourcez
  //   if(sphereRadius > 0)
  //   {
  //     // G4cout << "Using a sphere as a gamma source" << G4endl;
  //     // G4cout << "Sphere radius = " << sphereRadius << " mm" << G4endl;
  //     float sphereCenterX = sourcex;
  //     float sphereCenterY = sourcey;
  //     float sphereCenterZ = sourcez;
  //     //generate a random point in the sphere
  //     float alpha = G4UniformRand()*2*CLHEP::pi;
  //     float ran = G4UniformRand();
  //     float beta = acos(2.0*ran-1.0);
  //     float ranradius = cbrt(G4UniformRand());
  //     float EventRadius =  sphereRadius*ranradius;
  //     //now redefine the sourcex,sourcey,sourcez on the basis of this
  //     EventSourcex = sphereCenterX + EventRadius*cos(beta);
  //     EventSourcey = sphereCenterY + EventRadius*sin(beta)*cos(alpha);
  //     EventSourcez = sphereCenterZ + EventRadius*sin(beta)*sin(alpha);
  //     //find the target point in the matrix, to derive the direction vector
  //     //crystals are always centered in (0,0,0), so they cover a volume defined by they physical coordinates times their number
  //     // float ScintillatorsLengthX = (crystalx + esrThickness) * ncrystalx;
  //     // float ScintillatorsLengthY = (crystaly + esrThickness) * ncrystaly;
  //     // float ScintillatorsLengthZ = crystalz;
  //     EventTargetx = (G4UniformRand()*ScintillatorsLengthX) - (ScintillatorsLengthX/2.0);
  //     EventTargety = (G4UniformRand()*ScintillatorsLengthY) - (ScintillatorsLengthY/2.0);
  //     EventTargetz = (G4UniformRand()*ScintillatorsLengthZ) - (ScintillatorsLengthZ/2.0);
  //   }
  //   else
  //   {
  //     EventSourcex = sourcex;
  //     EventSourcey = sourcey;
  //     EventSourcez = sourcez;
  //   }
  //
  //   //save source coordinates
  //   // G4cout << EventSourcex << " " << EventSourcey << " " << EventSourcez << G4endl;
  //   // CreateTree::Instance()->SourceX = (Float_t) EventSourcex;
  //   // CreateTree::Instance()->SourceY = (Float_t) EventSourcey;
  //   // CreateTree::Instance()->SourceZ = (Float_t) EventSourcez;
  //
  //
  //
  //   //now direction
  //   G4ThreeVector directionVector;
  //   if(direction == 0)
  //   {
  //     directionVector = G4ThreeVector(0.,0.,1.);
  //   }
  //   else
  //   {
  //     if(sphereRadius == 0)
  //     {
  //       G4double halfDiagonal;
  //       if(userIrradiationDiagonal == 0)
  //       halfDiagonal = sqrt(pow((crystalx + esrThickness) * ncrystalx,2.0) + pow( (crystaly + esrThickness) * ncrystaly ,2.0)) / 2.0;
  //       else
  //       halfDiagonal = userIrradiationDiagonal/2.0;
  //       G4double angleLimit = atan(halfDiagonal / distance);
  //       //theta = (G4UniformRand() * 2.0*angleLimit) - angleLimit; //WRONG!!!! this would make cos(theta) uniform, not sin(theta)d(theta)
  //       //find the limit for acos
  //       double acosMin = cos(angleLimit);
  //       //so acos will have to be generated uniformely between acosMin and +1
  //       double randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
  //       theta = acos(randomNum);
  //       phi = G4UniformRand() * 2.0 * CLHEP::pi;
  //       directionVector = G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta) );
  //     }
  //     else
  //     {
  //       // G4ThreeVector source = G4ThreeVector(EventSourcex,EventSourcey,EventSourcez);
  //       // G4ThreeVector target = G4ThreeVector(EventTargetx,EventTargety,EventTargetz);
  //
  //       G4ThreeVector tempVector = G4ThreeVector(EventTargetx-EventSourcex,EventTargety-EventSourcey,EventTargetz-EventSourcez) ;
  //       float module = sqrt(pow(tempVector.getX(),2)+pow(tempVector.getY(),2)+pow(tempVector.getX(),2));
  //       directionVector = tempVector /module;
  //     }
  //   }
  //
  //   //TODO save direction vector
  //   // CreateTree::Instance()->SourceMomentumX = (Float_t) directionVector.getX();
  //   // CreateTree::Instance()->SourceMomentumY = (Float_t) directionVector.getY();
  //   // CreateTree::Instance()->SourceMomentumZ = (Float_t) directionVector.getZ();
  //
  //   // CreateTree::Instance()->SourceY = EventSourcey;
  //   // CreateTree::Instance()->SourceZ = EventSourcez;
  //   fParticleGun->SetParticlePosition(G4ThreeVector(EventSourcex,EventSourcey,EventSourcez));
  //   fParticleGun->SetParticleMomentumDirection(directionVector);
  //   //G4cout << "Theta = " << theta << G4endl;
  //   //G4cout << "Phi = " << phi << G4endl;
  //
  //   // if(direction == 0) // parallel to z axis, towards matrix
  //   // {
  //   //   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //   // }
  //   // else if(sphereRadius == 0) // from pointlike source, with a given aperture
  //   // {
  //   //   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta) ));
  //   // }
  //   // else // from spherical source, towards the matrix
  //   // {
  //   //   G4ThreeVector source = G4ThreeVector(EventSourcex,EventSourcey,EventSourcez);
  //   //   G4ThreeVector target = G4ThreeVector(EventTargetx,EventTargety,EventTargetz);
  //   //   directionVector = target - source;
  //   //
  //   // }
  //   fParticleGun->GeneratePrimaryVertex(anEvent);
  //
  //   // DEBUG output to check source position and target position
  //   // std::cout << EventSourcex << " "
  //   //           << EventSourcey << " "
  //   //           << EventSourcez << G4endl;
  //
  // }
  // std::ofstream myfile;
  // myfile.open ("SourcePosition.txt",std::ios::app);
  // myfile << EventSourcex << " " << EventSourcey << " " << EventSourcez << " " << EventTargetx << " " << EventTargety << " " << EventTargetz << std::endl;
  // myfile.close();


}
