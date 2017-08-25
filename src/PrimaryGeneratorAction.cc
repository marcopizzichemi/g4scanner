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

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4SPSRandomGenerator.hh"
#include <assert.h>     /* assert */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(ConfigFile& config)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{

  //READ FROM CONFIG FILE

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

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //random number of the source that will produce gammas
  //!! ONLY ONE SOURCE AT A TIME.. TO FIXME?
  int nSources = sourcex.size();
  int i = round(G4UniformRand()*nSources-0.5);

  //position of the source: randomize around chosen center in sphere volume
  radius = 1.5; //mm
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

  //debug
//   G4cout << "Number of source: " << i
//          << "\nPosition of source: " << (fParticleGun->GetParticlePosition()).getX() << "\t"
//          << (fParticleGun->GetParticlePosition()).getY() << "\t"
//          << (fParticleGun->GetParticlePosition()).getZ() << "\t" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
