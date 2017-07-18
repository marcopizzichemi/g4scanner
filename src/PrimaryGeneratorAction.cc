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
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(ConfigFile& config)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
//    ,fConfig(config)
{
  //G4int n_particle = 1;
  fParticleGun = new G4GeneralParticleSource();

  //create a messenger for this class
  //fGunMessenger = new PrimaryGeneratorMessenger(this);

  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

  fParticleGun->GetCurrentSource()->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);

  //source
  //x position
  sourcex = config.read<double>("sourcex",0);
  G4cout << "Gamma source x position [mm]: " << sourcex << G4endl;
  //y position of source
  sourcey = config.read<double>("sourcey",0);
  G4cout << "Gamma source y position [mm]: " << sourcey << G4endl;
  //x position of source
  sourcez = config.read<double>("sourcez",0);
  G4cout << "Gamma source z position [mm]: " << sourcez << G4endl;
//   distance = config.read<double>("distance");
//   G4cout << "Distance of source from back ESR [mm]: " << distance << G4endl;

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
//   greaseBack = config.read<double>("greaseBack");
//   glassBack = config.read<double>("glassBack");
//   airBack = config.read<double>("airBack");
  //multiply for units [mm]
  crystalz   = crystalz*mm;
//   greaseBack = greaseBack*mm;
//   glassBack  = glassBack*mm;
//   airBack    = airBack*mm;
//   G4double fakeAir = 0.1*mm;  //fixed to 0.1, it has no physical meaning in our simulation
//   G4double sourcez = -(distance + (crystalz/2.0) + greaseBack + glassBack + airBack + fakeAir);
//   G4cout << "Gamma source z position [mm]: " << sourcez << G4endl;
  //energy of gammas
  energy = config.read<double>("energy");
  G4cout << "Energy of gamma source [keV]: " << energy << G4endl;
  //gamma direction
//   direction = config.read<int>("direction");
//   if(direction == 0)
//   {
//     G4cout << "Gammas shoot parallel to the z axis" << G4endl;
//   }
//   else
//   {
//     G4cout << "Gammas shoot randomly towards the matrix" << G4endl;
//   }
  //find the x and y of matrix
  crystalx = config.read<double>("crystalx");
  crystaly = config.read<double>("crystaly");
  ncrystalx = config.read<int>("ncrystalx");
  ncrystaly = config.read<int>("ncrystaly");
  esrThickness = config.read<double>("esrThickness");

  //SOURCE = 1mm diameter sphere

  // set energy distribution
  G4SPSEneDistribution *eneDist = fParticleGun->GetCurrentSource()->GetEneDist() ;
  eneDist->SetEnergyDisType("Mono"); // or gauss
  eneDist->SetMonoEnergy(energy*keV);

  // set position distribution
  G4SPSPosDistribution *posDist = fParticleGun->GetCurrentSource()->GetPosDist();
  posDist->SetPosDisType("Volume");  // or Point,Plane,Volume,Beam
  posDist->SetPosDisShape("Sphere");
  posDist->SetRadius(0.5*mm);
  posDist->SetCentreCoords(G4ThreeVector(sourcex*mm,sourcey*mm,sourcez*mm));

  // set angular distribution
  G4SPSAngDistribution *angDist = fParticleGun->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection( G4ThreeVector(0., 0., 1.) );


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4double halfDiagonal =  sqrt(pow((crystalx + esrThickness) * ncrystalx,2.0) + pow( (crystaly + esrThickness) * ncrystaly ,2.0)) / 2.0;
  //we assume for this test that the 2 modules have the same distance from the center of setup
  double distance = sqrt(plate_x[0]*plate_x[0] + plate_z[0]*plate_z[0]) - crystalz/2.0;
  G4double angleLimit = atan(halfDiagonal / distance);
  //find the limit for acos
  double acosMin = cos(angleLimit);
  //so acos will have to be generated uniformely between acosMin and +1
  double randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
  //double randomNum =  G4UniformRand()*2.0  -1.0; //random num between -1 and 1
  theta = acos(randomNum);
  phi = G4UniformRand() * 2.0 * CLHEP::pi;

  G4SPSAngDistribution *angDist = fParticleGun->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta)));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  theta = theta + CLHEP::pi;
  angDist->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*sin(phi),sin(theta)*cos(phi),cos(theta)));
  fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
