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
#include "G4GeneralParticleSourceCustom.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"

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


  fParticleGun = new G4GeneralParticleSourceCustom();

  //definition of particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");
  fParticleGun->SetParticleTime(0.0*ns);

  //LOOP ON SOURCES: DEFAULT KINEMATICS
  //fix shape and dimensions, now they are all the same
  for(unsigned int i=0; i<sourcex.size(); i++)
  {
    fParticleGun->AddaSource(1);
    fParticleGun->GetCurrentSource()->SetParticleDefinition(particle);

    // set energy distribution
    G4SPSEneDistribution *eneDist = fParticleGun->GetCurrentSource()->GetEneDist() ;
    eneDist->SetEnergyDisType("Mono");
    eneDist->SetMonoEnergy(energy*keV);

    // set position distribution
    G4SPSPosDistribution *posDist = fParticleGun->GetCurrentSource()->GetPosDist();
    posDist->SetPosDisType("Volume");
    posDist->SetPosDisShape("Sphere");
    //FIXME
    posDist->SetRadius(0.5*mm);
    posDist->SetCentreCoords(G4ThreeVector(sourcex[i]*mm,sourcey[i]*mm,sourcez[i]*mm));

    // set angular distribution
    G4SPSAngDistribution *angDist = fParticleGun->GetCurrentSource()->GetAngDist();
    angDist->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    theta.push_back(0);
    phi.push_back(0);
  }

  //hack to make the loop on the sources
  fParticleGun->DeleteaSource(0);

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

  fParticleGun->SetCurrentSourceto(i);

  G4double halfDiagonal =  sqrt(pow((crystalx + esrThickness) * ncrystalx,2.0) + pow( (crystaly + esrThickness) * ncrystaly ,2.0)) / 2.0;
  //we take into account that the source might not be in the center
  distance_plate0.push_back(sqrt(plate_x[0]*plate_x[0] + (plate_z[0]-sourcez[i])*(plate_z[0]-sourcez[i])) - crystalz/2.0);
  distance_plate1.push_back(sqrt(plate_x[1]*plate_x[1] + (plate_z[1]-sourcez[i])*(plate_z[1]-sourcez[i])) - crystalz/2.0);

  //choose the minimum distance, so we are sure the gamma will hit all the crystals
  G4double angleLimit = atan(halfDiagonal / std::min(distance_plate0[i],distance_plate1[i]));
  //find the limit for acos
  double acosMin = cos(angleLimit);
  //so acos will have to be generated uniformely between acosMin and +1
  double randomNum =  G4UniformRand()*(1.0 - acosMin)  + (acosMin);
  theta[i] = acos(randomNum);

  phi[i] = G4UniformRand() * 2.0 * CLHEP::pi;

  G4SPSAngDistribution *angDist = fParticleGun->GetCurrentSource()->GetAngDist();
  //generate first gamma
  angDist->SetParticleMomentumDirection(G4ThreeVector(sin(theta[i])*sin(phi[i]),sin(theta[i])*cos(phi[i]),cos(theta[i])));
  fParticleGun->GeneratePrimaryVertex(anEvent);

  //generate second gamma (back to back)
  theta[i] = theta[i] + CLHEP::pi;
  angDist->SetParticleMomentumDirection(G4ThreeVector(sin(theta[i])*sin(phi[i]),sin(theta[i])*cos(phi[i]),cos(theta[i])));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
