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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "ConfigFile.h"
#include <vector>

class G4GeneralParticleSource;
class G4Event;
class PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(ConfigFile& config);
    virtual ~PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event*);

//     void SetOptPhotonPolar();
//     void SetOptPhotonPolar(G4double);

  private:
    G4GeneralParticleSource* fParticleGun;
    PrimaryGeneratorMessenger* fGunMessenger;
    std::vector<G4double> sourcex, sourcey, sourcez;
    std::vector<G4double> distance_plate0,distance_plate1;
    G4double energy;
    G4double direction;

    std::vector<G4double> theta, phi;
    G4double esrThickness;
    G4double crystalx;
    G4double crystaly;
    G4int ncrystalx;
    G4int ncrystaly;
    G4int modules,plates;
    G4double crystalz;
    G4double greaseBack;
    G4double glassBack;
    G4double airBack;
    std::string plate_x_s,plate_y_s,plate_z_s,rotation_s;
    std::vector<std::string> plate_x_f,plate_y_f,plate_z_f,rotation_f;
    std::vector<G4double>    plate_x,plate_y,plate_z,rotation;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
