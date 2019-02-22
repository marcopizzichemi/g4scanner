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
// $Id: g4matrixSteppingAction.cc 71007 2013-06-09 16:14:59Z maire $
//
/// \file g4matrixSteppingAction.cc
/// \brief Implementation of the g4matrixSteppingAction class

#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
// #include <sstream>
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "EventAction.hh"
// #include "B4Analysis.hh"
#include "CreateTree.hh"
// #include "struct.hh"

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


SteppingAction::SteppingAction(
  ConfigFile& config)
: G4UserSteppingAction()
//fEventAction(eventAction)
{
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
  fOpticalPhotonsStopped = 0;
  fGammaStopped = 0;
  fElectronStopped = 0;

  quantumEff = config.read<double>("quantumEff");


  //ch[] = {0};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4int eventNumber = G4RunManager::GetRunManager()->
  GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
    G4cout << "Number of Scintillation Photons in previous event: "
    << fScintillationCounter << G4endl;
    G4cout << "Number of Cerenkov Photons in previous event: "
    << fCerenkovCounter << G4endl;
    fEventNumber = eventNumber;
    fScintillationCounter = 0;
    fCerenkovCounter = 0;
  }



  G4Track* track = step->GetTrack();

  // if(track->GetTrackStatus())

  G4String ParticleName = track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  G4String materialName = step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName();

  G4VPhysicalVolume* thePrePV  = step->GetPreStepPoint() ->GetPhysicalVolume();
  G4VPhysicalVolume* thePostPV = step->GetPostStepPoint()->GetPhysicalVolume();

  G4String PreMaterialName ;
  G4String PostMaterialName ;
  G4String PreVolumeName;
  G4String PostVolumeName;

  if(thePrePV)
  {
    PreMaterialName = thePrePV->GetLogicalVolume()->GetMaterial()->GetName();
    PreVolumeName   = thePrePV->GetName();
  }
  if(thePostPV) {
    PostMaterialName = thePostPV->GetLogicalVolume()->GetMaterial()->GetName();
    PostVolumeName   = thePostPV->GetName();
  }


  //getting information from the particles

  //for opticalphoton, we want to know where they stopped
  if (ParticleName == "opticalphoton") //if it's an opticalphoton
  {
    //to avoid exitFound to be propagated to next optical photon, which would happen if in this photon we previously found an exit but then the photon dies, leaving the variable exitFound = true in CreateTree, we set it back to false if the photon is not alive after this step.
    if(track->GetTrackStatus() != 0)
    {
      // std::cout << "DEBUG = " << track->GetTrackStatus() << std::endl;
      CreateTree::Instance()->exitFound = false;
      CreateTree::Instance()->exitFace = 0;
    }

    //BEGIN of test part
    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

    //find the boundary process only once
    if(!boundary)
    {
      G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int i;
      for( i=0;i<nprocesses;i++)
      {
        if((*pv)[i]->GetProcessName()=="OpBoundary")
        {
          boundary = (G4OpBoundaryProcess*)(*pv)[i];
          break;
        }
      }
    }
    boundaryStatus = boundary->GetStatus();

    //each time the optical leaves the crystal, check if it was first time then record if it left from back or front face
    //exit face legend:
    // 0 - no exit was found among the next possibilities
    // 1 - exit from front face (face coupled to MPPC)
    // 2 - exit from back face (face opposite to MPPC)
    if(!CreateTree::Instance()->exitFound)//if this photon never exited from a crystal
    {

      if((boundaryStatus == FresnelRefraction) && (PreMaterialName == "LYSO"))
      {
        // G4cout <<  PostVolumeName << G4endl;
        if( ((PostVolumeName == "GreaseFrontCryToGlass") | (PostVolumeName == "GlassFront") | (PostVolumeName == "GreaseFrontGlassToMPPC") | (PostVolumeName == "AirGapCrystalMPPC")))
        {
          CreateTree::Instance()->exitFace = 1; //FIXME maybe enum? - 0 is front, 1 is back (and 2 lateral?)
          CreateTree::Instance()->exitFound = true; //first exit has been found
        }
        if( ((PostVolumeName == "GreaseBackCryToGlass") | (PostVolumeName == "GlassBack") | (PostVolumeName == "AirBack")))
        {
          CreateTree::Instance()->exitFace = 2; //FIXME maybe enum? - 0 is front, 1 is back (and 2 lateral?)
          CreateTree::Instance()->exitFound = true; //first exit has been found
        }
      }
    }
    if(/*(PreMaterialName == "Epoxy") && */(PostMaterialName == "SilicioMPPC") && (boundaryStatus == FresnelRefraction)) //whatever enters a detector volume
    {
      //G4cout << PreMaterialName << " " << PostMaterialName << G4endl;
      //HACK with names=numbers. better solutions will be welcome...
      G4String detectorName = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
      int numb;
      std::istringstream ( detectorName ) >> numb;

      //take into account quantum efficiency
      G4double rand = G4UniformRand();
      //record the photon, if QE allows it
      if(rand < quantumEff) //absorb the optical photon, save the relevant data
      {

        CreateTree::Instance()->DetectorHit[numb]++;

        G4ThreeVector OnDetectorPosition = track->GetStep()->GetPostStepPoint()->GetPosition(); //get the position vector
        G4ThreeVector PreOnDetectorMomentum = track->GetStep()->GetPreStepPoint()->GetMomentumDirection(); //get the momentum, unit vector
        G4ThreeVector PostOnDetectorMomentum = track->GetStep()->GetPostStepPoint()->GetMomentumDirection(); //get the momentum, unit vector
        G4double globalTime = track->GetGlobalTime();

        // CreateTree::Instance()->PositionX.push_back(OnDetectorPosition.getX());
        // CreateTree::Instance()->PositionY.push_back(OnDetectorPosition.getY());
        // CreateTree::Instance()->PositionZ.push_back(OnDetectorPosition.getZ());
        //
        // CreateTree::Instance()->PreMomentumX.push_back(PreOnDetectorMomentum.getX()/*/CLHEP::nm*/);
        // CreateTree::Instance()->PreMomentumY.push_back(PreOnDetectorMomentum.getY()/*/CLHEP::nm*/);
        // CreateTree::Instance()->PreMomentumZ.push_back(PreOnDetectorMomentum.getZ()/*/CLHEP::nm*/);
        //
        // CreateTree::Instance()->PostMomentumX.push_back(PostOnDetectorMomentum.getX()/*/CLHEP::nm*/);
        // CreateTree::Instance()->PostMomentumY.push_back(PostOnDetectorMomentum.getY()/*/CLHEP::nm*/);
        // CreateTree::Instance()->PostMomentumZ.push_back(PostOnDetectorMomentum.getZ()/*/CLHEP::nm*/);

        // get crystal of origin
        //HACK with names=numbers. better solutions will be welcome...
        // int NoDaughters = track->GetLogicalVolumeAtVertex()->GetNoDaughters();
        // G4cout << "NoDaughters = " << NoDaughters << G4endl;
        std::string CrystalName = (std::string) track->GetLogicalVolumeAtVertex()->GetName();
        // G4cout << "CrystalName = " << CrystalName << G4endl;
        //find crystal i and j
        int iCry,jCry;
        std::size_t foundFirst = CrystalName.find_first_of("_");
        std::size_t foundLast = CrystalName.find_last_of("_");
        std::string iString = CrystalName.substr(foundFirst+1,foundLast-foundFirst-1);
        std::string jString = CrystalName.substr(foundLast+1,CrystalName.size()-foundLast-1);
        // G4cout << "iString = " << iString << " , jString = " << jString << G4endl;
        std::istringstream ( iString ) >> iCry;
        std::istringstream ( jString ) >> jCry;
        // G4cout << "iCry = " << iCry << " , jCry = " << jCry << G4endl;
        // Short_t numbCrystal;
        // std::istringstream ( CrystalName ) >> numbCrystal;
        // CreateTree::Instance()->OriginCrystalI.push_back(iCry);
        // CreateTree::Instance()->OriginCrystalJ.push_back(jCry);
        //BEGIN of DEBUG for negative PostMomentumZ
        // 	if(PostOnDetectorMomentum.getZ() < 0)
        // 	{
        // 	  G4cout << "----------------------------------------------------------------------------------------------------------------" << G4endl;
        // 	  G4cout << "PreMaterialName = " << PreMaterialName << G4endl;
        // 	  G4cout << "PostMaterialName = " << PostMaterialName << G4endl;
        // 	  G4cout << "Position = " << OnDetectorPosition.getX() << " " << OnDetectorPosition.getY() << " " << OnDetectorPosition.getZ() << G4endl;
        // 	  G4cout << "----------------------------------------------------------------------------------------------------------------" << G4endl;
        // 	}
        //END of DEBUG for negative PostMomentumZ

        optPhot tempPhoton;
        //save the process that created the photon
        if(track->GetCreatorProcess()->GetProcessName() == "Scintillation")
        {
          // CreateTree::Instance()->PhotonType.push_back(0);
          tempPhoton.PhotonType = 0;
        }

        else if(track->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        {
          // CreateTree::Instance()->PhotonType.push_back(1);
          tempPhoton.PhotonType = 1;
        }
        else
        {
          // CreateTree::Instance()->PhotonType.push_back(2);
          tempPhoton.PhotonType = 2;
        }
        // CreateTree::Instance()->GlobalTime.push_back(globalTime/CLHEP::ns);
        // CreateTree::Instance()->PhotonEnergy.push_back(track->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV);




        tempPhoton.PositionX = OnDetectorPosition.getX();
        tempPhoton.PositionY = OnDetectorPosition.getY();
        tempPhoton.PositionZ = OnDetectorPosition.getZ();
        tempPhoton.PreMomentumX = PreOnDetectorMomentum.getX();
        tempPhoton.PreMomentumY = PreOnDetectorMomentum.getY();
        tempPhoton.PreMomentumZ = PreOnDetectorMomentum.getZ();
        tempPhoton.PostMomentumX = PostOnDetectorMomentum.getX();
        tempPhoton.PostMomentumY = PostOnDetectorMomentum.getY();
        tempPhoton.PostMomentumZ = PostOnDetectorMomentum.getZ();


        tempPhoton.GlobalTime = globalTime/CLHEP::ns;
        // tempPhoton.PhotonType =
        tempPhoton.PhotonEnergy = track->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV;
        tempPhoton.OriginCrystalI = iCry;
        tempPhoton.OriginCrystalJ = jCry;
        if(CreateTree::Instance()->exitFound)
        {
          tempPhoton.ExitFace = CreateTree::Instance()->exitFace;
        }
        else
        {
          tempPhoton.ExitFace = 0; //back to not fixed
        }
        CreateTree::Instance()->photons.push_back(tempPhoton);


      }
      CreateTree::Instance()->exitFound = false; //we pass to a new track
      CreateTree::Instance()->exitFace = 0;
      //kill the photon
      track->SetTrackStatus(fStopAndKill);
    }
    //END of test part


    //BEGIN of standard method
    //     if(track->GetTrackStatus()==fStopAndKill) //if it just died here
    //     {
    //       if(PostMaterialName == "SilicioMPPC" )
    //       {
    // 	//G4cout << PreMaterialName << " " << PostMaterialName << G4endl;
    // 	G4String detectorName = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    // 	int numb;
    // 	std::istringstream ( detectorName ) >> numb;
    //
    // 	//take into account quantum efficiency
    // 	G4double rand = G4UniformRand();
    //
    // 	if(rand < quantumEff) //absorb the optical photon, save the relevant data
    // 	{
    // 	  CreateTree::Instance()->DetectorHit[numb]++;
    //
    // 	  G4ThreeVector OnDetectorPosition = track->GetStep()->GetPostStepPoint()->GetPosition(); //get the position vector
    // 	  G4ThreeVector OnDetectorMomentum = track->GetStep()->GetPostStepPoint()->GetMomentumDirection(); //get the momentum, unit vector
    // 	  G4double globalTime = track->GetGlobalTime();
    //
    // 	  CreateTree::Instance()->PositionX.push_back(OnDetectorPosition.getX());
    //           CreateTree::Instance()->PositionY.push_back(OnDetectorPosition.getY());
    // 	  CreateTree::Instance()->PositionZ.push_back(OnDetectorPosition.getZ());
    //
    // 	  CreateTree::Instance()->MomentumX.push_back(OnDetectorMomentum.getX()/*/CLHEP::nm*/);
    //           CreateTree::Instance()->MomentumY.push_back(OnDetectorMomentum.getY()/*/CLHEP::nm*/);
    // 	  CreateTree::Instance()->MomentumZ.push_back(OnDetectorMomentum.getZ()/*/CLHEP::nm*/);
    //
    // 	  CreateTree::Instance()->GlobalTime.push_back(globalTime/CLHEP::ns);
    // 	}
    //       }
    //     }
    //END of standard method
  }
  else //there's only gammas and electrons. for them, we want to know where they left energy (opticalphoton don't really leave any)
  {
    if(track->GetTrackID() == 1) // primary gamma
    {
      if(!CreateTree::Instance()->FoundGammaVertexMomentum)
      {
        CreateTree::Instance()->SourceMomentumX = (Float_t) track->GetVertexMomentumDirection().getX();
        CreateTree::Instance()->SourceMomentumY = (Float_t) track->GetVertexMomentumDirection().getY();
        CreateTree::Instance()->SourceMomentumZ = (Float_t) track->GetVertexMomentumDirection().getZ();

        std::ofstream myfile;
        myfile.open ("SourceMomentum.txt",std::ios::app);
        myfile << track->GetVertexMomentumDirection().getX() << " "
               << track->GetVertexMomentumDirection().getY() << " "
               << track->GetVertexMomentumDirection().getZ() << std::endl;
        myfile.close();
        CreateTree::Instance()->FoundGammaVertexMomentum = true;
      }
    }

    //DEBUG
    // G4cout << ParticleName << " " << track->GetTrackID() << G4endl;
    if(PreMaterialName == "LYSO" && track->GetTrackID() == 1) // primary gamma reached the matrix
    {
      if(!CreateTree::Instance()->FoundGammaFirstEntrance)
      {
        // G4cout << track->GetStep()->GetPreStepPoint()->GetPosition()->getX() << " " << track->GetStep()->GetPreStepPoint()->GetPosition()->getY() << std::endl;
        std::ofstream myfile;
        myfile.open ("CryEntranceXYZ.txt",std::ios::app);
        myfile << track->GetStep()->GetPreStepPoint()->GetPosition().getX() << " "
               << track->GetStep()->GetPreStepPoint()->GetPosition().getY() << " "
               << track->GetStep()->GetPreStepPoint()->GetPosition().getZ() << std::endl;
        myfile.close();
        CreateTree::Instance()->FoundGammaFirstEntrance = true;
      }
    }
    //---DEBUG
    G4double edep;
    edep = step->GetTotalEnergyDeposit()/CLHEP::MeV; //check if there was an energy deposition

    if(edep != 0) //if there was en dep, save the data
    {

      if(materialName == "LYSO") //if the electron is interacting with the detector, it does a huge number of cherenkov
      {

        //get the deposition process name
        // 	G4OpBoundaryProcessStatus boundaryStatus=Undefined;
        // 	static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

        //find the boundary process only once
        // 	if(!boundary)
        // 	{
        // 	  G4ProcessManager* pm = track->GetDefinition()->GetProcessManager();
        // 	  G4int nprocesses = pm->GetProcessListLength();
        // 	  G4ProcessVector* pv = pm->GetProcessList();
        // 	  G4int i;

        // 	  for( i=0;i<nprocesses;i++)
        // 	  {
        // 	    if((*pv)[i]->GetProcessName()=="OpBoundary")
        // 	    G4cout << (*pv)[i]->GetProcessType() << " " << (*pv)[i]->GetProcessName() << G4endl;
        // 	    {
        // 	      boundary = (G4OpBoundaryProcess*)(*pv)[i];
        // 	      break;
        // 	    }
        // 	  }
        // 	}
        // 	boundaryStatus = boundary->GetStatus();


        //add total energy deposited
        CreateTree::Instance()->totalEnergyDeposited += edep;
        //take crystal name from physical volume
        G4String crystalNamePV = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
        //convert it to int
        int numb;
        std::istringstream ( crystalNamePV ) >> numb;

        // i and j from logical volume (like for opticals)
        std::string CrystalName = track->GetStep()->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
        // G4cout << "CrystalName En Dep = " << CrystalName << G4endl;
        //find crystal i and j
        int iCry,jCry;
        //now also plate num is in the PV name. Like "Crystal_0_1_2" i.e. plate = 0, i = 1, j = 2
        //so find the first _

        std::size_t foundPre = CrystalName.find_first_of("_");
        std::string subStr1 =  CrystalName.substr(foundPre+1,CrystalName.size()-foundPre-1); //from just after the first "_"

        std::size_t foundFirst = subStr1.find_first_of("_");
        std::size_t foundLast = subStr1.find_last_of("_");

        std::string iString = subStr1.substr(foundFirst+1,foundLast-foundFirst-1);
        std::string jString = subStr1.substr(foundLast+1,subStr1.size()-foundLast-1);
//         G4cout << "EN DEP: iString = " << iString << " , jString = " << jString << G4endl;
        std::istringstream ( iString ) >> iCry;
        std::istringstream ( jString ) >> jCry;
//         G4cout << "EN DEP: iCry = " << iCry << " , jCry = " << jCry << G4endl;

        //add energy deposited in this step to this crystal
        // CreateTree::Instance()->CryEnergyDeposited[numb].push_back(edep);
        //add position of deposited energy
        G4ThreeVector positionVector = track->GetStep()->GetPreStepPoint()->GetPosition(); //get the position vector
        //call the methods to fill the position vectors
        // CreateTree::Instance()->PosXEnDep[numb].push_back(positionVector.getX());
        // CreateTree::Instance()->PosYEnDep[numb].push_back(positionVector.getY());
        // CreateTree::Instance()->PosZEnDep[numb].push_back(positionVector.getZ());

        G4double globalTime = track->GetGlobalTime();


        if(!CreateTree::Instance()->FoundGammaFirstDeposition)
        {
          // G4cout << track->GetStep()->GetPreStepPoint()->GetPosition()->getX() << " " << track->GetStep()->GetPreStepPoint()->GetPosition()->getY() << std::endl;
          std::ofstream myfile;
          myfile.open ("EnDepositionXYZ.txt",std::ios::app);
          myfile << track->GetStep()->GetPreStepPoint()->GetPosition().getX() << " "
                 << track->GetStep()->GetPreStepPoint()->GetPosition().getY() << " "
                 << track->GetStep()->GetPreStepPoint()->GetPosition().getZ() << " "
                 << iCry << " "
                 << jCry << std::endl;

          myfile.close();
          CreateTree::Instance()->FoundGammaFirstDeposition = true;
        }

        enDep energyDeposition;

        energyDeposition.CrystalI        = iCry;
        energyDeposition.CrystalJ        = jCry;
        energyDeposition.CrystalID       = numb;
        energyDeposition.EnergyDeposited = (Float_t) edep;
        energyDeposition.DepositionTime  = (Float_t) globalTime;
        energyDeposition.DepositionX     = (Float_t) positionVector.getX();
        energyDeposition.DepositionY     = (Float_t) positionVector.getY();
        energyDeposition.DepositionZ     = (Float_t) positionVector.getZ();

        CreateTree::Instance()->energyDeposition.push_back(energyDeposition);
      }
      else //if there is energy deposition in silicium, we want to discard the outcome
           //these spikes are seen in real detectors but anyway they are useless and take a lot of CPU time, and space
      {
        track->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }

  if (ParticleName == "opticalphoton") return;
  const std::vector<const G4Track*>* secondaries =
  step->GetSecondaryInCurrentStep();
  if (secondaries->size()>0) {
    for(unsigned int i=0; i<secondaries->size(); ++i) {
      if (secondaries->at(i)->GetParentID()>0) {
        if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
          == G4OpticalPhoton::OpticalPhotonDefinition()){
          if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
            == "Scintillation")
          {
            CreateTree::Instance()->NumOptPhotons++;
            fScintillationCounter++;
          }
          if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
            == "Cerenkov")
          {
            fCerenkovCounter++;
            CreateTree::Instance()->NumCherenkovPhotons++;
          }
          }
      }
    }
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
