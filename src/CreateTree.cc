#include "CreateTree.hh"
#include <sstream>
#include "TROOT.h"
// #include "struct.hh"

CreateTree* CreateTree::fInstance = NULL;

using namespace std;

CreateTree::CreateTree(TString name, int x, int y, int z)
:
nDetectorsX(x),
nDetectorsY(y),
nPlates(z)
{

  //G4cout << "DETECTORS "  << " " << nDetectorsX << " " << nDetectorsY << G4endl; ;
  gROOT->ProcessLine("#include <vector>"); //this is needed otherwise ROOT will complain about not knowing what a std::vector is...


  // gROOT->ProcessLine(".L ./structDictionary.C+");
  //create the vectors
  DetectorHit         = new Short_t [nDetectorsX*nDetectorsY*nPlates];
  // CryEnergyDeposited  = new std::vector<float> [nCrystalsX*nCrystalsY];
  // pCryEnergyDeposited = new std::vector<float>* [nCrystalsX*nCrystalsY];
  // PosXEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  // pPosXEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  // PosYEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  // pPosYEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  // PosZEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  // pPosZEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];


  if(fInstance)
  {
    return;
  }

  this->HITS=true;
  this->ABSORPTIONS=true;

  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree("tree","name");


  this->GetTree()->Branch("Seed",&this->Seed,"Seed/L");
  this->GetTree()->Branch("Run",&this->Run,"Run/I");
  this->GetTree()->Branch("Event",&this->Event,"Event/I");
  this->GetTree()->Branch("SourceX",&this->SourceX,"SourceX/F");
  this->GetTree()->Branch("SourceY",&this->SourceY,"SourceY/F");
  this->GetTree()->Branch("SourceZ",&this->SourceZ,"SourceZ/F");
  this->GetTree()->Branch("SourceMomentumX",&this->SourceMomentumX,"SourceMomentumX/F");
  this->GetTree()->Branch("SourceMomentumY",&this->SourceMomentumY,"SourceMomentumY/F");
  this->GetTree()->Branch("SourceMomentumZ",&this->SourceMomentumZ,"SourceMomentumZ/F");
  this->GetTree()->Branch("totalEnergyDeposited",&this->totalEnergyDeposited,"totalEnergyDeposited/F");
  this->GetTree()->Branch("NumOptPhotons",&this->NumOptPhotons,"NumOptPhotons/I");
  this->GetTree()->Branch("NumCherenkovPhotons",&this->NumCherenkovPhotons,"NumCherenkovPhotons/I");
  this->GetTree()->Branch("energyDeposition","std::vector<enDep>",&energyDeposition);
  this->GetTree()->Branch("optical","std::vector<optPhot>",&photons);
  for (int i = 0 ; i < nDetectorsX*nDetectorsY*nPlates ; i++)
  {
    std::stringstream snames,stypes;
    snames << "detector" << i;
    stypes << "detector" << i << "/S";
    this->GetTree()->Branch(snames.str().c_str(),&DetectorHit[i],stypes.str().c_str());
  }
  exitFound = false;
  exitFace = 0;
  FoundGammaFirstEntrance = false;
  FoundGammaVertexMomentum = false;
  // default to 0
  //exit face legend:
  // 0 - no exit was found among the next possibilities
  // 1 - exit from front face (face coupled to MPPC)
  // 2 - exit from back face (face opposite to MPPC)


  // this->GetTree()->Branch("exitFound",&this->exitFound,"exitFound/B");

  // for(int i = 0; i < nCrystalsX*nCrystalsY ; i++)
  // {
  //   pCryEnergyDeposited[i] = &CryEnergyDeposited[i];
  //   pPosXEnDep[i] = &PosXEnDep[i];
  //   pPosYEnDep[i] = &PosYEnDep[i];
  //   pPosZEnDep[i] = &PosZEnDep[i];
  // }



  // for (int i = 0 ; i < nCrystalsX*nCrystalsY ; i++)
  // {
  //   std::stringstream snames;
  //   snames << "cry" << i;
  //   this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pCryEnergyDeposited[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosXEnDep";
  //   this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosXEnDep[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosYEnDep";
  //   this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosYEnDep[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosZEnDep";
  //   this->GetTree()->Branch(snames.str().c_str(),"std::vector<float>",&pPosZEnDep[i]);
  // }




  // 	pTransmissionX = &TransmissionX;
  // 	pTransmissionY = &TransmissionY;
  // 	pTransmissionZ = &TransmissionZ;
  //
  // 	this->GetTree()->Branch("TransmissionX","std::vector<float>",&pTransmissionX);
  // 	this->GetTree()->Branch("TransmissionY","std::vector<float>",&pTransmissionY);
  // 	this->GetTree()->Branch("TransmissionZ","std::vector<float>",&pTransmissionZ);
  //


  //
  //
  // pPositionX = &PositionX;
  // pPositionY = &PositionY;
  // pPositionZ = &PositionZ;
  //
  // this->GetTree()->Branch("PositionX","std::vector<float>",&pPositionX);
  // this->GetTree()->Branch("PositionY","std::vector<float>",&pPositionY);
  // this->GetTree()->Branch("PositionZ","std::vector<float>",&pPositionZ);
  //
  // pPreMomentumX = &PreMomentumX;
  // pPreMomentumY = &PreMomentumY;
  // pPreMomentumZ = &PreMomentumZ;
  //
  // this->GetTree()->Branch("PreMomentumX","std::vector<float>",&pPreMomentumX);
  // this->GetTree()->Branch("PreMomentumY","std::vector<float>",&pPreMomentumY);
  // this->GetTree()->Branch("PreMomentumZ","std::vector<float>",&pPreMomentumZ);
  //
  // pPostMomentumX = &PostMomentumX;
  // pPostMomentumY = &PostMomentumY;
  // pPostMomentumZ = &PostMomentumZ;
  //
  // this->GetTree()->Branch("PostMomentumX","std::vector<float>",&pPostMomentumX);
  // this->GetTree()->Branch("PostMomentumY","std::vector<float>",&pPostMomentumY);
  // this->GetTree()->Branch("PostMomentumZ","std::vector<float>",&pPostMomentumZ);
  //
  // pGlobalTime = &GlobalTime;
  // this->GetTree()->Branch("GlobalTime","std::vector<float>",&pGlobalTime);
  //
  // pPhotonType = &PhotonType;
  // this->GetTree()->Branch("PhotonType","std::vector<int>",&pPhotonType);
  //
  // pPhotonEnergy = &PhotonEnergy;
  // this->GetTree()->Branch("PhotonEnergy","std::vector<float>",&pPhotonEnergy);
  //
  // pOriginCrystalI = &OriginCrystalI;
  // pOriginCrystalJ = &OriginCrystalJ;
  // this->GetTree()->Branch("OriginCrystalI","std::vector<int>",&pOriginCrystalI);
  // this->GetTree()->Branch("OriginCrystalJ","std::vector<int>",&pOriginCrystalJ);

  this->Clear();
}

CreateTree::~CreateTree()
{
  delete DetectorHit;
  // delete CryEnergyDeposited;
  // delete pCryEnergyDeposited;
  // delete PosXEnDep ;
  // delete pPosXEnDep;
  // delete PosYEnDep;
  // delete pPosYEnDep;
  // delete PosZEnDep;
  // delete pPosZEnDep;
  // delete photons;
}

Bool_t CreateTree::Write()
{
  TString filename = this->GetName();
  filename+=".root";
  TFile* file = new TFile(filename,"RECREATE");
  this->GetTree()->Write();
  file->Write();
  file->Close();
  return true;
}

void CreateTree::Clear()
{
  Run=0;
  Event=0;
  SourceX=0;
  SourceY=0;
  SourceZ=0;
  SourceMomentumX=0;
  SourceMomentumY=0;
  SourceMomentumZ=0;

  NumOptPhotons=0;
  NumCherenkovPhotons=0;
  totalEnergyDeposited=0;

  // for(int i = 0; i < nCrystalsX*nCrystalsY ; i++)
  // {
  //   CryEnergyDeposited[i].clear();
  //   PosXEnDep[i].clear();
  //   PosYEnDep[i].clear();
  //   PosZEnDep[i].clear();
  // }

  for (int i = 0 ; i < nDetectorsX*nDetectorsY*nPlates ; i++)//
  {
    DetectorHit[i] = 0;
  }

  // 	TransmissionX.clear();
  // 	TransmissionY.clear();
  // 	TransmissionZ.clear();

  // PositionX.clear();
  // PositionY.clear();
  // PositionZ.clear();
  //
  // PreMomentumX.clear();
  // PreMomentumY.clear();
  // PreMomentumZ.clear();
  //
  // PostMomentumX.clear();
  // PostMomentumY.clear();
  // PostMomentumZ.clear();
  //
  // GlobalTime.clear();
  // PhotonType.clear();
  // PhotonEnergy.clear();
  // OriginCrystalI.clear();
  // OriginCrystalJ.clear();
  exitFound = false;
  exitFace = 0;
  photons.clear();
  energyDeposition.clear();
}
