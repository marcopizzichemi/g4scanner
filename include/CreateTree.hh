#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "struct.hh"
#include "globals.hh"

class CreateTree
{
private:

  TTree*              ftree;
  TString             fname;

  Bool_t              HITS;
  Bool_t              ABSORPTIONS;

  static const Int_t  MaxNum = 2000000;
  static const Int_t  MaxNumPro = 100;

  int nDetectorsX, nDetectorsY, nPlates;

public:

  CreateTree(TString name,int x, int y, int z);
  ~CreateTree();


  TTree*              GetTree() const { return ftree; };
  TString             GetName() const { return fname;};
  Int_t               Fill() { return this->GetTree()->Fill(); };
  Bool_t              Write();
  void                Clear();
  static CreateTree*  Instance() { return fInstance; };
  static CreateTree*  fInstance;

  Bool_t              Hits() const { return this->HITS; };
  Bool_t              Absorptions() const { return this->ABSORPTIONS; };

  //Short_t             tempOriginCrystal;


  // Global part
  long int            Seed;
  Int_t               Run;
  Int_t               Event;
  Float_t             SourceX;
  Float_t             SourceY;
  Float_t             SourceZ;
  Float_t             SourceMomentumX;
  Float_t             SourceMomentumY;
  Float_t             SourceMomentumZ;
  Float_t             totalEnergyDeposited;
  Int_t               NumOptPhotons;
  Int_t               NumCherenkovPhotons;

  //Hits on detectors
  Short_t*            DetectorHit;

  // Optical Photons part
  std::vector<optPhot>    photons;

  //energy deposition part
  std::vector<enDep> energyDeposition;

  bool exitFound;
  int exitFace;
  bool FoundGammaFirstEntrance;
  bool FoundGammaFirstDeposition;
  bool FoundGammaVertexMomentum;

  //exit face legend:
  // 0 - no exit was found among the next possibilities
  // 1 - exit from front face (face coupled to MPPC)
  // 2 - exit from back face (face opposite to MPPC)


  // std::vector<float>*  CryEnergyDeposited;
  // std::vector<float>** pCryEnergyDeposited;
  //
  // std::vector<float>*  PosXEnDep;
  // std::vector<float>** pPosXEnDep;
  // std::vector<float>*  PosYEnDep;
  // std::vector<float>** pPosYEnDep;
  // std::vector<float>*  PosZEnDep;
  // std::vector<float>** pPosZEnDep;

  // std::vector<float> PositionX;
  // std::vector<float> *pPositionX;
  // std::vector<float> PositionY;
  // std::vector<float> *pPositionY;
  // std::vector<float> PositionZ;
  // std::vector<float> *pPositionZ;
  //
  // std::vector<float> PreMomentumX;
  // std::vector<float> *pPreMomentumX;
  // std::vector<float> PreMomentumY;
  // std::vector<float> *pPreMomentumY;
  // std::vector<float> PreMomentumZ;
  // std::vector<float> *pPreMomentumZ;
  //
  // std::vector<float> PostMomentumX;
  // std::vector<float> *pPostMomentumX;
  // std::vector<float> PostMomentumY;
  // std::vector<float> *pPostMomentumY;
  // std::vector<float> PostMomentumZ;
  // std::vector<float> *pPostMomentumZ;
  //
  // std::vector<float> GlobalTime;
  // std::vector<float> *pGlobalTime;
  //
  // std::vector<int> PhotonType;
  // std::vector<int> *pPhotonType;
  //
  // std::vector<float> PhotonEnergy;
  // std::vector<float> *pPhotonEnergy;
  //
  // std::vector<int> OriginCrystalI;
  // std::vector<int> *pOriginCrystalI;
  // std::vector<int> OriginCrystalJ;
  // std::vector<int> *pOriginCrystalJ;


};
