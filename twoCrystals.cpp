// compile with (from directory build)
// g++ -o ../build/twoCrystals twoCrystals.cpp `root-config --cflags --glibs`
// syntax
// twoCrystals input1.root [input2.root ...]

#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TF2.h"
#include "TH3I.h"


#include "../code/include/struct.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>    // std::sort
#include <numeric>      // std::accumulate
#include <iomanip>      // std::setprecision
#include <vector>

//struct of averaged crystal events
struct avgCryEnergyDep
{
  int id;
  float x;
  float y;
  float z;
  float sx;
  float sy;
  float sz;
  float energy;
  float time;
};

//struct of crystal related data
struct CrystalData
{
  int id;
  double x;
  double y;
  int mppci;
  int mppcj;
  TGraph* wz; //w(z) graph for this crystal
  TGraphDelaunay*** gd; //pointers to the pi_(w,E) for this crystal
  TF2*** quadratic; //quadratic function describing the pi_(z,p_tot), fitted by matlab on the calibration data
  Float_t correction; //
};

//function to compare deposition event struct vectors using the field time
bool compareByTime(const enDep &a,const enDep  &b)
{
  return a.DepositionTime < b.DepositionTime;
}


//MAIN
int main (int argc, char** argv)
{
  if(argc < 3) //check input, provide usage explainations
  {
    std::cout << "USAGE:\t\t\t twoCrystals  inputFile length" << std::endl;
    std::cout << "inputFile \t\t input root file name (output of simulation)" << std::endl;
    std::cout << "length    \t\t length of scintillators [mm]" << std::endl;
    std::cout << std::endl;
    return 1;
  }


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        ROOT STUFF                                    //
  //                                                                      //
  //----------------------------------------------------------------------//
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors
  //HACK to use the dictionary easily
  // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
  std::string fullFileName = "";
  std::string path = "";
  pid_t pid = getpid();
  char buf[20] = {0};
  sprintf(buf,"%d",pid);
  std::string _link = "/proc/";
  _link.append( buf );
  _link.append( "/exe");
  char proc[512];
  int ch = readlink(_link.c_str(),proc,512);
  if (ch != -1) {
    proc[ch] = 0;
    path = proc;
    std::string::size_type t = path.find_last_of("/");
    path = path.substr(0,t);
  }
  fullFileName = path + std::string("/");
  std::string command = ".L " + fullFileName + "structDictionary.C+";
  gROOT->ProcessLine(command.c_str());
  // gStyle->SetOptStat(0);
  //----------------------------------------------------------------------//
  //                                                                      //
  //                        INPUT FILES                                   //
  //                                                                      //
  //----------------------------------------------------------------------//
  TChain *tree =  new TChain("tree"); // read input files
  float length = atof(argv[1]);
  for (int i = 2 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[1] << std::endl;
    tree->Add(argv[2]);
  }

  // find the number of channels directly from the tchain file
  // before creating the variables
  // first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++){
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  std::string det_prefix("detector");
  for(int i = 0 ; i < nLeaves ; i++){
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix)) numOfCh++;
  }
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  //FIXME hardcoded for now
  int crystalsPerArray = 1;

  float backPosition = 150.0 + (length/2.0); // mm
  float lightSpeed = 3e11; // mm/s
  float refIndex = 1.8;


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        INPUT TTREES                                  //
  //                                                                      //
  //----------------------------------------------------------------------//
  //create the branches in the input ttree and connect to the variables
  // ----Global variables
  // these are 1 number per TTree entry - so 1 number per gamma shot
  Long64_t Seed;                      // seed of the simulation (read every time, but always the same)
  int Run;                            // run id (usually just 1)(read every time, but always the same)
  int Event;                          // event id
  float totalEnergyDeposited;         // total energy deposited in this event, in all the matrix
  Float_t sourcex;
  Float_t sourcey;
  Float_t sourcez;
  Float_t sourceMomentumX;
  Float_t sourceMomentumY;
  Float_t sourceMomentumZ;
  int NumOptPhotons;                  // number of optical photons generated in this event, in the entire matrix
  int NumCherenkovPhotons;            // number of Cherenkov photons generated in this event, in the entire matrix
  // energy deposition, each gamma 511 event has a std::vector of struct (type enDep) with all the data of each energy deposition
  std::vector<enDep> *energyDeposition = 0;
  // Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];
  // optical photons. for each gamma 511 event, every optical photon detected is a struct of type optPhot. a std::vector<optPhot> is saved for each gamma 511
  std::vector<optPhot> *photons = 0;
  // ----Brach addresses
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("SourceX",&sourcex);
  tree->SetBranchAddress("SourceY",&sourcey);
  tree->SetBranchAddress("SourceZ",&sourcez);
  tree->SetBranchAddress("SourceMomentumX",&sourceMomentumX);
  tree->SetBranchAddress("SourceMomentumY",&sourceMomentumY);
  tree->SetBranchAddress("SourceMomentumZ",&sourceMomentumZ);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  tree->SetBranchAddress("optical",&photons);
  tree->SetBranchAddress("energyDeposition",&energyDeposition);
  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }







  //----------------------------------------------------------------------//
  //                                                                      //
  //                        OUTPUT HISTOGRAMS                             //
  //                                                                      //
  //----------------------------------------------------------------------//
  long int totalDepositionCounter = 0;
  long int singleCounter = 0;
  long int doubleCounter = 0;
  long int tripleCounter = 0;
  long int quadrupleCounter = 0;
  long int multipleCounter = 0;
  long int globalReturned = 0;
  long int HighestCrystalRight = 0;
  long int HighestCrystalWrong = 0;



  //----------------------------------------------------------------------//
  //                                                                      //
  //                        LOOP ON EVENTS                                //
  //                                                                      //
  //----------------------------------------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  long int goodCounter = 0;
  long int goodChi = 0;
  long int goodCounterOnlyCompton = 0;
  long int goodCounterNoPhotEl = 0;
  std::cout << "nEntries = " << nEntries << std::endl;
  // int cEvent = 5197;
  // int cEvent = 400;
  // int cEvent = 398;

  // int cEvent = 3620;
  // start and stop events
  int startEvent = 0;
  int stopEvent = nEntries;
  // if(!allAnalysis)
  // {
  //   startEvent = cEvent;
  //   stopEvent = cEvent+1;
  // }

  TH1F *histo = new TH1F("histo","",200,-100,100);
  histo->GetXaxis()->SetTitle("Delay [ps]");
  histo->GetYaxis()->SetTitle("Normalized Counts");

  for(int iEvent = startEvent; iEvent < stopEvent ; iEvent++)
  // for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  // for(int iEvent = cEvent; iEvent < cEvent + 1 ; iEvent++)
  {
    tree->GetEvent(iEvent);

    // "average" what happened into each crystal
    // stackingaction in g4 is last-in-first-out so we need to sort energyDeposition before we check what crystal was hit first
    // check crystal sequence, check for returning gammas, calculate an average xyz per crystal
    std::sort(energyDeposition->begin(), energyDeposition->end(), compareByTime);
    std::vector<std::vector < enDep > > separatedEnDep;

    int CurrentID = -1;
    float RealX = 0.0;
    float RealY = 0.0;
    float RealZ = 0.0;
    std::vector < enDep > CrystalEnDepCollection;
    for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
    {
      //this for cycles on all the endep events. it stores them (each one is a struct) in a std::vector of structs, until the gamma changes crystal
      //when we enter a new crystal, the std::vector of that crystal is pushed_back into (guess what?) a std::vector of std::vector and another std::vector is created
      // so to sumarize:
      // energy deposition event -> struct
      // all energy deposition events in a given crystal -> std::vector<struct> -> i call this a "collection" below
      // all crystals touched by the gamma -> std::vector< std:: vector <struct> > this would be the std::vector of "collections"
      //read the crystal where energy was deposited
      int cry = energyDeposition->at(eEvent).CrystalID;
      if(eEvent == 0) CurrentID = cry; //needed for the first step
      //create temp enDep variable and copy this eEvent into it
      enDep tempCrystalEnDep;
      tempCrystalEnDep = energyDeposition->at(eEvent);
      if(cry != CurrentID) // if this endep event is happening in a new crystal wrt the one before
      {
        separatedEnDep.push_back(CrystalEnDepCollection); // save the collection of this crystal into the std::vector of collections
        CrystalEnDepCollection.clear(); //clear this collection
        CurrentID = cry; //change the current id
      }

      CrystalEnDepCollection.push_back(tempCrystalEnDep); // save this enDep event into the collection of this crystal
      if(eEvent == energyDeposition->size() -1)
      {
        separatedEnDep.push_back(CrystalEnDepCollection);
      }
    }

    // now for each crystal average what happened inside
    // what happened in a crystal cannot be seen in split details by the detector so we need some useful
    // variables to gain the chance to filter the simulation dataset, if we want to understand
    // what's going on
    std::vector<avgCryEnergyDep> averageDepEvents;
    // now the en dep events are collected by crystal
    // run on each collection and find average
    for(int iColl = 0 ; iColl < separatedEnDep.size(); iColl++)
    {
      avgCryEnergyDep tempAvgEnDep;

      //initialize
      tempAvgEnDep.id = separatedEnDep.at(iColl).at(0).CrystalID;
      tempAvgEnDep.x = 0;
      tempAvgEnDep.y = 0;
      tempAvgEnDep.z = 0;
      tempAvgEnDep.time = 0;
      tempAvgEnDep.energy = 0;
      tempAvgEnDep.sx = 0;
      tempAvgEnDep.sy = 0;
      tempAvgEnDep.sz = 0;
      //now run on the energy deposition events in this crystal and calculate the averages
      for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
      {
        tempAvgEnDep.energy += separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.x += separatedEnDep.at(iColl).at(iEndep).DepositionX * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.y += separatedEnDep.at(iColl).at(iEndep).DepositionY * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.z += separatedEnDep.at(iColl).at(iEndep).DepositionZ * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.time += separatedEnDep.at(iColl).at(iEndep).DepositionTime * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
      }
      tempAvgEnDep.x = tempAvgEnDep.x / tempAvgEnDep.energy;
      tempAvgEnDep.y = tempAvgEnDep.y / tempAvgEnDep.energy;
      tempAvgEnDep.z = tempAvgEnDep.z / tempAvgEnDep.energy;
      tempAvgEnDep.time = tempAvgEnDep.time / tempAvgEnDep.energy;
      float varx = 0.0; // variance (needed for stdev afterwards)
      float vary = 0.0;
      float varz = 0.0;
      for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
      {
        varx += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionX  - tempAvgEnDep.x,2)) / tempAvgEnDep.energy;
        vary += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionY  - tempAvgEnDep.y,2)) / tempAvgEnDep.energy;
        varz += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionZ  - tempAvgEnDep.z,2)) / tempAvgEnDep.energy;
      }
      tempAvgEnDep.sx = sqrt(varx);
      tempAvgEnDep.sy = sqrt(vary);
      tempAvgEnDep.sz = sqrt(varz);

      //detection threshold
      if(tempAvgEnDep.energy > 0.05)
      {  //save into the std::vector of averages
        averageDepEvents.push_back(tempAvgEnDep);
      }
    }

    //two crystals or more hit, but a crystal is hit more than once (like, 28 -> 36 -> 28)
    std::vector <int> checkReturning; // this std::vector has same size of averageDepEvents if there is no returning, smaller if there is
    int returned = 0;
    bool isReturned = false;
    for (int iAvg = 0 ; iAvg < averageDepEvents.size() ; iAvg++)
    {
      int cry = averageDepEvents.at(iAvg).id;
      bool sameID = false;
      for(int i = 0 ; i < checkReturning.size(); i++)
      {
        if(cry == checkReturning[i])
        {
          sameID = true;
          isReturned = true;
          returned++;
        }
      }
      if(!sameID) checkReturning.push_back(cry);
    }
    globalReturned += returned;

    // now let's start filtering the dataset.
    // we start by default
    // bool isCandidate = true;
    if(totalEnergyDeposited < 1.000) //filter out events with non total depostion in two crystals
    {
      // isCandidate = false; // cut on totalEnergyDeposited in this event
    }
    else
    {
      totalDepositionCounter++;

      if(averageDepEvents.size() == 1) singleCounter++; // only one crystal hit, but impossible if totalEnergyDeposited > 1 MeV (we shoot only 2 x 511KeV back to back..)
      if(averageDepEvents.size() == 2)  // two crystals hit in total, it has to be 511 in one and 511 in the other, since totalEnergyDeposited > 1 MeV
      {
        //just check that the 2 are in two different arrays... as it should be
        int array0 = averageDepEvents.at(0).id / crystalsPerArray;
        int array1 = averageDepEvents.at(1).id / crystalsPerArray;


        if(array0 == array1 ) // check if the 2 depo are in the same array (should not be possible)
        {
          std::cout << "same array for both 511 KeV depositions???" << std::endl;
        }
        else // do what you want with these events...
        {
          float d0,d1;
          float z0,z1;
          //swap averageDepEvents if they are in reversed order
          if(array0 == 0 && array1 == 1)
          {
            // ok
            d0 = sqrt(pow(averageDepEvents.at(0).x,2)+pow(averageDepEvents.at(0).y,2)+pow(averageDepEvents.at(0).z,2));
            d1 = sqrt(pow(averageDepEvents.at(1).x,2)+pow(averageDepEvents.at(1).y,2)+pow(averageDepEvents.at(1).z,2));
            z0 = averageDepEvents.at(0).z;
            z1 = averageDepEvents.at(1).z;
          }
          else
          {
            // reverse
            d1 = sqrt(pow(averageDepEvents.at(0).x,2)+pow(averageDepEvents.at(0).y,2)+pow(averageDepEvents.at(0).z,2));
            d0 = sqrt(pow(averageDepEvents.at(1).x,2)+pow(averageDepEvents.at(1).y,2)+pow(averageDepEvents.at(1).z,2));
            z1 = averageDepEvents.at(0).z;
            z0 = averageDepEvents.at(1).z;

          }
          doubleCounter++;


          //gamma 0
          float gamma0dist = fabs(d0);
          float gamma0t = gamma0dist / (lightSpeed);
          float opt0dist = backPosition - fabs(z0);
          float opt0t = opt0dist / (lightSpeed/refIndex);
          float t0 = gamma0t + opt0t;

          //gamma 1
          float gamma1dist = fabs(d1);
          float gamma1t = gamma1dist / (lightSpeed);
          float opt1dist = backPosition - fabs(z1);
          float opt1t = opt1dist / (lightSpeed/refIndex);
          float t1 = gamma1t + opt1t;

          histo->Fill((t0-t1)*1e12); // in ps
        }





      }
      if(averageDepEvents.size() > 2)
      {
        int hit[2] = {0,0};
        std::vector<float> hitTime[2];
        std::vector<float> hitEnergy[2];
        for(int iAverage = 0 ; iAverage < averageDepEvents.size(); iAverage++)
        {
          hit[(averageDepEvents.at(iAverage).id / crystalsPerArray)]++;
          hitTime[(averageDepEvents.at(iAverage).id / crystalsPerArray)].push_back(averageDepEvents.at(iAverage).time);
          hitEnergy[(averageDepEvents.at(iAverage).id / crystalsPerArray)].push_back(averageDepEvents.at(iAverage).energy);
        }

        if( ((hit[0] == 1) && (hit[1] == 2)) || ((hit[0] == 2) && (hit[1] == 1))  )
        {
          tripleCounter++;
        }

        if( ((hit[0] == 1) && (hit[1] == 2)) || ((hit[0] == 2) && (hit[1] == 1)) || ( (hit[0] == 2) && (hit[1] == 2) ) )
        {
          for(int iHit = 0; iHit < 2 ;iHit++)
          {
            if(hit[iHit] == 2)
            {
              if(hitTime[iHit].at(0) > hitTime[iHit].at(1))
              {
                if(hitEnergy[iHit].at(0) > hitEnergy[iHit].at(1) )
                  HighestCrystalWrong++;
                else
                  HighestCrystalRight++;
              }
              else
              {
                if(hitEnergy[iHit].at(0) > hitEnergy[iHit].at(1) )
                  HighestCrystalRight++;
                else
                  HighestCrystalWrong++;
              }
            }
          }

        }




        if(hit[0] == 2)
        {
          if(hit[1] == 2)
          {
            quadrupleCounter++;
          }
        }

        if( (hit[0] > 2) | (hit[1] >2))
        {
          multipleCounter++;
        }



      }
      // if(averageDepEvents.size() > 3)  multipleCounter++;
    }



    counter++;

    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      // std::cout << counter << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "Total events with > 1 MeV deposition = "    << totalDepositionCounter << std::endl;
  std::cout << "1 cry [511 KeV deposition] events = "    << singleCounter << std::endl;
  std::cout << "1+1 511 KeV deposition events = "    << doubleCounter << std::endl;
  std::cout << "1+2 or 2+1 511 KeV deposition events = "    << tripleCounter << std::endl;
  std::cout << "2+2 511 KeV deposition events = "<< quadrupleCounter << std::endl;
  std::cout << "multiple 511 KeV deposition events = "<< multipleCounter << std::endl;
  std::cout << "Returned = " << globalReturned << std::endl;
  std::cout << "-------------------------------------------------------" << std::endl;
  std::cout << "Highest Energy Right = "<< HighestCrystalRight << std::endl;
  std::cout << "Highest Energy Wrong = "<< HighestCrystalWrong << std::endl;

  TCanvas *c_histo = new TCanvas("c_histo","c_histo",1200,800);
  c_histo->SetGrid();

  histo->GetYaxis()->SetNdivisions(111);
  histo->Scale(1.0/histo->GetMaximum());
  histo->Draw("HIST");
  std::stringstream sstream;
  sstream << "plots_" << argv[1] << ".root";
  std::string outputFileName = sstream.str();

  TFile* fOut = new TFile(outputFileName.c_str(),"recreate");
  c_histo->Write();
  histo->Write();
  fOut->Close();



  return 0;
}
