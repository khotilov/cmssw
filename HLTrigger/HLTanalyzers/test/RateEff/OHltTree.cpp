#define OHltTree_cxx

#include "OHltTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TFormula.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

void OHltTree::Loop(OHltRateCounter *rc,OHltConfig *cfg,OHltMenu *menu,int procID,float &Den,TH1F* &h1,TH1F* &h2,TH1F* &h3,TH1F* &h4) 
{
  cout<<"Start looping on sample "<<procID<<endl;
  if (fChain == 0) {cerr<<"Error: no tree!"<<endl; return;}

  Long64_t nentries = (Long64_t)cfg->nEntries; 
  if (cfg->nEntries <= 0)
    nentries = fChain->GetEntries();
  cout<<"Entries to be processed: "<<nentries<<endl;

  // Only for experts:
  // Select certain branches to speed up code.
  // Modify only if you know what you do!
  if (cfg->doSelectBranches) {
    fChain->SetBranchStatus("*",kFALSE);
    fChain->SetBranchStatus("MCmu*",kTRUE); // for ppMuX
    fChain->SetBranchStatus("MCel*",kTRUE); // for ppEleX
    if (cfg->selectBranchL1) {
      fChain->SetBranchStatus("L1_*",kTRUE);
    }
    if (cfg->selectBranchHLT) {
      fChain->SetBranchStatus("HLT_*",kTRUE);
      fChain->SetBranchStatus("AlCa_*",kTRUE);
    }
    if (cfg->selectBranchOpenHLT) {
      fChain->SetBranchStatus("*oh*",kTRUE);
      fChain->SetBranchStatus("*recoJet*",kTRUE);
      fChain->SetBranchStatus("*recoMet*",kTRUE);
    }
    if (cfg->selectBranchReco) {
      fChain->SetBranchStatus("*reco*",kTRUE);
    }
    if (cfg->selectBranchL1extra) {
      fChain->SetBranchStatus("*L1*",kTRUE);
    }
    if (cfg->selectBranchMC) {
      fChain->SetBranchStatus("*MC*",kTRUE);
    }
  } else {
    fChain->SetBranchStatus("*",kTRUE);
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    fChain->GetEntry(jentry);

    if (jentry%cfg->nPrintStatusEvery == 0)
      cout<<"Processing entry "<<jentry<<"/"<<nentries<<"\r"<<flush<<endl;

    if ( cfg->pdomucuts[procID] && MCmu3!=0 ) continue;
    if ( cfg->pdoecuts[procID] && MCel3!=0 ) continue;
 
    SetOpenL1Bits(); 

    RemoveEGOverlaps();
   
    // 1. Loop to check which Bit fired
    // Triggernames are assigned to trigger cuts in unambigous way!
    // If you define a new trigger also define a new unambigous name!
    if(menu->DoL1preLoop() && menu->IsHltMenu()) {
      ApplyL1Prescales(menu);
    }

    //SetMapL1BitOfStandardHLTPath(menu);
    SetMapL1BitOfStandardHLTPathUsingLogicParser(menu,(int)jentry);
    SetL1MuonQuality();
	  
    //////////////////////////////////////////////////////////////////
    // Get Denominator (normalization and acc. definition) for efficiency evaluation
    //////////////////////////////////////////////////////////////////


    if (cfg->pnames[procID]=="zee"||cfg->pnames[procID]=="zmumu"){
      int accMCMu=0;
      int accMCEle=0;
      for(int iMCpart = 0; iMCpart < NMCpart; iMCpart ++){
	if((MCpid[iMCpart]==13||MCpid[iMCpart]==-13) && MCstatus[iMCpart]==1 && (MCeta[iMCpart] < 2.1 && MCeta[iMCpart] > -2.1) && (MCpt[iMCpart]>3))accMCMu=accMCMu+1;
	if((MCpid[iMCpart]==11||MCpid[iMCpart]==-11 )&& MCstatus[iMCpart]==1 && (MCeta[iMCpart] < 2.5 && MCeta[iMCpart] > -2.5) && (MCpt[iMCpart]>5))accMCEle=accMCEle+1;
      }
      if((cfg->pnames[procID]=="zee" && accMCEle>1)||(cfg->pnames[procID] == "zmumu" && accMCMu>1)){ Den=Den+1;}
      else{continue;}
    }


    //////////////////////////////////////////////////////////////////
    // Make efficiency curves
    //////////////////////////////////////////////////////////////////
    TString hlteffmode;
    TString ohltobject;
    hlteffmode="GEN";
    //    hlteffmode="L1";
    //    hlteffmode="RECO";
    ohltobject="None";
    if (cfg->pnames[procID]=="zee")ohltobject="electron";
    if (cfg->pnames[procID]=="zmumu")ohltobject="muon";
    PlotOHltEffCurves(hlteffmode,ohltobject,h1,h2,h3,h4);


    //////////////////////////////////////////////////////////////////
    // Loop over trigger paths and do rate counting
    //////////////////////////////////////////////////////////////////
    for (int i = 0; i < nTrig; i++){
      triggerBit[i] = false;
      previousBitsFired[i] = false;
      allOtherBitsFired[i] = false;

      //////////////////////////////////////////////////////////////////
      // Standard paths
      TString st = menu->GetTriggerName(i);
      if (st.BeginsWith("HLT_") || st.BeginsWith("L1_") || st.BeginsWith("AlCa_")) {
	// Prefixes reserved for Standard HLT&L1	
	if ( (map_BitOfStandardHLTPath.find(st)->second==1) ) {	
	  if (map_L1BitOfStandardHLTPath.find(st)->second>0) {
	    if (GetIntRandom() % menu->GetPrescale(i) == 0) { triggerBit[i] = true; }
	  }
	}
      } else {
	CheckOpenHlt(cfg,menu,i);
      }
    }

    /* ******************************** */
    // 2. Loop to check overlaps
    for (int it = 0; it < nTrig; it++){
      if (triggerBit[it]) {
	rc->iCount[it]++;
	for (int it2 = 0; it2 < nTrig; it2++){
	  if (triggerBit[it2]) {
	    rc->overlapCount[it][it2] += 1;
	    if (it2<it)
	      previousBitsFired[it] = true;
	    if (it2!=it)
	      allOtherBitsFired[it] = true;
	  }
	}
	if (not previousBitsFired[it])
	  rc->sPureCount[it]++;
	if (not allOtherBitsFired[it])
	  rc->pureCount[it]++;
      }
    }
    /* ******************************** */

    
  }
  
}

void OHltTree::SetLogicParser(std::string l1SeedsLogicalExpression) {

  //if (l1SeedsLogicalExpression != "") {
    
    //std::cout<<"@@@ L1 condition: "<<l1SeedsLogicalExpression<<std::endl;
    // check also the logical expression - add/remove spaces if needed
    m_l1AlgoLogicParser.push_back(new L1GtLogicParser(l1SeedsLogicalExpression));
    
  //}
};

