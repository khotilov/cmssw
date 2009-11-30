#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <string>
#include <vector>
#include <iostream>

// important: define root tree variable list outside of main function to make it accessible from subfunctions.
#include "VariableList.h"


// declare some functions. These are implemented at the end of this file
void increaseCountersByOne(int &BsJpsiPhiSignal, int& BsJpsiKKSignal, int& BdJpsiKstar, int& BsOther, int& BdOther,  int& other  );
void fillHistograms(double value, TH1F* hBsJpsiPhiSignal, TH1F* hBsOther, TH1F * hOther, TH1F *hBsJPsiKK = 0);
void setBranchAddresses(TTree *tree);
void bookHistos();
void writeHistos(TString outputfilename);
void identifyEventType(bool & isGenBsJpsiPhiMuMuKKEvent, bool &isGenBsJpsiKKEvent, bool &isGenBdJpsiKstarMuMuKpiEvent, bool &isGenBsEvent_, bool &isGenBdEvent_);
void printDecayTree();

// define cuts:
const double minProbability = 0.02;
const double KaonPtCut = 0.9;
const double nominalPhiMass = 1.0195 ;
const double phiMassWindow = 0.007;
const double decayLengthCut = 3;
const double pointingCut = 0.8;

// maximum number of event to loop over. -1 means all events.
const int maxNumEvents = -1;

// ***************main function*****************************
void BsAnalysisNew() {
  
  TFile *f = new TFile("BtoJPsiMuMu_Nov30.root");
  //    TFile *f = new TFile("ppmumuX.root");
  //  TFile *f = new TFile("ppmuX.root");
  //   TFile *f = new TFile("JPsiMuMu_14Oct09.root");
  
  TString outputhistofilename = "BtoJPsiMuMuHistoFile.root";
  //    TString outputhistofilename = "ppmumuXHistoFile.root";
  //        TString outputhistofilename = "ppmuXHistoFile.root";
  /// TString outputhistofilename = "JPsiMuMuHistoFile.root";
  

  // booking histograms
  bookHistos();  

 
  // get tree from file
  TTree *tree = (TTree*)f->Get("BsTree");
  
  // call function to set branch addresses
  setBranchAddresses(tree);

  int nentries = tree->GetEntries();
  std::cout<<"Number of entries in tree = " << nentries << std::endl;

  int maxEvents;
  if(maxNumEvents < 0) maxEvents = nentries;
  else maxEvents = maxNumEvents;
  //***********************************
  // STARTING event loop here *********
  //***********************************
  for(int i=0; i < maxEvents; i++){
    entrycounter++;
    if( (entrycounter % 10000) == 0) std::cout<<"processing entry number " << entrycounter << std::endl;
    tree->GetEntry(i);
    
    // identify event type
    identifyEventType(isGenBsJpsiPhiMuMuKKEvent_, isGenBsJpsiKKEvent_, isGenBdJpsiKstarMuMuKpiEvent_, isGenBsEvent_, isGenBdEvent_);

    // print decay tree
    //    printDecayTree();

   
    // count number of events of each type
    increaseCountersByOne(iBsJPsiPhiSignalEvents, iBsJPsiKKSignalEvents, iBdJPsiKstarSignalEvents, iBsOtherEvents, iBdOtherEvents, iOtherEvents);
    
    if(triggerbit_HLTdoubleMu3_ != 1) continue;      // IS TRIGGERED_HLTdoubleMu3**************
    HLT++;
    
    // count number of triggered events of each type
    increaseCountersByOne(iTriggered_DoubleMu3_BsJPsiPhiSignalEvents,  iTriggered_DoubleMu3_BsJpsiKK, 
			  iTriggered_DoubleMu3_BdJpsiKstar,            iTriggered_DoubleMu3_BsOtherEvents, 
			  iTriggered_DoubleMu3_BdOtherEvents,          iTriggered_DoubleMu3_OtherEvents     );
    
     
    if(iPassedCutIdent_ < 7) continue;   // passed  pre KinFit mass cuts
    increaseCountersByOne(iBsJPsiPhiSignalEventsPreKinFit, iBsJPsiKKSignalEventsPreKinFit, iBdJPsiKstarSignalEventsPreKinFit,
			  iBsOtherEventsPreKinFit, iBdOtherEventsPreKinFit, iOtherEventsPreKinFit);
    
    
    if(iPassedCutIdent_ < 13) continue; // passed kinematic fit  
    increaseCountersByOne(iBsJPsiPhiSignalEventsOfflineSel1, iBsJPsiKKSignalEventsOfflineSel1, iBdJPsiKstarSignalEventsOfflineSel1,
			  iBsOtherEventsOfflineSel1, iBdOtherEventsOfflineSel1, iOtherEventsOfflineSel1 );
    
  
    // successful fit, fill histograms
    fillHistograms(JpsiMass_after_,     hJPsiMass_BsJpsiSignal, hJPsiMass_BsOther, hJPsiMass_Other); 
    fillHistograms(PhiMass_after_,      hPhiMass_BsJpsiSignal, hPhiMass_BsOther, hPhiMass_Other); 
    fillHistograms(Dist3d_,             hDist3D_BsJpsiSignal, hDist3D_BsOther, hDist3D_Other); 
    fillHistograms(chi2_Bs_,            hChi2_BsJpsiSignal, hChi2_BsOther, hChi2_Other);
    fillHistograms(chi2_Bs_ / ndof_Bs_ , hChi2Ndof_BsJpsiSignal, hChi2Ndof_BsOther, hChi2Ndof_Other);
    fillHistograms(BsVtxProb_,          hBsVtxProb_BsJpsiSignal,  hBsVtxProb_BsOther,  hBsVtxProb_Other );
    fillHistograms(K1Pt_after_,         hK1Pt_BsJpsiSignal, hK1Pt_BsOther, hK1Pt_Other);
    fillHistograms(K2Pt_after_,         hK2Pt_BsJpsiSignal, hK2Pt_BsOther, hK2Pt_Other);
    fillHistograms(Dist3d_ / dDist3d_ , hDistSign3D_BsJpsiSignal, hDistSign3D_BsOther, hDistSign3D_Other);
    fillHistograms(Dist_ / dDist_    ,  hDistSign1D_BsJpsiSignal, hDistSign1D_BsOther, hDistSign1D_Other);
    fillHistograms(AngleBsDecayLength_, hPointingAngle_BsJpsiSignal, hPointingAngle_BsOther, hPointingAngle_Other );


    if(BsVtxProb_ < minProbability) continue; // passed minimum vertex probability cut
    increaseCountersByOne(iBsJPsiPhiSignalEventsProbVertex, iBsJPsiKKSignalEventsProbVertex, iBdJPsiKstarSignalEventsProbVertex,
			  iBsOtherEventsProbVertex, iBdOtherEventsProbVertex, iOtherEventsProbVertex);
    
    if(K1Pt_after_ < KaonPtCut  || K2Pt_after_ < KaonPtCut) continue;   // passed Kaon pt cut
    increaseCountersByOne(iBsJPsiPhiSignalEventsKaonPtCut, iBsJPsiKKSignalEventsKaonPtCut, iBdJPsiKstarSignalEventsKaonPtCut, 
			  iBsOtherEventsKaonPtCut, iBdOtherEventsKaonPtCut, iOtherEventsKaonPtCut);
    
    if( fabs(PhiMass_after_ - nominalPhiMass) > phiMassWindow) continue;   // passed tight phi mass cut
    increaseCountersByOne(iBsJPsiPhiSignalEventsPhiMassCut, iBsJPsiKKSignalEventsPhiMassCut, iBdJPsiKstarSignalEventsPhiMassCut, 
			  iBsOtherEventsPhiMassCut, iBdOtherEventsPhiMassCut, iOtherEventsPhiMassCut);
    
    if(Dist_/ dDist_ < decayLengthCut) continue;  // passed transverse decay length cut
    increaseCountersByOne(iBsJPsiPhiSignalEventsDecayLengthCut, iBsJPsiKKSignalEventsDecayLengthCut, iBdJPsiKstarSignalEventsDecayLengthCut, 
			  iBsOtherEventsDecayLengthCut, iBdOtherEventsDecayLengthCut, iOtherEventsDecayLengthCut);
    
    if(AngleBsDecayLength_ < pointingCut) continue;  // passed pointing constraint
    increaseCountersByOne(iBsJPsiPhiSignalEventsPointingCut, iBsJPsiKKSignalEventsPointingCut, iBdJPsiKstarSignalEventsPointingCut, 
			  iBsOtherEventsPointingCut, iBdOtherEventsPointingCut, iOtherEventsPointingCut);
    
    // make plots after all cuts
    fillHistograms(PhiMass_after_,    hPhiMassFinal_BsJpsiSignal, hPhiMassFinal_BsOther, hPhiMassFinal_Other);
    fillHistograms(  BsMass_after_ ,  hBsMassFinal_BsJpsiSignal, hBsMassFinal_BsOther, hBsMassFinal_Other);
    fillHistograms( BfitM_KK_,        hBsMassFinalAfterFit_BsJpsiSignal, hBsMassFinalAfterFit_BsOther, 
		                      hBsMassFinalAfterFit_Other, hBsMassFinalAfterFit_BsJpsiKK);
  }
  //****************************
  // finished event loop *******
  //****************************
  
  
  // print out results
  std::cout<< "There are " << entrycounter << " events in the tree." << std::endl;
 std::cout<< "Events with only HLT: " << HLT << std::endl << std::endl;
 
 std::cout<< "                          | Bs->JpsiPhi | Bd->Jpsi K* | Bs->Jpsi KK | Other Bs | Other Bd | Other "<<  std::endl;
 std::cout<< "events num                |    " << iBsJPsiPhiSignalEvents << "     |    " << iBdJPsiKstarSignalEvents  
	  << "    |    " << iBsJPsiKKSignalEvents <<"     |   " << iBsOtherEvents << " |   " << iBdOtherEvents << " |   " << iOtherEvents<< std::endl;
 
 std::cout<< "Triggered  events         |     " << iTriggered_DoubleMu3_BsJPsiPhiSignalEvents << "     |     " 
	  << iTriggered_DoubleMu3_BdJpsiKstar  << "    |     " << iTriggered_DoubleMu3_BsJpsiKK << "     |     " 
	  << iTriggered_DoubleMu3_BsOtherEvents << " |    " << iTriggered_DoubleMu3_BdOtherEvents <<  " |    " << iTriggered_DoubleMu3_OtherEvents << std::endl;
 
 std::cout<< "Pre   Kinem  fit*         |      " << iBsJPsiPhiSignalEventsPreKinFit << "     |      " 
	  << iBdJPsiKstarSignalEventsPreKinFit << "    |     " << iBsJPsiKKSignalEventsPreKinFit << "     |     " 
	  << iBsOtherEventsPreKinFit  << " |    " << iBdOtherEventsPreKinFit <<  " |    " << iOtherEventsPreKinFit << std::endl;
 
 std::cout<< "After Kinem. fit.         |      " << iBsJPsiPhiSignalEventsOfflineSel1 << "     |      " << iBdJPsiKstarSignalEventsOfflineSel1 
	  << "    |      " << iBsJPsiKKSignalEventsOfflineSel1 << "     |      " << iBsOtherEventsOfflineSel1 << " |     " 
	  << iBdOtherEventsOfflineSel1 <<  " |     "   << iOtherEventsOfflineSel1 << std::endl; 

 std::cout<< "Prob vertex(KK)>2 percent |      " << iBsJPsiPhiSignalEventsProbVertex <<  "     |       " <<  iBdJPsiKstarSignalEventsProbVertex 
	  << "    |       " << iBsJPsiKKSignalEventsProbVertex << "     |      " << iBsOtherEventsProbVertex << " |      " 
	  << iBdOtherEventsProbVertex     << " |      "   << iOtherEventsProbVertex     <<std::endl;

  std::cout<< "Kaon pt > 0.9            |      " << iBsJPsiPhiSignalEventsKaonPtCut <<  "     |       " <<  iBdJPsiKstarSignalEventsKaonPtCut 
	  << "    |       " << iBsJPsiKKSignalEventsKaonPtCut << "     |      " << iBsOtherEventsKaonPtCut << " |      " 
	  << iBdOtherEventsKaonPtCut     << " |      "   << iOtherEventsKaonPtCut     <<std::endl;

 std::cout<< "deltaPhiMass< 7MeV        |      " << iBsJPsiPhiSignalEventsPhiMassCut <<  "     |       " <<  iBdJPsiKstarSignalEventsPhiMassCut 
	  << "    |       " << iBsJPsiKKSignalEventsPhiMassCut << "     |      " << iBsOtherEventsPhiMassCut << " |      " 
	  << iBdOtherEventsPhiMassCut     << " |      "   << iOtherEventsPhiMassCut     <<std::endl;

 std::cout<< "decayLengthSign > 3       |      " << iBsJPsiPhiSignalEventsDecayLengthCut <<  "     |       " <<  iBdJPsiKstarSignalEventsDecayLengthCut 
	  << "    |       " << iBsJPsiKKSignalEventsDecayLengthCut << "     |      " << iBsOtherEventsDecayLengthCut << " |      " 
	  << iBdOtherEventsDecayLengthCut     << " |      "   << iOtherEventsDecayLengthCut     <<std::endl;

 std::cout<< "Pointing Cut             |      " << iBsJPsiPhiSignalEventsPointingCut <<  "     |       " <<  iBdJPsiKstarSignalEventsPointingCut 
	  << "    |       " << iBsJPsiKKSignalEventsPointingCut << "     |      " << iBsOtherEventsPointingCut << " |      " 
	  << iBdOtherEventsPointingCut     << " |      "   << iOtherEventsPointingCut     <<std::endl;


 std::cout<< "  "<<  std::endl;
 std::cout<< "  "<<  std::endl;

 std::cout<<"alternative (simple) event category definition:"<< std::endl;
   std::cout<< "                          | Bs->JpsiPhi |  Other Bs  | Other "<<  std::endl;
 std::cout<< "events num                |    " << iBsJPsiPhiSignalEvents << "     |    " << iBsJPsiKKSignalEvents + iBsOtherEvents  
	  << " |   " << iBdOtherEvents + iOtherEvents + iBdJPsiKstarSignalEvents << std::endl;
 
 std::cout<< "Triggered  events         |     " << iTriggered_DoubleMu3_BsJPsiPhiSignalEvents << "     |     " 
	  << iTriggered_DoubleMu3_BsOtherEvents + iTriggered_DoubleMu3_BsJpsiKK<< " |    " 
	  << iTriggered_DoubleMu3_BdOtherEvents+  iTriggered_DoubleMu3_BdJpsiKstar + iTriggered_DoubleMu3_OtherEvents << std::endl;
 
 std::cout<< "Pre   Kinem  fit*         |      " << iBsJPsiPhiSignalEventsPreKinFit << "     |      " 
	  << iBsJPsiKKSignalEventsPreKinFit + iBsOtherEventsPreKinFit  << " |    " 
	  << iBdOtherEventsPreKinFit + iBdJPsiKstarSignalEventsPreKinFit+ iOtherEventsPreKinFit << std::endl;
 
 std::cout<< "After Kinem. fit.         |      " << iBsJPsiPhiSignalEventsOfflineSel1 << "     |      " 
	  << iBsJPsiKKSignalEventsOfflineSel1 + iBsOtherEventsOfflineSel1 << " |     " 
	  << iBdOtherEventsOfflineSel1 + iOtherEventsOfflineSel1 + iBdJPsiKstarSignalEventsOfflineSel1<< std::endl; 

 std::cout<< "Prob vertex(KK)>2 percent |      " << iBsJPsiPhiSignalEventsProbVertex <<  "     |       "  
	  << iBsJPsiKKSignalEventsProbVertex + iBsOtherEventsProbVertex << " |      " 
	  << iBdOtherEventsProbVertex   + iOtherEventsProbVertex +iBdJPsiKstarSignalEventsProbVertex    <<std::endl;

 std::cout<< "Kaon pt > 0.9            |      " << iBsJPsiPhiSignalEventsKaonPtCut <<  "     |       "  
	  << iBsJPsiKKSignalEventsKaonPtCut + iBsOtherEventsKaonPtCut << " |      " 
	  << iBdOtherEventsKaonPtCut +iOtherEventsKaonPtCut + iBdJPsiKstarSignalEventsKaonPtCut    <<std::endl;

 std::cout<< "deltaPhiMass< 7MeV        |      " << iBsJPsiPhiSignalEventsPhiMassCut <<  "     |       "  
	  << iBsJPsiKKSignalEventsPhiMassCut + iBsOtherEventsPhiMassCut << " |      " 
	  << iBdOtherEventsPhiMassCut + iOtherEventsPhiMassCut +iBdJPsiKstarSignalEventsPhiMassCut    <<std::endl;

 std::cout<< "decayLengthSign > 3       |      " << iBsJPsiPhiSignalEventsDecayLengthCut <<  "     |       "   
	  << iBsJPsiKKSignalEventsDecayLengthCut + iBsOtherEventsDecayLengthCut << " |      " 
	  << iBdOtherEventsDecayLengthCut  + iOtherEventsDecayLengthCut  +iBdJPsiKstarSignalEventsDecayLengthCut   <<std::endl;

 std::cout<< "Pointing Cut             |      " << iBsJPsiPhiSignalEventsPointingCut <<  "     |       " 
	  << iBsJPsiKKSignalEventsPointingCut + iBsOtherEventsPointingCut << " |      " 
	  << iBdOtherEventsPointingCut  + iOtherEventsPointingCut  + iBdJPsiKstarSignalEventsPointingCut    <<std::endl;



 // write histograms to output file
 writeHistos(outputhistofilename);

}   // finished main function



// define counting function
void increaseCountersByOne(int &BsJpsiPhiSignal, int& BsJpsiKKSignal, int& BdJpsiKstar, int& BsOther, int& BdOther,  int& other  ){

  if(isGenBsJpsiPhiMuMuKKEvent_ == 1)               BsJpsiPhiSignal++;
  else if (isGenBsJpsiKKEvent_ == 1)                BsJpsiKKSignal++;
       else if (isGenBdJpsiKstarMuMuKpiEvent_ == 1) BdJpsiKstar++;
            else if (isGenBsEvent_ == 1)            BsOther++;
	         else if (isGenBdEvent_ == 1)       BdOther++;		 
		      else                          other++; 

}  // end definition of counting function



// histogram filling function
void fillHistograms(double value, TH1F* hBsJpsiPhiSignal, TH1F* hBsOther, TH1F * hOther, TH1F *hBsJPsiKK ){

  if(isGenBsJpsiPhiMuMuKKEvent_ == 1) hBsJpsiPhiSignal->Fill(value);
  else if(isGenBsEvent_==1) hBsOther->Fill(value);
  else hOther->Fill(value);
  
  // special case to check distributions of BsJPsiKK events
  if(isGenBsJpsiKKEvent_ == 1)  if(hBsJPsiKK) hBsJPsiKK->Fill(value);
} // end histogram filling function



// find type of events based on MC info arrays
void identifyEventType(bool & isGenBsJpsiPhiMuMuKKEvent, bool &isGenBsJpsiKKEvent, bool &isGenBdJpsiKstarMuMuKpiEvent, bool &isGenBsEvent, bool &isGenBdEvent){

  //set to default values
  isGenBsJpsiPhiMuMuKKEvent = isGenBsJpsiKKEvent = isGenBdJpsiKstarMuMuKpiEvent = isGenBsEvent = isGenBdEvent = 0;

 

  //loop over B mesons
  for(int i=0; i < GenNumberOfBdecays_; i++){
    int B_id = abs(BmesonsId_[i]);
    if(B_id== 531) isGenBsEvent = 1;
    if(B_id == 511) isGenBdEvent = 1;
    
    //loop over daughters
    // helper flags
    bool BsJpsi = 0;
    bool BsPhi = 0;
    bool BsJpsiMuPlus = 0;
    bool BsJpsiMuMinus = 0;
    bool BsPhiKplus = 0;
    bool BsPhiKminus = 0;
    bool BsKplus = 0;
    bool BsKminus = 0;

    bool BdJpsi = 0;
    bool BdJpsiMuPlus = 0;
    bool BdJpsiMuMinus = 0;
    bool BdKstar = 0;
    bool BdKstarK = 0;
    bool BdKstarPi = 0;

   

    for(int j=0; j < GenNumberOfDaughters_[i]; j++){
      int Dau_id = BDauIdMC_[i][j];
      int absDau_id = abs(Dau_id);
      if(B_id == 531 && absDau_id == 443) BsJpsi =1;
      if(B_id == 531 && absDau_id == 333) BsPhi = 1;
      if(B_id == 531 && Dau_id == 321)    BsKplus = 1;
      if(B_id == 531 && Dau_id == -321)   BsKminus = 1;
      
      if(B_id == 511 && absDau_id == 443) BdJpsi = 1;
      if(B_id == 511 && absDau_id == 313) BdKstar = 1;   //// ATTENTION: add other excited Kaons??

      // loop over daughters daughters
      for(int k=0; k < GenNumberOfDaughtersDaughters_[i][j]; k++){
	int DauDau_id = BDauDauIdMC_[i][j][k];
	int absDauDau_id= abs(DauDau_id);
	if(B_id == 531 && absDau_id == 443 && DauDau_id == 13) BsJpsiMuMinus = 1;
	if(B_id == 531 && absDau_id == 443 && DauDau_id == -13) BsJpsiMuPlus = 1;
	if(B_id == 531 && absDau_id == 333 && DauDau_id == 321) BsPhiKplus = 1;
	if(B_id == 531 && absDau_id == 333 && DauDau_id == -321) BsPhiKminus = 1;

	if(B_id == 511 && absDau_id == 443 && DauDau_id == 13) BdJpsiMuMinus = 1;
	if(B_id == 511 && absDau_id == 443 && DauDau_id == -13) BdJpsiMuPlus = 1;
	if(B_id == 511 && absDau_id == 313 && absDauDau_id == 321) BdKstarK = 1;      ////add other excited Kaons???
	if(B_id == 511 && absDau_id == 313 && absDauDau_id == 211) BdKstarPi = 1;
      }
    }
 
    if(BsJpsiMuMinus && BsJpsiMuPlus && BsPhiKplus && BsPhiKminus) isGenBsJpsiPhiMuMuKKEvent = 1;
    if(BsJpsiMuMinus && BsJpsiMuPlus && BsKplus && BsKminus)       isGenBsJpsiKKEvent = 1;  
    if(BdJpsiMuMinus && BdJpsiMuPlus && BdKstarPi && BdKstarK )    isGenBdJpsiKstarMuMuKpiEvent = 1;
  }
   
}

//----------------
void printDecayTree(){
  for(int i=0; i<GenNumberOfBdecays_; i++){
    std::cout<<"B pdg id = " << BmesonsId_[i] << std::endl;
    for(int j=0; j<GenNumberOfDaughters_[i]; j++){
      std::cout<<"  Dau pdg id = " << BDauIdMC_[i][j] << std::endl;
      for(int k=0; k<GenNumberOfDaughtersDaughters_[i][j]; k++){
	std::cout<<"      DauDau pdg id = " << BDauDauIdMC_[i][j][k] << std::endl;
      }
    }
  }
}
//-------------------

// function to set branch addresses
void setBranchAddresses(TTree *tree){
  
   tree->SetBranchAddress("triggerbit_HLTmu3",                                 &triggerbit_HLTmu3_);    
   tree->SetBranchAddress("triggerbit_HLTmu5",				       &triggerbit_HLTmu5_);                                                  
   tree->SetBranchAddress("triggerbit_HLTmu9",				       &triggerbit_HLTmu9_);                                                  
   tree->SetBranchAddress("triggerbit_HLTdoubleIsoMu3",			       &triggerbit_HLTdoubleIsoMu3_);                                         
   tree->SetBranchAddress("triggerbit_HLTdoubleMu3",			       &triggerbit_HLTdoubleMu3_);                                            
   tree->SetBranchAddress("triggerbit_HLTdoubleMu3_JPsi",		       &triggerbit_HLTdoubleMu3_JPsi_);                                       
   tree->SetBranchAddress("BeamspotX" ,					       &BSx_) ;                                                               
   tree->SetBranchAddress("BeamspotY" ,					       &BSy_ );                                                               
   tree->SetBranchAddress("BeamspotZ" ,					       &BSz_ );                                                               
   tree->SetBranchAddress("PVx" ,					       &PVx_ );                                                               
   tree->SetBranchAddress("PVy" ,					       &PVy_ );                                                               
   tree->SetBranchAddress("PVz" ,					       &PVz_ );                                                               
   tree->SetBranchAddress("PVerrx" ,					       &PVerrx_ );                                                            
   tree->SetBranchAddress("PVerry" ,					       &PVerry_ );                                                            
   tree->SetBranchAddress("PVerrz" ,					       &PVerrz_ );                                                            
   tree->SetBranchAddress("JpsiVtxProb",				       &JpsiVtxProb_);                                                        
   tree->SetBranchAddress("JpsiM_alone",				       &JpsiM_alone_);                                                        
   tree->SetBranchAddress("JpsiPhi_alone",				       &JpsiPhi_alone_);                                                      
   tree->SetBranchAddress("JpsiEta_alone",				       &JpsiEta_alone_);                                                      
   tree->SetBranchAddress("JpsiPt_alone",				       &JpsiPt_alone_);                                                       
   tree->SetBranchAddress("JpsiMu1Pt",					       &JpsiMu1Pt_);                                                          
   tree->SetBranchAddress("JpsiMu2Pt",					       &JpsiMu2Pt_);                                                          
   tree->SetBranchAddress("JpsiMu1Phi",					       &JpsiMu1Phi_);                                                         
   tree->SetBranchAddress("JpsiMu2Phi",					       &JpsiMu2Phi_);                                                         
   tree->SetBranchAddress("JpsiMu1Eta",					       &JpsiMu1Eta_);                                                         
   tree->SetBranchAddress("JpsiMu2Eta",					       &JpsiMu2Eta_);                                                         
   tree->SetBranchAddress("JpsiMuon1Cat",				       &JpsiMuon1Cat_);                                                       
   tree->SetBranchAddress("JpsiMuon2Cat",				       &JpsiMuon2Cat_);                                                       
   tree->SetBranchAddress("BsMass_before",				       &BsMass_before_);                                                      
   tree->SetBranchAddress("BsPhi_before",				       &BsPhi_before_);                                                       
   tree->SetBranchAddress("BsEta_before",				       &BsEta_before_);                                                       
   tree->SetBranchAddress("BsPt_before",				       &BsPt_before_);                                                        
   tree->SetBranchAddress("BsPz_before",				       &BsPz_before_);                                                        
   tree->SetBranchAddress("JpsiMass_before",				       &JpsiMass_before_);                                                    
   tree->SetBranchAddress("JpsiPhi_before",				       &JpsiPhi_before_);                                                     
   tree->SetBranchAddress("JpsiEta_before",				       &JpsiEta_before_);                                                     
   tree->SetBranchAddress("JpsiPt_before",				       &JpsiPt_before_);                                                      
   tree->SetBranchAddress("JpsiPz_before",				       &JpsiPz_before_);                                                      
   tree->SetBranchAddress("PhiMass_before",				       &PhiMass_before_);                                                     
   tree->SetBranchAddress("PhiPhi_before",				       &PhiPhi_before_);                                                      
   tree->SetBranchAddress("PhiEta_before",				       &PhiEta_before_);                                                      
   tree->SetBranchAddress("PhiPt_before",				       &PhiPt_before_);                                                       
   tree->SetBranchAddress("PhiPz_before",				       &PhiPz_before_);                                                       
   tree->SetBranchAddress("K1Pt_before",				       &K1Pt_before_);                                                        
   tree->SetBranchAddress("K1Pz_before",				       &K1Pz_before_);                                                        
   tree->SetBranchAddress("K1Eta_before",				       &K1Eta_before_);                                                       
   tree->SetBranchAddress("K1Phi_before",				       &K1Phi_before_);                                                       
   tree->SetBranchAddress("K2Eta_before",				       &K2Eta_before_);                                                       
   tree->SetBranchAddress("K2Pt_before",				       &K2Pt_before_);                                                        
   tree->SetBranchAddress("K2Pz_before",				       &K2Pz_before_);                                                        
   tree->SetBranchAddress("K2Phi_before",				       &K2Phi_before_);                                                       
   tree->SetBranchAddress("chi2_Bs",					       &chi2_Bs_);                                                            
   tree->SetBranchAddress("ndof_Bs",					       &ndof_Bs_);                                                            
   tree->SetBranchAddress("BsVtxProb",					       &BsVtxProb_);                                                          
   tree->SetBranchAddress("BsVtxProbKpi",				       &BsVtxProbKpi_);                                                       
   tree->SetBranchAddress("BsVtxProbpipi",				       &BsVtxProbpipi_);                                                      
   tree->SetBranchAddress("BFitM_KK",					       &BfitM_KK_);                                                           
   tree->SetBranchAddress("BFitM_Kpi",					       &BfitM_Kpi_);                                                          
   tree->SetBranchAddress("BFitM_pipi",					       &BfitM_pipi_);                                                         
   tree->SetBranchAddress("BsVtx_x",					       &BsVtx_x_);                                                            
   tree->SetBranchAddress("BsVtx_y",					       &BsVtx_y_);                                                            
   tree->SetBranchAddress("BsVtx_z",					       &BsVtx_z_);                                                            
   tree->SetBranchAddress("BsMass_after",				       &BsMass_after_);                                                       
   tree->SetBranchAddress("BsPhi_after",				       &BsPhi_after_);                                                        
   tree->SetBranchAddress("BsEta_after",				       &BsEta_after_);                                                        
   tree->SetBranchAddress("BsPt_after",					       &BsPt_after_);                                                         
   tree->SetBranchAddress("BsPz_after",					       &BsPz_after_);                                                         
   tree->SetBranchAddress("JpsiMass_after",				       &JpsiMass_after_);                                                     
   tree->SetBranchAddress("JpsiPhi_after",				       &JpsiPhi_after_);                                                      
   tree->SetBranchAddress("JpsiEta_after",				       &JpsiEta_after_);                                                      
   tree->SetBranchAddress("JpsiPt_after",				       &JpsiPt_after_);                                                       
   tree->SetBranchAddress("JpsiPz_after",				       &JpsiPz_after_);                                                       
   tree->SetBranchAddress("PhiMass_after",				       &PhiMass_after_);                                                      
   tree->SetBranchAddress("PhiPhi_after",				       &PhiPhi_after_);                                                       
   tree->SetBranchAddress("PhiEta_after",				       &PhiEta_after_);                                                       
   tree->SetBranchAddress("PhiPt_after",				       &PhiPt_after_);                                                        
   tree->SetBranchAddress("PhiPz_after",				       &PhiPz_after_);                                                        
   tree->SetBranchAddress("K1Pt_after",					       &K1Pt_after_);                                                         
   tree->SetBranchAddress("K1Pz_after",					       &K1Pz_after_);                                                         
   tree->SetBranchAddress("K1Eta_after",				       &K1Eta_after_);                                                        
   tree->SetBranchAddress("K1Phi_after",				       &K1Phi_after_);                                                        
   tree->SetBranchAddress("K2Eta_after",				       &K2Eta_after_);                                                        
   tree->SetBranchAddress("K2Pt_after",					       &K2Pt_after_);                                                         
   tree->SetBranchAddress("K2Pz_after",					       &K2Pz_after_);                                                         
   tree->SetBranchAddress("K2Phi_after",				       &K2Phi_after_);                                                        
   tree->SetBranchAddress("K1Chi2",					       &K1Chi2_);                                                             
//    tree->SetBranchAddress("K1nHits",					       &K1nHits_);                                                            
//    tree->SetBranchAddress("K2Chi2",					       &K2Chi2_);                                                             
//    tree->SetBranchAddress("K2nHits",					       &K2nHits_);                                                            
//    tree->SetBranchAddress("K1pixH",					       &K1pixH_);                                                             
//    tree->SetBranchAddress("K1trkH",					       &K1trkH_);                                                             
//    tree->SetBranchAddress("K2pixH",					       &K2pixH_);                                                             
//    tree->SetBranchAddress("K2trkH",					       &K2trkH_);                                                             
//    tree->SetBranchAddress("Mu1Chi2",					       &Mu1Chi2_);                                                            
//    tree->SetBranchAddress("Mu1nHits",					       &Mu1nHits_);                                                           
//    tree->SetBranchAddress("Mu2Chi2",					       &Mu2Chi2_);                                                            
//    tree->SetBranchAddress("Mu2nHits",					       &Mu2nHits_);                                                           
//    tree->SetBranchAddress("Mu1pixH",					       &Mu1pixH_);                                                            
//    tree->SetBranchAddress("Mu1trkH",					       &Mu1trkH_);                                                            
//    tree->SetBranchAddress("Mu2pixH",					       &Mu2pixH_);                                                            
//    tree->SetBranchAddress("Mu2trkH",					       &Mu2trkH_);                                                            
   tree->SetBranchAddress("costheta",					       &costheta_);                                                           
   tree->SetBranchAddress("phi",					       &phi_);                                                                
   tree->SetBranchAddress("cospsi",					       &cospsi_);                                                             
   tree->SetBranchAddress("AngleBsDecayLength",				       &AngleBsDecayLength_);                                                 
   tree->SetBranchAddress("isMatched",					       &isMatched_);                                                          
   tree->SetBranchAddress("isMatchedBd",				       &isMatchedBd_);                                                        
   tree->SetBranchAddress("BLxy",					       &BLxy_);                                                               
   tree->SetBranchAddress("BLxy2",					       &BLxy2_);                                                              
   tree->SetBranchAddress("BerrX",					       &BerrX_);                                                              
   tree->SetBranchAddress("BerrY",					       &BerrY_);                                                              
   tree->SetBranchAddress("BerrXY",					       &BerrXY_);                                                             
   tree->SetBranchAddress("Bsct1",					       &Bsct1_);                                                              
   tree->SetBranchAddress("Bsct2",					       &Bsct2_);                                                              
   tree->SetBranchAddress("K1trkLay",					       &K1trkLay_);                                                           
   tree->SetBranchAddress("K1muDTh",					       &K1muDTh_);                                                            
   tree->SetBranchAddress("K1muCSCh",					       &K1muCSCh_);                                                           
   tree->SetBranchAddress("K1muRPCh",					       &K1muRPCh_);                                                           
   tree->SetBranchAddress("K2trkLay",					       &K2trkLay_);                                                           
   tree->SetBranchAddress("K2muDTh",					       &K2muDTh_);                                                            
   tree->SetBranchAddress("K2muCSCh",					       &K2muCSCh_);                                                           
   tree->SetBranchAddress("K2muRPCh",					       &K2muRPCh_);                                                           
   tree->SetBranchAddress("Mu1trkLay",					       &Mu1trkLay_);                                                          
   tree->SetBranchAddress("Mu1muDTh",					       &Mu1muDTh_);                                                           
   tree->SetBranchAddress("Mu1muCSCh",					       &Mu1muCSCh_);                                                          
   tree->SetBranchAddress("Mu1muRPCh",					       &Mu1muRPCh_);                                                          
   tree->SetBranchAddress("Mu2trkLay",					       &Mu2trkLay_);                                                          
   tree->SetBranchAddress("Mu2muDTh",					       &Mu2muDTh_);                                                           
   tree->SetBranchAddress("Mu2muCSCh",					       &Mu2muCSCh_);                                                          
   tree->SetBranchAddress("Mu2muRPCh",					       &Mu2muRPCh_);                                                          
   tree->SetBranchAddress("K1mcId",					       &K1mcId_);                                                             
   tree->SetBranchAddress("K1momId",					       &K1momId_);                                                            
   tree->SetBranchAddress("K1gmomId",					       &K1gmomId_);                                                           
   tree->SetBranchAddress("K2mcId",					       &K2mcId_);                                                             
   tree->SetBranchAddress("K2momId",					       &K2momId_);                                                            
   tree->SetBranchAddress("K2gmomId",					       &K2gmomId_);                                                           
   tree->SetBranchAddress("Mu1mcId",					       &Mu1mcId_);                                                            
   tree->SetBranchAddress("Mu1momId",					       &Mu1momId_);                                                           
   tree->SetBranchAddress("Mu1gmomId",					       &Mu1gmomId_);                                                          
   tree->SetBranchAddress("Mu2mcId",					       &Mu2mcId_);                                                            
   tree->SetBranchAddress("Mu2momId",					       &Mu2momId_);                                                           
   tree->SetBranchAddress("Mu2gmomId",					       &Mu2gmomId_);                                                          
   tree->SetBranchAddress("K1Truth",					       &K1Truth_);                                                            
   tree->SetBranchAddress("K2Truth",					       &K2Truth_);                                                            
   tree->SetBranchAddress("Mu1Truth",					       &Mu1Truth_);                                                           
   tree->SetBranchAddress("Mu2Truth",					       &Mu2Truth_);                                                           
   tree->SetBranchAddress("Dist3d",					       &Dist3d_);                                                             
   tree->SetBranchAddress("dDist3d",					       &dDist3d_);                                                            
   tree->SetBranchAddress("Time3d",					       &Time3d_);                                                             
   tree->SetBranchAddress("dTime3d",					       &dTime3d_);                                                            
   tree->SetBranchAddress("Dist",					       &Dist_);                                                               
   tree->SetBranchAddress("dDist",					       &dDist_);                                                              
   tree->SetBranchAddress("Time",					       &Time_);                                                               
   tree->SetBranchAddress("dTime",					       &dTime_);                                                              
   tree->SetBranchAddress("dedxTrk",					       &dedxTrk_);                                                            
   tree->SetBranchAddress("errdedxTrk",					       &errdedxTrk_);                                                         
   tree->SetBranchAddress("numdedxTrk",					       &numdedxTrk_);                                                         
   tree->SetBranchAddress("iPassedCutIdent",				       &iPassedCutIdent_);                                                    
   tree->SetBranchAddress("iPassedCutIdentBd",				       &iPassedCutIdentBd_);                                                  
 


 tree->SetBranchAddress("K1_par0KK", &K1_kk_par0_ );
  tree->SetBranchAddress("K1_par1KK", &K1_kk_par1_);
  tree->SetBranchAddress("K1_par2KK", &K1_kk_par2_);
  tree->SetBranchAddress("K1_par3KK", &K1_kk_par3_);
  tree->SetBranchAddress("K1_par4KK", &K1_kk_par4_);
  tree->SetBranchAddress("K1_par5KK", &K1_kk_par5_);
  tree->SetBranchAddress("K1_par6KK", &K1_kk_par6_);

  tree->SetBranchAddress("K2_par0KK", &K2_kk_par0_);
  tree->SetBranchAddress("K2_par1KK", &K2_kk_par1_);
  tree->SetBranchAddress("K2_par2KK", &K2_kk_par2_);
  tree->SetBranchAddress("K2_par3KK", &K2_kk_par3_);
  tree->SetBranchAddress("K2_par4KK", &K2_kk_par4_);
  tree->SetBranchAddress("K2_par5KK", &K2_kk_par5_);
  tree->SetBranchAddress("K2_par6KK", &K2_kk_par6_);

  tree->SetBranchAddress("K1_par0Kpi", &K1_kpi_par0_);
  tree->SetBranchAddress("K1_par1Kpi", &K1_kpi_par1_);
  tree->SetBranchAddress("K1_par2Kpi", &K1_kpi_par2_);
  tree->SetBranchAddress("K1_par3Kpi", &K1_kpi_par3_);
  tree->SetBranchAddress("K1_par4Kpi", &K1_kpi_par4_);
  tree->SetBranchAddress("K1_par5Kpi", &K1_kpi_par5_);
  tree->SetBranchAddress("K1_par6Kpi", &K1_kpi_par6_);

  tree->SetBranchAddress("K2_par0Kpi", &K2_kpi_par0_);
  tree->SetBranchAddress("K2_par1Kpi", &K2_kpi_par1_);
  tree->SetBranchAddress("K2_par2Kpi", &K2_kpi_par2_);
  tree->SetBranchAddress("K2_par3Kpi", &K2_kpi_par3_);
  tree->SetBranchAddress("K2_par4Kpi", &K2_kpi_par4_);
  tree->SetBranchAddress("K2_par5Kpi", &K2_kpi_par5_);
  tree->SetBranchAddress("K2_par6Kpi", &K2_kpi_par6_);

  tree->SetBranchAddress("K1_par0pipi", &K1_pipi_par0_);
  tree->SetBranchAddress("K1_par1pipi", &K1_pipi_par1_);
  tree->SetBranchAddress("K1_par2pipi", &K1_pipi_par2_);
  tree->SetBranchAddress("K1_par3pipi", &K1_pipi_par3_);
  tree->SetBranchAddress("K1_par4pipi", &K1_pipi_par4_);
  tree->SetBranchAddress("K1_par5pipi", &K1_pipi_par5_);
  tree->SetBranchAddress("K1_par6pipi", &K1_pipi_par6_);

  tree->SetBranchAddress("K2_par0pipi", &K2_pipi_par0_);
  tree->SetBranchAddress("K2_par1pipi", &K2_pipi_par1_);
  tree->SetBranchAddress("K2_par2pipi", &K2_pipi_par2_);
  tree->SetBranchAddress("K2_par3pipi", &K2_pipi_par3_);
  tree->SetBranchAddress("K2_par4pipi", &K2_pipi_par4_);
  tree->SetBranchAddress("K2_par5pipi", &K2_pipi_par5_);
  tree->SetBranchAddress("K2_par6pipi", &K2_pipi_par6_);




                                                  
   tree->SetBranchAddress("K1_kk_sigX_",					       &K1_kk_sigX_);                                                         
   tree->SetBranchAddress("K1_kk_sigY_",					       &K1_kk_sigY_);                                                         
   tree->SetBranchAddress("K1_kk_sigZ_",					       &K1_kk_sigZ_);                                                         
   tree->SetBranchAddress("K1_kpi_sigX_",				       &K1_kpi_sigX_);                                                        
   tree->SetBranchAddress("K1_kpi_sigY_",				       &K1_kpi_sigY_);                                                        
   tree->SetBranchAddress("K1_kpi_sigZ_",				       &K1_kpi_sigZ_);                                                        
   tree->SetBranchAddress("K1_pipi_sigX_",				       &K1_pipi_sigX_);                                                       
   tree->SetBranchAddress("K1_pipi_sigY_",				       &K1_pipi_sigY_);                                                       
   tree->SetBranchAddress("K1_pipi_sigZ_",				       &K1_pipi_sigZ_);                                                       
   tree->SetBranchAddress("K2_kk_sigX_",					       &K2_kk_sigX_);                                                         
   tree->SetBranchAddress("K2_kk_sigY_",					       &K2_kk_sigY_);                                                         
   tree->SetBranchAddress("K2_kk_sigZ_",					       &K2_kk_sigZ_);                                                         
   tree->SetBranchAddress("K2_kpi_sigX_",				       &K2_kpi_sigX_);                                                        
   tree->SetBranchAddress("K2_kpi_sigY_",				       &K2_kpi_sigY_);                                                        
   tree->SetBranchAddress("K2_kpi_sigZ_",				       &K2_kpi_sigZ_);                                                        
   tree->SetBranchAddress("K2_pipi_sigX_",				       &K2_pipi_sigX_);                                                       
   tree->SetBranchAddress("K2_pipi_sigY_",				       &K2_pipi_sigY_);                                                       
   tree->SetBranchAddress("K2_pipi_sigZ_",				       &K2_pipi_sigZ_);                                                       
   tree->SetBranchAddress("K1_kk_sigPX_",				       &K1_kk_sigPX_);                                                        
   tree->SetBranchAddress("K1_kk_sigPY_",				       &K1_kk_sigPY_);                                                        
   tree->SetBranchAddress("K1_kk_sigPZ_",				       &K1_kk_sigPZ_);                                                        
   tree->SetBranchAddress("K1_kpi_sigPX_",				       &K1_kpi_sigPX_);                                                       
   tree->SetBranchAddress("K1_kpi_sigPY_",				       &K1_kpi_sigPY_);                                                       
   tree->SetBranchAddress("K1_kpi_sigPZ_",				       &K1_kpi_sigPZ_);                                                       
   tree->SetBranchAddress("K1_pipi_sigPX_",				       &K1_pipi_sigPX_);                                                      
   tree->SetBranchAddress("K1_pipi_sigPY_",				       &K1_pipi_sigPY_);                                                      
   tree->SetBranchAddress("K1_pipi_sigPZ_",				       &K1_pipi_sigPZ_);                                                      
   tree->SetBranchAddress("K2_kk_sigPX_",				       &K2_kk_sigPX_);                                                        
   tree->SetBranchAddress("K2_kk_sigPY_",				       &K2_kk_sigPY_);                                                        
   tree->SetBranchAddress("K2_kk_sigPZ_",				       &K2_kk_sigPZ_);                                                        
   tree->SetBranchAddress("K2_kpi_sigPX_",				       &K2_kpi_sigPX_);                                                       
   tree->SetBranchAddress("K2_kpi_sigPY_",				       &K2_kpi_sigPY_);                                                       
   tree->SetBranchAddress("K2_kpi_sigPZ_",				       &K2_kpi_sigPZ_);                                                       
   tree->SetBranchAddress("K2_pipi_sigPX_",				       &K2_pipi_sigPX_);                                                      
   tree->SetBranchAddress("K2_pipi_sigPY_",				       &K2_pipi_sigPY_);                                                      
   tree->SetBranchAddress("K2_pipi_sigPZ_",				       &K2_pipi_sigPZ_);                                                      
   tree->SetBranchAddress("sigmaK1Pt",					       &K1Pt_error_);                                                         
   tree->SetBranchAddress("sigmaK2Pt",					       &K2Pt_error_);                                                         
   tree->SetBranchAddress("GenNumberOfBdecays",				       &GenNumberOfBdecays_);                                                 
   tree->SetBranchAddress("BmesonsId",					       BmesonsId_);                                                          
   tree->SetBranchAddress("BDauIdMC",					       BDauIdMC_);                                                           
   tree->SetBranchAddress("BDauDauIdMC",				       BDauDauIdMC_);                                                        
   tree->SetBranchAddress("GenNumberOfDaughters",			       GenNumberOfDaughters_);                                               
   tree->SetBranchAddress("GenNumberOfDaughtersDaughters",		       GenNumberOfDaughtersDaughters_);                                      
   tree->SetBranchAddress("BDauMMC",					       BDauMMC_);                                                            
   tree->SetBranchAddress("BDauPtMC",					       BDauPtMC_);                                                           
   tree->SetBranchAddress("BDauPzMC",					       BDauPzMC_);                                                           
   tree->SetBranchAddress("BDauEtaMC",					       BDauEtaMC_);                                                          
   tree->SetBranchAddress("BDauPhiMC",					       BDauPhiMC_);                                                          
   tree->SetBranchAddress("BDauDauMMC",					       BDauDauMMC_);                                                         
   tree->SetBranchAddress("BDauDauPtMC",				       BDauDauPtMC_);                                                        
   tree->SetBranchAddress("BDauDauPzMC",				       BDauDauPzMC_);                                                        
   tree->SetBranchAddress("BDauDauEtaMC",				       BDauDauEtaMC_);                                                       
   tree->SetBranchAddress("BDauDauPhiMC",				       BDauDauPhiMC_);                                                       
   tree->SetBranchAddress("BMMC",					       BMMC_);                                                               
   tree->SetBranchAddress("BPtMC",					       BPtMC_);                                                              
   tree->SetBranchAddress("BPzMC",					       BPzMC_);                                                              
   tree->SetBranchAddress("BEtaMC",					       BEtaMC_);                                                             
   tree->SetBranchAddress("BPhiMC",					       BPhiMC_);                                                             

   tree->SetBranchAddress("genBsVtx_z",                                       &genBsVtx_z_);
   tree->SetBranchAddress("genBsVtx_y",                                       &genBsVtx_y_);
   tree->SetBranchAddress("genBsVtx_x",                                       &genBsVtx_x_) ;                              
   tree->SetBranchAddress("genBsSVtx_z",  	                               &genBsSVtx_z_) ;                           
   tree->SetBranchAddress("genBsSVtx_y",                                      &genBsSVtx_y_);
   tree->SetBranchAddress("genBsSVtx_x",                                      &genBsSVtx_x_) ; 
   tree->SetBranchAddress("isGenJpsiEvent_",                                    &isGenJpsiEvent_);                                                     
  

}


void writeHistos(TString outputfilename){
  TFile *histofile = new TFile(outputfilename , "RECREATE");
 hJPsiMass_BsJpsiSignal            ->Write();
 hJPsiMass_BsOther                 ->Write();
 hJPsiMass_Other                   ->Write();
                               
 hPhiMass_BsJpsiSignal             ->Write();
 hPhiMass_BsOther                  ->Write();
 hPhiMass_Other                    ->Write();
                          
 hPhiMassFinal_BsJpsiSignal        ->Write();
 hPhiMassFinal_BsOther             ->Write();
 hPhiMassFinal_Other               ->Write();
                                
 hBsMassFinal_BsJpsiSignal         ->Write();
 hBsMassFinal_BsOther              ->Write();
 hBsMassFinal_Other                ->Write();
                                 
 hBsMassFinalAfterFit_BsJpsiSignal ->Write();
 hBsMassFinalAfterFit_BsOther      ->Write();
 hBsMassFinalAfterFit_Other        ->Write();
 hBsMassFinalAfterFit_BsJpsiKK     ->Write();
                            
 hK1Pt_BsJpsiSignal                ->Write();
 hK1Pt_BsOther                     ->Write();
 hK1Pt_Other                       ->Write();
                        
 hK2Pt_BsJpsiSignal                ->Write();
 hK2Pt_BsOther                     ->Write();
 hK2Pt_Other                       ->Write();
                             
 hDist3D_BsJpsiSignal              ->Write();
 hDist3D_BsOther                   ->Write();
 hDist3D_Other                     ->Write();
                       
 hDistSign3D_BsJpsiSignal          ->Write();
 hDistSign3D_BsOther               ->Write();
 hDistSign3D_Other                 ->Write();
                              
 hDistSign1D_BsJpsiSignal          ->Write();
 hDistSign1D_BsOther               ->Write();
 hDistSign1D_Other                 ->Write();
                                
 hPointingAngle_BsJpsiSignal        ->Write();
 hPointingAngle_BsOther             ->Write();
 hPointingAngle_Other               ->Write();
 
 hNEvents_vsGenMuonPt              ->Write();
 hNTrigEvents_vsGenMuonPt          ->Write();
 
 hNEvents_vsGenJPsiPt              ->Write();
 hNTrigEvents_vsGenJPsiPt          ->Write();
 
 hChi2_BsJpsiSignal                 ->Write();
 hChi2Ndof_BsJpsiSignal             ->Write();
 hBsVtxProb_BsJpsiSignal            ->Write();
 
 hChi2_BsOther                      ->Write();
 hChi2Ndof_BsOther                  ->Write();
 hBsVtxProb_BsOther                 ->Write();
 
 hChi2_Other                        ->Write();
 hChi2Ndof_Other                    ->Write();
 hBsVtxProb_Other                   ->Write();

 histofile->Write();
 histofile->Close();

}


//**************
// function to book histograms
void bookHistos(){
  
  // reconstructed jpsi mass before strong cuts
  hJPsiMass_BsJpsiSignal            = new TH1F ("hJPsiMass_BsJpsiSignal", "hJPsiMass_BsJpsiSignal", 100, 2.8, 3.4 );
  hJPsiMass_BsOther                 = new TH1F ("hJPsiMass_BsOther", "hJPsiMass_BsOther", 100, 2.8, 3.4 );
  hJPsiMass_Other                   = new TH1F ("hJPsiMass_Other", "hJPsiMass_Other", 100, 2.8, 3.4 );
  
  // reconstructed phi mass before strong cuts
  hPhiMass_BsJpsiSignal             = new TH1F ("hPhiMass_BsJpsiSignal", "hPhiMass_BsJpsiSignal", 100, 1, 1.035 );
  hPhiMass_BsOther                  = new TH1F ("hPhiMass_BsOther", "hPhiMass_BsOther", 100, 1, 1.035 );
  hPhiMass_Other                    = new TH1F ("hPhiMass_Other", "hPhiMass_Other", 100, 1, 1.035 );
  
  // reconstructed phi mass AFTER strong cuts
  hPhiMassFinal_BsJpsiSignal        = new TH1F ("hPhiMassFinal_BsJpsiSignal", "hPhiMassFinal_BsJpsiSignal", 100, 1, 1.035 );
  hPhiMassFinal_BsOther             = new TH1F ("hPhiMassFinal_BsOther", "hPhiMassFinal_BsOther", 100, 1, 1.035 );
  hPhiMassFinal_Other               = new TH1F ("hPhiMassFinal_Other", "hPhiMassFinal_Other", 100, 1, 1.035 );
  
  // mass of Bs candidate AFTER strong cuts
  hBsMassFinal_BsJpsiSignal         = new TH1F ("hBsMassFinal_BsJpsiSignal", "hBsMassFinal_BsJpsiSignal", 100, 5.2, 5.7 );
  hBsMassFinal_BsOther              = new TH1F ("hBsMassFinal_BsOther", "hBsMassFinal_BsOther", 100,  5.2, 5.7 );
  hBsMassFinal_Other                = new TH1F ("hBsMassFinal_Other", "hBsMassFinal_Other", 100,  5.2, 5.7 );
  
  // mass of fitted Bs AFTER strong cuts
  hBsMassFinalAfterFit_BsJpsiSignal = new TH1F ("hBsMassFinalAfterFit_BsJpsiSignal", "hBsMassFinalAfterFit_BsJpsiSignal", 100,  5.2, 5.7 );
  hBsMassFinalAfterFit_BsOther      = new TH1F ("hBsMassFinalAfterFit_BsOther", "hBsMassFinalAfterFit_BsOther", 100,  5.2, 5.7 );
  hBsMassFinalAfterFit_Other        = new TH1F ("hBsMassFinalAfterFit_Other", "hBsMassFinalAfterFit_Other", 100,  5.2, 5.7 );
  hBsMassFinalAfterFit_BsJpsiKK     = new TH1F ("hBsMassFinalAfterFit_BsJpsiKK", "hBsMassFinalAfterFit_BsJpsiKK", 100,  5.2, 5.7 );
  
  // Kaon pt
  hK1Pt_BsJpsiSignal                = new TH1F("hK1Pt_BsJpsiSignal", "hK1Pt_BsJpsiSignal", 100, 0, 5);
  hK1Pt_BsOther                     = new TH1F("hK1Pt_BsOther", "hK1Pt_BsOther", 100, 0, 5);
  hK1Pt_Other                       = new TH1F("hK1Pt_Other", "hK1Pt_Other", 100, 0, 5);
  
  hK2Pt_BsJpsiSignal                = new TH1F("hK2Pt_BsJpsiSignal", "hK2Pt_BsJpsiSignal", 100, 0, 5);
  hK2Pt_BsOther                     = new TH1F("hK2Pt_BsOther", "hK2Pt_BsOther", 100, 0, 5);
  hK2Pt_Other                       = new TH1F("hK2Pt_Other", "hK2Pt_Other", 100, 0, 5);
  
  // decay length 3D
  hDist3D_BsJpsiSignal              = new TH1F("hDist3D_BsJpsiSignal", "hDist3D_BsJpsiSignal", 100, 0, 1 );
  hDist3D_BsOther                   = new TH1F("hDist3D_BsOther", "hDist3D_BsOther", 100, 0, 1 );
  hDist3D_Other                     = new TH1F("hDist3D_Other", "hDist3D_Other", 100, 0, 1 );
  
  // decay length significance 3D
  hDistSign3D_BsJpsiSignal          = new TH1F("hDistSign3D_BsJpsiSignal", "hDistSign3D_BsJpsiSignal", 100, 0, 50 );
  hDistSign3D_BsOther               = new TH1F("hDistSign3D_BsOther", "hDistSign3D_BsOther", 100, 0, 50 );
  hDistSign3D_Other                 = new TH1F("hDistSign3D_Other", "hDistSign3D_Other", 100, 0, 50 );
  
  // transverse decay length significance
  hDistSign1D_BsJpsiSignal          = new TH1F("hDistSign1D_BsJpsiSignal", "hDistSign1D_BsJpsiSignal", 100, 0, 50 );
  hDistSign1D_BsOther               = new TH1F("hDistSign1D_BsOther", "hDistSign1D_BsOther", 100, 0, 50 );
  hDistSign1D_Other                 = new TH1F("hDistSign1D_Other", "hDistSign1D_Other", 100, 0, 50 );
  
  // pointing angle (cosine angle between decay length and Bs momentum)
  hPointingAngle_BsJpsiSignal        = new TH1F("hPointingAngle_BsJpsiSignal", "hPointingAngle_BsJpsiSignal", 100, -1, 1);
  hPointingAngle_BsOther             = new TH1F("hPointingAngle_BsOther", "hPointingAngle_BsOther", 100, -1, 1);
  hPointingAngle_Other               = new TH1F("hPointingAngle_Other", "hPointingAngle_Other", 100, -1, 1);
  
  // for trigger efficiency plots
  hNEvents_vsGenMuonPt              = new TH1F("hNEvents_vsGenMuonPt", "hNEvents_vsGenMuonPt", 200, 0, 20);
  hNTrigEvents_vsGenMuonPt          = new TH1F("hNTrigEvents_vsGenMuonPt", "hNTrigEvents_vsGenMuonPt", 200, 0 ,20);
  
  hNEvents_vsGenJPsiPt              = new TH1F("hNEvents_vsGenJPsiPt", "hNEvents_vsGenJPsiPt", 200, 0, 40);
  hNTrigEvents_vsGenJPsiPt          = new TH1F("hNTrigEvents_vsGenJPsiPt", "hNTrigEvents_vsGenJPsiPt", 200, 0 ,40);
  
  //kinematic fit chi^2 and fit probability
  hChi2_BsJpsiSignal                 = new TH1F("hChi2_BsJpsiSignal", "hChi2_BsJpsiSignal", 100, 0,30);
  hChi2Ndof_BsJpsiSignal             = new TH1F("hChi2Ndof_BsJpsiSignal", "hChi2Ndof_BsJpsiSignal", 100, 0, 30);
  hBsVtxProb_BsJpsiSignal            = new TH1F("hBsVtxProb_BsJpsiSignal", "hBsVtxProb_BsJpsiSignal", 100, -0.5, 1.5);
  
  hChi2_BsOther                      = new TH1F("hChi2_BsOther", "hChi2_BsOther", 100, 0,30);
  hChi2Ndof_BsOther                  = new TH1F("hChi2Ndof_BsOther", "hChi2Ndof_BsOther", 100, 0, 30);
  hBsVtxProb_BsOther                 = new TH1F("hBsVtxProb_BsOther", "hBsVtxProb_BsOther", 100, -0.5, 1.5);

  hChi2_Other                        = new TH1F("hChi2_Other", "hChi2_Other", 100, 0,30);
  hChi2Ndof_Other                    = new TH1F("hChi2Ndof_Other", "hChi2Ndof_Other", 100, 0, 30);
  hBsVtxProb_Other                   = new TH1F("hBsVtxProb_Other", "hBsVtxProb_Other", 100, -0.5, 1.5);
  
}


