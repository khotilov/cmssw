// $Id: TtHadEvtSolutionMaker.cc,v 1.4 2007/10/31 17:22:12 mfhansen Exp $

#include "TopQuarkAnalysis/TopEventProducers/interface/TtHadEvtSolutionMaker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/Utilities/interface/DeltaR.h"
#include "AnalysisDataFormats/TopObjects/interface/TtHadEvtSolution.h"
#include "TopQuarkAnalysis/TopTools/interface/JetPartonMatching.h"
#include "TopQuarkAnalysis/TopKinFitter/interface/TtHadKinFitter.h"
#include "TopQuarkAnalysis/TopJetCombination/interface/TtHadSimpleBestJetComb.h"
#include "TopQuarkAnalysis/TopJetCombination/interface/TtHadLRJetCombObservables.h"
#include "TopQuarkAnalysis/TopJetCombination/interface/TtHadLRJetCombCalc.h"
#include "TopQuarkAnalysis/TopEventSelection/interface/TtHadLRSignalSelObservables.h"
#include "TopQuarkAnalysis/TopEventSelection/interface/TtHadLRSignalSelCalc.h"

#include <memory>


/// constructor
TtHadEvtSolutionMaker::TtHadEvtSolutionMaker(const edm::ParameterSet & iConfig) {
  // configurables
  lJetSrc_         = iConfig.getParameter<edm::InputTag>    ("lJetSource");
  bJetSrc_         = iConfig.getParameter<edm::InputTag>    ("bJetSource");
  doKinFit_        = iConfig.getParameter<bool>             ("doKinFit");
  addLRSignalSel_  = iConfig.getParameter<bool>             ("addLRSignalSel");
  lrSignalSelObs_  = iConfig.getParameter<std::vector<int> >("lrSignalSelObs");
  lrSignalSelFile_ = iConfig.getParameter<std::string>      ("lrSignalSelFile");
  addLRJetComb_    = iConfig.getParameter<bool>             ("addLRJetComb");
  lrJetCombObs_    = iConfig.getParameter<std::vector<int> >("lrJetCombObs");
  lrJetCombFile_   = iConfig.getParameter<std::string>      ("lrJetCombFile");
  maxNrIter_       = iConfig.getParameter<int>              ("maxNrIter");
  maxDeltaS_       = iConfig.getParameter<double>           ("maxDeltaS");
  maxF_            = iConfig.getParameter<double>           ("maxF");
  jetParam_        = iConfig.getParameter<int>              ("jetParametrisation");
  constraints_     = iConfig.getParameter<std::vector<int> >("constraints");
  matchToGenEvt_   = iConfig.getParameter<bool>             ("matchToGenEvt");

  // define kinfitter
  if(doKinFit_){
    myKinFitter = new TtHadKinFitter(jetParam_, maxNrIter_, maxDeltaS_, maxF_, constraints_);
  }
  
  
  // define jet combinations related calculators
  mySimpleBestJetComb                    = new TtHadSimpleBestJetComb();
  myLRSignalSelObservables               = new TtHadLRSignalSelObservables();
  myLRJetCombObservables                 = new TtHadLRJetCombObservables();
  if (addLRJetComb_)   myLRJetCombCalc   = new TtHadLRJetCombCalc(lrJetCombFile_, lrJetCombObs_);
  
  // instantiate signal selection calculator
  if (addLRSignalSel_) myLRSignalSelCalc = new TtHadLRSignalSelCalc(lrSignalSelFile_, lrSignalSelObs_);
 
  // define what will be produced
  produces<std::vector<TtHadEvtSolution> >();
  
}


/// destructor
TtHadEvtSolutionMaker::~TtHadEvtSolutionMaker()
{
  if (doKinFit_) {
    delete myKinFitter;
  }
  delete mySimpleBestJetComb;
  delete myLRSignalSelObservables;
  delete myLRJetCombObservables;
  if(addLRSignalSel_) delete myLRSignalSelCalc;
  if(addLRJetComb_)   delete myLRJetCombCalc;
}


void TtHadEvtSolutionMaker::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  // TopObject Selection
  // Select Jets (TopJet vector is sorted on recET, so four first elements in both the lJets and bJets vector are the same )

  bool jetsFound = false;
  edm::Handle<std::vector<TopJet> > lJets;
  iEvent.getByLabel(lJetSrc_, lJets);
  edm::Handle<std::vector<TopJet> > bJets;
  iEvent.getByLabel(bJetSrc_, bJets);

  if (lJets->size() >= 6) jetsFound = true;
  
  // Build Event solutions according to the ambiguity in the jet combination
  // Note, hardcoded to only run through the 6 most energetic jets - could be changed ....
  
  std::vector<TtHadEvtSolution> * evtsols = new std::vector<TtHadEvtSolution>();
  if(jetsFound){
    //std::cout<<"constructing solutions"<<std::endl;
    for (unsigned int p=0; p<3; p++) {  // loop over light jet p
      for (unsigned int q=p+1; q<4; q++) { // loop over light jet q
	for (unsigned int j=q+1; j<5; j++) {  // loop over light jet j
	  for (unsigned int k=j+1; k<6; k++) { // loop over light jet k
	    for (unsigned int bh=0; bh!=bJets->size(); bh++) { //loop over hadronic b-jet1
	      if(!(bh==p || bh==q || bh==j || bh==k)) { 
		for (unsigned int bbarh=0; bbarh!=bJets->size(); bbarh++) { //loop over hadronic b-jet2
		  if (!(bbarh==p || bbarh==q || bbarh==j || bbarh==k) && bbarh==bh) {
		    // Make event solutions for all possible combinations of the 4 light
		    // jets and 2 possible b-jets, not including the option of the b's being swapped. 
		    // Hadp,Hadq is one pair, Hadj,Hadk the other
		    vector<TtHadEvtSolution> asol;
		    asol.resize(3);
		    //[p][q][b] and [j][k][bbar]
		    asol[0].setHadp(lJets, p);
		    asol[0].setHadq(lJets, q);
		    asol[0].setHadj(lJets, j);
		    asol[0].setHadk(lJets, k);
		    asol[0].setHadb(bJets, bh);
		    asol[0].setHadbbar(bJets, bbarh);

		    //[p][j][b] and [q][k][bbar]
		    asol[1].setHadp(lJets, p);
		    asol[1].setHadq(lJets, j);
		    asol[1].setHadj(lJets, q);
		    asol[1].setHadk(lJets, k);
		    asol[1].setHadb(bJets, bh);
		    asol[1].setHadbbar(bJets, bbarh);

		    //[p][k][b] and [j][q][bbar]
		    asol[2].setHadp(lJets, p);
		    asol[2].setHadq(lJets, k);
		    asol[2].setHadj(lJets, j);
		    asol[2].setHadk(lJets, q);
		    asol[2].setHadb(bJets, bh);
		    asol[2].setHadbbar(bJets, bbarh);
		   
		    if(doKinFit_){
		      for(unsigned int i=0;i!=asol.size();i++){
			asol[i] = myKinFitter->addKinFitInfo(&(asol[i]));
			asol[i].setJetParametrisation(jetParam_);
		      }
		     
		    }else{
		      std::cout<<"Fitting needed to decide on best solution, enable fitting!"<<std::endl;
		    }
		    // these lines calculate the observables to be used in the TtHadSignalSelection LR

		    for(unsigned int i=0;i!=asol.size();i++){ 
		    (*myLRSignalSelObservables)(asol[i]);
		    }
		    // if asked for, calculate with these observable values the LRvalue and 
		    // (depending on the configuration) probability this event is signal
		    if(addLRSignalSel_){
		      for(unsigned int i=0;i!=asol.size();i++){
			(*myLRSignalSelCalc)(asol[i]);
		      }
		    }
		    
		    // these lines calculate the observables to be used in the TtHadJetCombination LR
		    for(unsigned int i=0;i!=asol.size();i++){
		    (*myLRJetCombObservables)(asol[i]);
		    } 
		    // if asked for, calculate with these observable values the LRvalue and 
		    // (depending on the configuration) probability a jet combination is correct
		    if(addLRJetComb_){
		      for(unsigned int i=0;i!=asol.size();i++){
			(*myLRJetCombCalc)(asol[i]);
		      }
		    }
		    //std::cout<<"SignalSelLRval = "<<asol.getLRSignalEvtLRval()<<"  JetCombProb = "<<asol.getLRSignalEvtProb()<<std::endl;
		    //std::cout<<"JetCombLRval = "<<asol.getLRJetCombLRval()<<"  JetCombProb = "<<asol.getLRJetCombProb()<<std::endl;
		    
		    // fill solution to vector with all possible solutions 
		    for(unsigned int i=0;i!=asol.size();i++){
		      evtsols->push_back(asol[i]);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    
    // add TtHadSimpleBestJetComb to solutions
    int simpleBestJetComb = (*mySimpleBestJetComb)(*evtsols);

    for(size_t s=0; s<evtsols->size(); s++){
      (*evtsols)[s].setSimpleBestJetComb(simpleBestJetComb);    
      // if asked for, match the event solutions to the gen Event
      if(matchToGenEvt_){
	int bestSolution = -999; 
	int bestSolutionChangeW1Q = -999;
	int bestSolutionChangeW2Q = -999;
	edm::Handle<TtGenEvent> genEvt;
	iEvent.getByLabel ("genEvt",genEvt); 
	vector<const reco::Candidate*> quarks;
	const reco::Candidate & genp  = *(genEvt->quarkFromTop());
	const reco::Candidate & genq  = *(genEvt->quarkFromTopBar());
	const reco::Candidate & genb  = *(genEvt->b());
	const reco::Candidate & genj  = *(genEvt->quarkFromAntiTop());
	const reco::Candidate & genk  = *(genEvt->quarkFromAntiTopBar());
	const reco::Candidate & genbbar = *(genEvt->bBar());

	quarks.push_back( &genp );       
	quarks.push_back( &genq );   
	quarks.push_back( &genb );
	quarks.push_back( &genj );       
	quarks.push_back( &genk );   
	quarks.push_back( &genbbar );
	vector<const reco::Candidate*> jets;         
	for(size_t s=0; s<evtsols->size(); s++) {
	  jets.clear();     
	  const reco::Candidate & jetp  = (*evtsols)[s].getRecHadp();
	  const reco::Candidate & jetq  = (*evtsols)[s].getRecHadq();
	  const reco::Candidate & jetbh = (*evtsols)[s].getRecHadb();
	  const reco::Candidate & jetj  = (*evtsols)[s].getRecHadj();
	  const reco::Candidate & jetk  = (*evtsols)[s].getRecHadk();
	  const reco::Candidate & jetbbar = (*evtsols)[s].getRecHadbbar();
	  jets.push_back( &jetp );      
	  jets.push_back( &jetq );        
	  jets.push_back( &jetbh );
	  jets.push_back( &jetj );
	  jets.push_back( &jetk );
	  jets.push_back( &jetbbar );
	  JetPartonMatching aMatch(jets,quarks,2);  // 1: SpaceAngle; 2: DeltaR  
	  (*evtsols)[s].setGenEvt(genEvt);   
	  (*evtsols)[s].setMCBestSumAngles(aMatch.getSumAngles());
	  (*evtsols)[s].setMCBestAngleHadp(aMatch.getAngleForParton(0));
	  (*evtsols)[s].setMCBestAngleHadq(aMatch.getAngleForParton(1));
	  (*evtsols)[s].setMCBestAngleHadb(aMatch.getAngleForParton(2));
	  (*evtsols)[s].setMCBestAngleHadb(aMatch.getAngleForParton(2));
	  (*evtsols)[s].setMCBestAngleHadj(aMatch.getAngleForParton(3));
	  (*evtsols)[s].setMCBestAngleHadk(aMatch.getAngleForParton(4));
	  (*evtsols)[s].setMCBestAngleHadbbar(aMatch.getAngleForParton(5));
	  
	  // Check match - checking if two light quarks are swapped wrt matched gen particle
	  if((aMatch.getMatchForParton(2) == 2 && aMatch.getMatchForParton(5) == 5)
	     || (aMatch.getMatchForParton(2) == 5 && aMatch.getMatchForParton(5) == 2)){ // check b-jets
	    
	    if(aMatch.getMatchForParton(3) == 3 && aMatch.getMatchForParton(4) == 4){ //check light jets
	      bestSolutionChangeW2Q = 0;
	      if(aMatch.getMatchForParton(0) == 0 && aMatch.getMatchForParton(1) == 1) { 
		bestSolution = s;
		bestSolutionChangeW1Q = 0;
	      }else{
		if(aMatch.getMatchForParton(0) == 1 && aMatch.getMatchForParton(1) == 0){
		  bestSolution = s;
		  bestSolutionChangeW1Q = 1;
		}
	      }
	    }else{
	      if(aMatch.getMatchForParton(2) == 3 && aMatch.getMatchForParton(3) == 2){ // or check if swapped 
		bestSolutionChangeW2Q = 1;
		if(aMatch.getMatchForParton(0) == 1 && aMatch.getMatchForParton(1) == 0){
		  bestSolution = s;
		  bestSolutionChangeW1Q = 1;
		}else{
		  if(aMatch.getMatchForParton(0) == 0 && aMatch.getMatchForParton(1) == 1) { 
		    bestSolution = s;
		    bestSolutionChangeW1Q = 0;
		  }
		}
	      }
	      if(aMatch.getMatchForParton(2) == 2 && aMatch.getMatchForParton(3) == 3){
		bestSolutionChangeW2Q = 0;
		if(aMatch.getMatchForParton(0) == 0 && aMatch.getMatchForParton(1) == 1) {
		  bestSolution = s; 
		  bestSolutionChangeW1Q = 0;
		} else if(aMatch.getMatchForParton(0) == 1 && aMatch.getMatchForParton(1) == 0) {
		  bestSolution = s;             
		  bestSolutionChangeW1Q = 1;           
		}
	      }
	    }
	  }
	  for(size_t s=0; s<evtsols->size(); s++) {
	    (*evtsols)[s].setMCBestJetComb(bestSolution);
	    (*evtsols)[s].setMCChangeW1Q(bestSolutionChangeW1Q);
	    (*evtsols)[s].setMCChangeW2Q(bestSolutionChangeW2Q);
	  }
	}
      } // end matchEvt
    }       
    //store the vector of solutions to the event  
 
    std::auto_ptr<std::vector<TtHadEvtSolution> > pOut(evtsols);
    iEvent.put(pOut);
  }else {     //end loop jet/MET found
    std::cout<<"No calibrated solutions built, because only "<<lJets->size()<<" were present";
     
    std::auto_ptr<std::vector<TtHadEvtSolution> > pOut(evtsols);
    iEvent.put(pOut);
  }
}  

