/*
 *  See header file for a description of this class.
 *
 *  $Date: 2009/10/08 10:14:36 $
 *  $Revision: 1.25 $
 *  \author F. Chlebana - Fermilab
 *          K. Hatakeyama - Rockefeller University
 */

#include "DQMOffline/JetMET/interface/JetMETAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
using namespace std;
using namespace edm;

#define DEBUG 0

// ***********************************************************
JetMETAnalyzer::JetMETAnalyzer(const edm::ParameterSet& pSet) {

  parameters = pSet;
  
  // Jet Collection Label
  theSCJetCollectionLabel       = parameters.getParameter<edm::InputTag>("SCJetsCollectionLabel");
  theICJetCollectionLabel       = parameters.getParameter<edm::InputTag>("ICJetsCollectionLabel");
  theJPTJetCollectionLabel      = parameters.getParameter<edm::InputTag>("JPTJetsCollectionLabel");
  thePFJetCollectionLabel       = parameters.getParameter<edm::InputTag>("PFJetsCollectionLabel");

  theTriggerResultsLabel        = parameters.getParameter<edm::InputTag>("TriggerResultsLabel");
  
  theJetAnalyzerFlag            = parameters.getUntrackedParameter<bool>("DoJetAnalysis",    true);
  theJetCleaningFlag            = parameters.getUntrackedParameter<bool>("DoJetCleaning",    true);
  theJetPtAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoJetPtAnalysis",  false);
  theJPTJetAnalyzerFlag         = parameters.getUntrackedParameter<bool>("DoJPTJetAnalysis", true);
  thePFJetAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoPFJetAnalysis",  true);
  theCaloMETAnalyzerFlag        = parameters.getUntrackedParameter<bool>("DoCaloMETAnalysis",true);
  theTcMETAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoTcMETAnalysis",  true);
  thePfMETAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoPfMETAnalysis",  true);
  theHTMHTAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoHTMHTAnalysis",  true);

  // --- do the analysis on the Jets
  if(theJetAnalyzerFlag) {
    theSCJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    theSCJetAnalyzer->setSource("SISConeJets");
    theICJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    theICJetAnalyzer->setSource("IterativeConeJets");  
  }
  if(theJetCleaningFlag) {
    theCleanedSCJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("CleanedjetAnalysis"));
    theCleanedSCJetAnalyzer->setSource("CleanedSISConeJets");
    theCleanedICJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("CleanedjetAnalysis"));
    theCleanedICJetAnalyzer->setSource("CleanedIterativeConeJets");
  }
  // Do Pt analysis
  if(theJetPtAnalyzerFlag ) {
    thePtSCJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("PtAnalysis"));
    thePtSCJetAnalyzer->setSource("PtAnalysisSISConeJets");
    thePtICJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("PtAnalysis"));
    thePtICJetAnalyzer->setSource("PtAnalysisIterativeConeJets");
  }

  // --- do the analysis on JPT Jets
  if(theJPTJetAnalyzerFlag) {
    theJPTJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("JPTJetAnalysis"));
    theJPTJetAnalyzer->setSource("JPTJets");
  }

  // --- do the analysis on the PFJets
  if(thePFJetAnalyzerFlag)
    thePFJetAnalyzer = new PFJetAnalyzer(parameters.getParameter<ParameterSet>("pfJetAnalysis"));

  LoJetTrigger = parameters.getParameter<std::string>("JetLo");
  HiJetTrigger = parameters.getParameter<std::string>("JetHi");

  // --- do the analysis on the MET
  if(theCaloMETAnalyzerFlag){
     theCaloMETAnalyzer       = new CaloMETAnalyzer(parameters.getParameter<ParameterSet>("caloMETAnalysis"));
     theCaloMETNoHFAnalyzer   = new CaloMETAnalyzer(parameters.getParameter<ParameterSet>("caloMETNoHFAnalysis"));
     theCaloMETHOAnalyzer     = new CaloMETAnalyzer(parameters.getParameter<ParameterSet>("caloMETHOAnalysis"));
     theCaloMETNoHFHOAnalyzer = new CaloMETAnalyzer(parameters.getParameter<ParameterSet>("caloMETNoHFHOAnalysis"));
  }
  if(theTcMETAnalyzerFlag){
     theTcMETAnalyzer = new METAnalyzer(parameters.getParameter<ParameterSet>("tcMETAnalysis"));
  }
  if(thePfMETAnalyzerFlag){
     thePfMETAnalyzer = new PFMETAnalyzer(parameters.getParameter<ParameterSet>("pfMETAnalysis"));
  }
  if(theHTMHTAnalyzerFlag){
     theHTMHTAnalyzer         = new HTMHTAnalyzer(parameters.getParameter<ParameterSet>("HTMHTAnalysis"));
  }

  processname_ = parameters.getParameter<std::string>("processname");

}

// ***********************************************************
JetMETAnalyzer::~JetMETAnalyzer() {   

  if(theJetAnalyzerFlag) {
    delete theSCJetAnalyzer;
    delete theICJetAnalyzer;
  }
  if(theJetCleaningFlag) {
    delete theCleanedSCJetAnalyzer;
    delete theCleanedICJetAnalyzer;
  }
  if(theJetPtAnalyzerFlag) {
    delete thePtSCJetAnalyzer;
    delete thePtICJetAnalyzer;
  }

  if(theJPTJetAnalyzerFlag)        delete theJPTJetAnalyzer;
  if(thePFJetAnalyzerFlag)         delete thePFJetAnalyzer;

  if(theCaloMETAnalyzerFlag){
    delete theCaloMETAnalyzer;
    delete theCaloMETNoHFAnalyzer;
    delete theCaloMETHOAnalyzer;
    delete theCaloMETNoHFHOAnalyzer;
  }
  if(theTcMETAnalyzerFlag)         delete theTcMETAnalyzer;
  if(thePfMETAnalyzerFlag)         delete thePfMETAnalyzer;
  if(theHTMHTAnalyzerFlag)         delete theHTMHTAnalyzer;

}

// ***********************************************************
void JetMETAnalyzer::beginJob(edm::EventSetup const& iSetup) {

  metname = "JetMETAnalyzer";

  LogTrace(metname)<<"[JetMETAnalyzer] Parameters initialization";
  dbe = edm::Service<DQMStore>().operator->();

  //
  //--- Jet
  if(theJetAnalyzerFlag) { 
    theSCJetAnalyzer->beginJob(iSetup, dbe);
    theICJetAnalyzer->beginJob(iSetup, dbe);
  }
  if(theJetCleaningFlag) {
    theCleanedSCJetAnalyzer->beginJob(iSetup, dbe);
    theCleanedICJetAnalyzer->beginJob(iSetup, dbe);
  }
  if(theJetPtAnalyzerFlag ) {
    thePtSCJetAnalyzer->beginJob(iSetup, dbe);
    thePtICJetAnalyzer->beginJob(iSetup, dbe);
  }

  if(theJPTJetAnalyzerFlag) theJPTJetAnalyzer->beginJob(iSetup, dbe);
  if(thePFJetAnalyzerFlag)  thePFJetAnalyzer->beginJob(iSetup, dbe);

  //
  //--- MET
  if(theCaloMETAnalyzerFlag){
    theCaloMETAnalyzer->beginJob(iSetup, dbe);
    theCaloMETNoHFAnalyzer->beginJob(iSetup, dbe);
    theCaloMETHOAnalyzer->beginJob(iSetup, dbe);
    theCaloMETNoHFHOAnalyzer->beginJob(iSetup, dbe);
  }
  if(theTcMETAnalyzerFlag) theTcMETAnalyzer->beginJob(iSetup, dbe);
  if(thePfMETAnalyzerFlag) thePfMETAnalyzer->beginJob(iSetup, dbe);
  if(theHTMHTAnalyzerFlag) theHTMHTAnalyzer->beginJob(iSetup, dbe);
  
  dbe->setCurrentFolder("JetMET");
  lumisecME = dbe->book1D("lumisec", "lumisec", 500, 0., 500.);

}

// ***********************************************************
void JetMETAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  //LogDebug("JetMETAnalyzer") << "beginRun, run " << run.id();

  //
  //--- htlConfig_
  //processname_="HLT";
  hltConfig_.init(processname_);
  if (!hltConfig_.init(processname_)) {
    processname_ = "FU";
    if (!hltConfig_.init(processname_)){
      LogDebug("JetMETAnalyzer") << "HLTConfigProvider failed to initialize.";
    }
  }

  if (hltConfig_.size()){
    dbe->setCurrentFolder("JetMET");
    hltpathME = dbe->book1D("hltpath", "hltpath", 300, 0., 300.);
  }

  for (unsigned int j=0; j!=hltConfig_.size(); ++j) {
    hltpathME->setBinLabel(j+1,hltConfig_.triggerName(j));
  }

  //
  //--- Jet

  //
  //--- MET
  if(theCaloMETAnalyzerFlag){
    theCaloMETAnalyzer->beginRun(iRun, iSetup);
    theCaloMETNoHFAnalyzer->beginRun(iRun, iSetup);
    theCaloMETHOAnalyzer->beginRun(iRun, iSetup);
    theCaloMETNoHFHOAnalyzer->beginRun(iRun, iSetup);
  }
  //if(theTcMETAnalyzerFlag) theTcMETAnalyzer->beginRun(iRun, iSetup);
  //if(thePfMETAnalyzerFlag) thePfMETAnalyzer->beginRun(iRun, iSetup);
  //if(theHTMHTAnalyzerFlag) theHTMHTAnalyzer->beginRun(iRun, iSetup);

}

// ***********************************************************
void JetMETAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

  //
  //--- Jet

  //
  //--- MET
  if(theCaloMETAnalyzerFlag){
    theCaloMETAnalyzer->endRun(iRun, iSetup, dbe);
    theCaloMETNoHFAnalyzer->endRun(iRun, iSetup, dbe);
    theCaloMETHOAnalyzer->endRun(iRun, iSetup, dbe);
    theCaloMETNoHFHOAnalyzer->endRun(iRun, iSetup, dbe);
  }

}

// ***********************************************************
void JetMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  LogTrace(metname)<<"[JetMETAnalyzer] Analysis of event # ";

  // *** Fill lumisection ME
  int myLuminosityBlock;
  myLuminosityBlock = iEvent.luminosityBlock();
  lumisecME->Fill(myLuminosityBlock);

  // **** Get the TriggerResults container
  edm::Handle<TriggerResults> triggerResults;
  iEvent.getByLabel(theTriggerResultsLabel, triggerResults);

  // *** Fill trigger results ME
  if (&triggerResults){
    edm::TriggerNames triggerNames;
    triggerNames.init(*(triggerResults.product()));
    for (unsigned int j=0; j!=hltConfig_.size(); ++j) {
      if (triggerResults->accept(j)){
        hltpathME->Fill(j);
      }
    }
  }

  if (DEBUG)  std::cout << "trigger label " << theTriggerResultsLabel << std::endl;

  Int_t JetLoPass = 0;
  Int_t JetHiPass = 0;

  if (triggerResults.isValid()) {

    if (DEBUG) std::cout << "trigger valid " << std::endl;
    edm::TriggerNames triggerNames;    // TriggerNames class
    triggerNames.init(*triggerResults);
    unsigned int n = triggerResults->size();
    for (unsigned int i=0; i!=n; i++) {

      if ( triggerNames.triggerName(i) == LoJetTrigger ) {
	JetLoPass =  triggerResults->accept(i);
	if (DEBUG) std::cout << "Found  HLT_Jet30" << std::endl;
      }
      if ( triggerNames.triggerName(i) == HiJetTrigger ) {
	JetHiPass =  triggerResults->accept(i);
      }
    }

  } else {

    //
    triggerResults = edm::Handle<TriggerResults>(); 

    if (DEBUG) std::cout << "trigger not valid " << std::endl;
    edm::LogInfo("JetMETAnalyzer") << "TriggerResults::HLT not found, "
      "automatically select events";

  }

  if (DEBUG) {
    std::cout << ">>> Trigger  Lo = " <<  JetLoPass
	      <<             " Hi = " <<  JetHiPass
	      << std::endl;
  }


  // **** Get the Calo Jet container
  edm::Handle<reco::CaloJetCollection> caloJets;

  // **** Get the SISCone Jet container
  iEvent.getByLabel(theSCJetCollectionLabel, caloJets);

  if(caloJets.isValid()){
    theCleanedSCJetAnalyzer->setJetHiPass(JetHiPass);
    theCleanedSCJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if( theJetAnalyzerFlag  ){ // KH Is this correct?
        LogTrace(metname)<<"[JetMETAnalyzer] Call to the SC Jet analyzer"; // KH Is this correct?
        if (cal == caloJets->begin()) {
          theCleanedSCJetAnalyzer->setNJets(caloJets->size());
          theCleanedSCJetAnalyzer->setLeadJetFlag(1);
        } else {
          theCleanedSCJetAnalyzer->setLeadJetFlag(0);
        }
        theCleanedSCJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }
  if(caloJets.isValid()){
    theSCJetAnalyzer->setJetHiPass(JetHiPass);
    theSCJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetCleaningFlag){ // KH Is this correct?
        LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned SC Jet analyzer"; // KH Is this correct?
        if (cal == caloJets->begin()) {
          theSCJetAnalyzer->setNJets(caloJets->size());
          theSCJetAnalyzer->setLeadJetFlag(1);
        } else {
          theSCJetAnalyzer->setLeadJetFlag(0);
        }
        theSCJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }
  if(caloJets.isValid()){
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetPtAnalyzerFlag){
        LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt SisCone analyzer";
        if (cal == caloJets->begin()) {
          thePtSCJetAnalyzer->setNJets(caloJets->size());
        }
        thePtSCJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }

  // **** Get the Iterative Cone Jet container
  iEvent.getByLabel(theICJetCollectionLabel, caloJets);

  if(caloJets.isValid()){
    theICJetAnalyzer->setJetHiPass(JetHiPass);
    theICJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the IC Jet analyzer";
	if (cal == caloJets->begin()) {
	  theICJetAnalyzer->setNJets(caloJets->size());
	  theICJetAnalyzer->setLeadJetFlag(1);
	} else {
	  theICJetAnalyzer->setLeadJetFlag(0);
	}
	//	theICJetAnalyzer->analyze(iEvent, iSetup, *triggerResults, *cal);	
	theICJetAnalyzer->analyze(iEvent, iSetup, *cal);	
      }
    }
  }
  if(caloJets.isValid()){
    theCleanedICJetAnalyzer->setJetHiPass(JetHiPass);
    theCleanedICJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetCleaningFlag){
        LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned IC Jet analyzer";
        if (cal == caloJets->begin()) {
          theCleanedICJetAnalyzer->setNJets(caloJets->size());
          theCleanedICJetAnalyzer->setLeadJetFlag(1);
        } else {
          theCleanedICJetAnalyzer->setLeadJetFlag(0);
        }
        //      theICJetAnalyzer->analyze(iEvent, iSetup, *triggerResults, *cal);
        theCleanedICJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }
  if(caloJets.isValid()){
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetPtAnalyzerFlag){
        LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt ICone analyzer";
        if (cal == caloJets->begin()) {
          thePtICJetAnalyzer->setNJets(caloJets->size());
        }
        thePtICJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }


// **** Get the JPT Jet container
  iEvent.getByLabel(theJPTJetCollectionLabel, caloJets);
  //jpt
  if(caloJets.isValid()){
    theJPTJetAnalyzer->setJetHiPass(JetHiPass);
    theJPTJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJPTJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the JPT Jet analyzer";
	if (cal == caloJets->begin()) {	  
	  theJPTJetAnalyzer->setNJets(caloJets->size());
	  theJPTJetAnalyzer->setLeadJetFlag(1);
	} else {
	  theJPTJetAnalyzer->setLeadJetFlag(0);
	}
	theJPTJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }

  // **** Get the PFlow Jet container
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel(thePFJetCollectionLabel, pfJets);

  if(pfJets.isValid()){
//     for (reco::PFJetCollection::const_iterator cal = pfJets->begin(); cal!=pfJets->end(); ++cal){
//       if(thePFJetAnalyzerFlag){
// 	LogTrace(metname)<<"[JetMETAnalyzer] Call to the PFJet analyzer";
// 	thePFJetAnalyzer->analyze(iEvent, iSetup, *cal);
//       }
//     }
    thePFJetAnalyzer->setJetHiPass(JetHiPass);
    thePFJetAnalyzer->setJetLoPass(JetLoPass);
    for (reco::PFJetCollection::const_iterator cal = pfJets->begin(); cal!=pfJets->end(); ++cal){
      if(thePFJetAnalyzerFlag){
	if (cal == pfJets->begin()) {	  
	  thePFJetAnalyzer->setLeadJetFlag(1);
	} else {
	  thePFJetAnalyzer->setLeadJetFlag(0);
	}
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the PFJet analyzer";
	thePFJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  } else {
    if (DEBUG) LogTrace(metname)<<"[JetMETAnalyzer] pfjets NOT VALID!!";
  }

  //
  // **** CaloMETAnalyzer **** //
  //
  if(theCaloMETAnalyzerFlag){

    theCaloMETAnalyzer->analyze(iEvent,       iSetup, *triggerResults);
    theCaloMETNoHFAnalyzer->analyze(iEvent,   iSetup, *triggerResults);
    theCaloMETHOAnalyzer->analyze(iEvent,     iSetup, *triggerResults);
    theCaloMETNoHFHOAnalyzer->analyze(iEvent, iSetup, *triggerResults);

  }

  //
  // **** TcMETAnalyzer **** //
  //
  if(theTcMETAnalyzerFlag){

    theTcMETAnalyzer->analyze(iEvent, iSetup, *triggerResults);

  }

  //
  // **** PfMETAnalyzer **** //
  //
  if(thePfMETAnalyzerFlag){
    
    thePfMETAnalyzer->analyze(iEvent, iSetup, *triggerResults);

  }

  //
  // **** HTMHTAnalyzer **** //
  //
  if(theHTMHTAnalyzerFlag){

    theHTMHTAnalyzer->analyze(iEvent, iSetup, *triggerResults);

  }

}

// ***********************************************************
void JetMETAnalyzer::endJob(void) {
  LogTrace(metname)<<"[JetMETAnalyzer] Saving the histos";
  bool outputMEsInRootFile   = parameters.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = parameters.getParameter<std::string>("OutputFileName");

  if(outputMEsInRootFile){
    dbe->showDirStructure();
    dbe->save(outputFileName);
  }

  if(theCaloMETAnalyzerFlag) theCaloMETAnalyzer->endJob();
  //if(theTcMETAnalyzerFlag) theTcMETAnalyzer->endJob();
  //if(thePfMETAnalyzerFlag) thePfMETAnalyzer->endJob();
  //if(theHTMHTAnalyzerFlag) theHTMHTAnalyzer->endJob();

}

