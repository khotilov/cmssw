/*
 *  See header file for a description of this class.
 *
 *  $Date: 2010/03/22 09:30:02 $
 *  $Revision: 1.57 $
 *  \author F. Chlebana - Fermilab
 *          K. Hatakeyama - Rockefeller University
 */

#include "DQMOffline/JetMET/interface/JetMETAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include <string>
using namespace std;
using namespace edm;

#define DEBUG 0

// ***********************************************************
JetMETAnalyzer::JetMETAnalyzer(const edm::ParameterSet& pSet) {

  parameters = pSet;
  
  // Jet Collection Label 
  theAKJetCollectionLabel       = parameters.getParameter<edm::InputTag>("AKJetsCollectionLabel");
  theSCJetCollectionLabel       = parameters.getParameter<edm::InputTag>("SCJetsCollectionLabel");
  theICJetCollectionLabel       = parameters.getParameter<edm::InputTag>("ICJetsCollectionLabel");
  theJPTJetCollectionLabel      = parameters.getParameter<edm::InputTag>("JPTJetsCollectionLabel");
  thePFJetCollectionLabel       = parameters.getParameter<edm::InputTag>("PFJetsCollectionLabel");

  theTriggerResultsLabel        = parameters.getParameter<edm::InputTag>("TriggerResultsLabel");
  
  theJetAnalyzerFlag            = parameters.getUntrackedParameter<bool>("DoJetAnalysis",    true);
  theJetCleaningFlag            = parameters.getUntrackedParameter<bool>("DoJetCleaning",    true);
  theIConeJetAnalyzerFlag       = parameters.getUntrackedParameter<bool>("DoIterativeCone",  true);
  theSConeJetAnalyzerFlag       = parameters.getUntrackedParameter<bool>("DoSisCone",  false);
  theJetPtAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoJetPtAnalysis",  false);
  theJetPtCleaningFlag          = parameters.getUntrackedParameter<bool>("DoJetPtCleaning",  false);
  theJPTJetAnalyzerFlag         = parameters.getUntrackedParameter<bool>("DoJPTJetAnalysis", true);
  theJPTJetCleaningFlag         = parameters.getUntrackedParameter<bool>("DoJPTJetCleaning", true);
  thePFJetAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoPFJetAnalysis",  true);
  thePFJetCleaningFlag          = parameters.getUntrackedParameter<bool>("DoPFJetCleaning",  true);
  //
  theDiJetSelectionFlag         = parameters.getUntrackedParameter<bool>("DoDiJetSelection", true);
  //
  theCaloMETAnalyzerFlag        = parameters.getUntrackedParameter<bool>("DoCaloMETAnalysis",true);
  theTcMETAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoTcMETAnalysis",  true);
  theMuCorrMETAnalyzerFlag      = parameters.getUntrackedParameter<bool>("DoMuCorrMETAnalysis",  true);
  thePfMETAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoPfMETAnalysis",  true);
  theHTMHTAnalyzerFlag          = parameters.getUntrackedParameter<bool>("DoHTMHTAnalysis",  true);

  // --- do the analysis on the Jets
  if(theJetAnalyzerFlag) {
    theAKJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    theAKJetAnalyzer->setSource("AntiKtJets");    
    if(theSConeJetAnalyzerFlag) {
      theSCJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
      theSCJetAnalyzer->setSource("SISConeJets");
    }
    if(theIConeJetAnalyzerFlag){
      theICJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
      theICJetAnalyzer->setSource("IterativeConeJets");  
    }
  }

  if(theJetCleaningFlag) {
    theCleanedAKJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("CleanedjetAnalysis"));
    theCleanedAKJetAnalyzer->setSource("CleanedAntiKtJets");
    if(theSConeJetAnalyzerFlag){
      theCleanedSCJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("CleanedjetAnalysis"));
      theCleanedSCJetAnalyzer->setSource("CleanedSISConeJets");
    }
    if(theIConeJetAnalyzerFlag){
      theCleanedICJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("CleanedjetAnalysis"));
      theCleanedICJetAnalyzer->setSource("CleanedIterativeConeJets");
    }
  }

  if(theDiJetSelectionFlag){
    theDiJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("DijetAnalysis"));
    theDiJetAnalyzer->setSource("DiJets");
  }

  // Do Pt analysis
  if(theJetPtAnalyzerFlag ) {
    thePtAKJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("PtAnalysis"));
    thePtAKJetAnalyzer->setSource("PtAnalysisAntiKtJets");
    if(theSConeJetAnalyzerFlag){
      thePtSCJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("PtAnalysis"));
      thePtSCJetAnalyzer->setSource("PtAnalysisSISConeJets");
    }
    if(theIConeJetAnalyzerFlag){
      thePtICJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("PtAnalysis"));
      thePtICJetAnalyzer->setSource("PtAnalysisIterativeConeJets");
    }
  }
  // Do Cleaned Pt analysis
  if(theJetPtCleaningFlag ) {
    theCleanedPtAKJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("CleanedPtAnalysis"));
    theCleanedPtAKJetAnalyzer->setSource("PtAnalysisCleanedAntiKtJets");
    if(theSConeJetAnalyzerFlag){
      theCleanedPtSCJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("CleanedPtAnalysis"));
      theCleanedPtSCJetAnalyzer->setSource("PtAnalysisCleanedSISConeJets");
    }
    if(theIConeJetAnalyzerFlag){
      theCleanedPtICJetAnalyzer  = new JetPtAnalyzer(parameters.getParameter<ParameterSet>("CleanedPtAnalysis"));
      theCleanedPtICJetAnalyzer->setSource("PtAnalysisCleanedIterativeConeJets");
    }
  }

  // --- do the analysis on JPT Jets
  if(theJPTJetAnalyzerFlag) {
    theJPTJetAnalyzer  = new JPTJetAnalyzer(parameters.getParameter<ParameterSet>("JPTJetAnalysis"));
  }
  // --- do the analysis on JPT Cleaned Jets
  if(theJPTJetCleaningFlag) {
    theCleanedJPTJetAnalyzer  = new JPTJetAnalyzer(parameters.getParameter<ParameterSet>("CleanedJPTJetAnalysis"));
  }

  // --- do the analysis on the PFJets
  if(thePFJetAnalyzerFlag)
    thePFJetAnalyzer = new PFJetAnalyzer(parameters.getParameter<ParameterSet>("pfJetAnalysis"));
  // --- do the analysis on the CleanedPFJets
  if(thePFJetCleaningFlag)
    theCleanedPFJetAnalyzer = new PFJetAnalyzer(parameters.getParameter<ParameterSet>("CleanedpfJetAnalysis"));

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
  if(theMuCorrMETAnalyzerFlag){
    theMuCorrMETAnalyzer = new METAnalyzer(parameters.getParameter<ParameterSet>("mucorrMETAnalysis"));
  }
  if(thePfMETAnalyzerFlag){
    thePfMETAnalyzer = new PFMETAnalyzer(parameters.getParameter<ParameterSet>("pfMETAnalysis"));
  }
  if(theHTMHTAnalyzerFlag){
    theHTMHTAnalyzer         = new HTMHTAnalyzer(parameters.getParameter<ParameterSet>("HTMHTAnalysis"));
  }

  _LSBegin     = parameters.getParameter<int>("LSBegin");
  _LSEnd       = parameters.getParameter<int>("LSEnd");

  processname_ = parameters.getParameter<std::string>("processname");

  //jet cleanup parameters
  theCleaningParameters = parameters.getParameter<ParameterSet>("CleaningParameters"),

  _hlt_PhysDec   = theCleaningParameters.getParameter<std::string>("HLT_PhysDec");

  _techTrigsAND  = theCleaningParameters.getParameter<std::vector<unsigned > >("techTrigsAND");
  _techTrigsOR   = theCleaningParameters.getParameter<std::vector<unsigned > >("techTrigsOR");
  _techTrigsNOT  = theCleaningParameters.getParameter<std::vector<unsigned > >("techTrigsNOT");

  _doPVCheck          = theCleaningParameters.getParameter<bool>("doPrimaryVertexCheck");
  _doHLTPhysicsOn     = theCleaningParameters.getParameter<bool>("doHLTPhysicsOn");

  _tightBHFiltering   = theCleaningParameters.getParameter<bool>("tightBHFiltering");
  _tightHcalFiltering = theCleaningParameters.getParameter<bool>("tightHcalFiltering");

  _theVertexLabel     = theCleaningParameters.getParameter<edm::InputTag>("vertexLabel");
  _theGTLabel         = theCleaningParameters.getParameter<edm::InputTag>("gtLabel");

  //Vertex requirements
  if (_doPVCheck) {
    _nvtx_min        = theCleaningParameters.getParameter<int>("nvtx_min");
    _nvtxtrks_min    = theCleaningParameters.getParameter<int>("nvtxtrks_min");
    _vtxndof_min     = theCleaningParameters.getParameter<int>("vtxndof_min");
    _vtxchi2_max     = theCleaningParameters.getParameter<double>("vtxchi2_max");
    _vtxz_max        = theCleaningParameters.getParameter<double>("vtxz_max");
  }


}

// ***********************************************************
JetMETAnalyzer::~JetMETAnalyzer() {   

  if(theJetAnalyzerFlag) {
    delete theAKJetAnalyzer;
    if(theSConeJetAnalyzerFlag) delete theSCJetAnalyzer;    
    if(theIConeJetAnalyzerFlag) delete theICJetAnalyzer;
  }
  if(theJetCleaningFlag) {
    delete theCleanedAKJetAnalyzer;
    if(theSConeJetAnalyzerFlag) delete theCleanedSCJetAnalyzer;
    if(theIConeJetAnalyzerFlag) delete theCleanedICJetAnalyzer;
  }
  if(theDiJetSelectionFlag) {
    delete theDiJetAnalyzer;
  }
  if(theJetPtAnalyzerFlag) {
    delete thePtAKJetAnalyzer;
    if(theSConeJetAnalyzerFlag) delete thePtSCJetAnalyzer;
    if(theIConeJetAnalyzerFlag) delete thePtICJetAnalyzer;
  }

  if(theJetPtCleaningFlag) {
    delete theCleanedPtAKJetAnalyzer;
    if(theSConeJetAnalyzerFlag) delete theCleanedPtSCJetAnalyzer;
    if(theIConeJetAnalyzerFlag) delete theCleanedPtICJetAnalyzer;
  }

  if(theJPTJetAnalyzerFlag)        delete theJPTJetAnalyzer;
  if(theJPTJetCleaningFlag)        delete theCleanedJPTJetAnalyzer;

  if(thePFJetAnalyzerFlag)         delete thePFJetAnalyzer;
  if(thePFJetCleaningFlag)         delete theCleanedPFJetAnalyzer;

  if(theCaloMETAnalyzerFlag){
    delete theCaloMETAnalyzer;
    delete theCaloMETNoHFAnalyzer;
    delete theCaloMETHOAnalyzer;
    delete theCaloMETNoHFHOAnalyzer;
  }
  if(theTcMETAnalyzerFlag)         delete theTcMETAnalyzer;
  if(theMuCorrMETAnalyzerFlag)     delete theMuCorrMETAnalyzer;
  if(thePfMETAnalyzerFlag)         delete thePfMETAnalyzer;
  if(theHTMHTAnalyzerFlag)         delete theHTMHTAnalyzer;

}

// ***********************************************************
void JetMETAnalyzer::beginJob(void) {

  metname = "JetMETAnalyzer";

  LogTrace(metname)<<"[JetMETAnalyzer] Parameters initialization";
  dbe = edm::Service<DQMStore>().operator->();

  //
  //--- Jet
  if(theJetAnalyzerFlag) { 
    theAKJetAnalyzer->beginJob(dbe);
    if(theSConeJetAnalyzerFlag) theSCJetAnalyzer->beginJob(dbe); 
    if(theIConeJetAnalyzerFlag) theICJetAnalyzer->beginJob(dbe);
  }
  if(theJetCleaningFlag) {
    theCleanedAKJetAnalyzer->beginJob(dbe);
    if(theSConeJetAnalyzerFlag) theCleanedSCJetAnalyzer->beginJob(dbe); 
    if(theIConeJetAnalyzerFlag) theCleanedICJetAnalyzer->beginJob(dbe);
  }
  if(theDiJetSelectionFlag) theDiJetAnalyzer->beginJob(dbe); 
  if(theJetPtAnalyzerFlag ) {
    thePtAKJetAnalyzer->beginJob(dbe);
    if(theSConeJetAnalyzerFlag) thePtSCJetAnalyzer->beginJob(dbe);
    if(theIConeJetAnalyzerFlag) thePtICJetAnalyzer->beginJob(dbe);
  }
  if(theJetPtCleaningFlag ) {
    theCleanedPtAKJetAnalyzer->beginJob(dbe);
    if(theSConeJetAnalyzerFlag) theCleanedPtSCJetAnalyzer->beginJob(dbe);
    if(theIConeJetAnalyzerFlag) theCleanedPtICJetAnalyzer->beginJob(dbe);
  }

  if(theJPTJetAnalyzerFlag) theJPTJetAnalyzer->beginJob(dbe);
  if(theJPTJetCleaningFlag) theCleanedJPTJetAnalyzer->beginJob(dbe);

  if(thePFJetAnalyzerFlag)  thePFJetAnalyzer->beginJob(dbe);
  if(thePFJetCleaningFlag)  theCleanedPFJetAnalyzer->beginJob(dbe);

  //
  //--- MET
  if(theCaloMETAnalyzerFlag){
    theCaloMETAnalyzer->beginJob(dbe);
    theCaloMETNoHFAnalyzer->beginJob(dbe);
    theCaloMETHOAnalyzer->beginJob(dbe);
    theCaloMETNoHFHOAnalyzer->beginJob(dbe);
  }
  if(theTcMETAnalyzerFlag) theTcMETAnalyzer->beginJob(dbe);
  if(theMuCorrMETAnalyzerFlag) theMuCorrMETAnalyzer->beginJob(dbe);
  if(thePfMETAnalyzerFlag) thePfMETAnalyzer->beginJob(dbe);
  if(theHTMHTAnalyzerFlag) theHTMHTAnalyzer->beginJob(dbe);
  
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
  bool changed(true);
  _hlt_initialized = hltConfig_.init(iRun,iSetup,processname_,changed);
  if (!_hlt_initialized) {
  //if (!hltConfig_.init(iRun,iSetup,processname_,changed)) {
    processname_ = "FU";
    _hlt_initialized = hltConfig_.init(iRun,iSetup,processname_,changed);
    if(!_hlt_initialized){
      //if (!hltConfig_.init(iRun,iSetup,processname_,changed)){
      LogDebug("JetMETAnalyzer") << "HLTConfigProvider failed to initialize.";
    }
  }

  hltpathME = 0;
  physdecME = 0;
  if (_hlt_initialized) {
  //if (hltConfig_.init(iRun,iSetup,processname_,changed)) {
    if (hltConfig_.size()){
      dbe->setCurrentFolder("JetMET");
      hltpathME = dbe->book1D("hltpath", "hltpath", 300, 0., 300.);
      physdecME = dbe->book1D("physdec", "physdec", 2,   0., 2.);
      if (physdecME) physdecME->setBinLabel(1,"All Events");
    }
    
    for (unsigned int j=0; j!=hltConfig_.size(); ++j) {
      if (hltpathME) hltpathME->setBinLabel(j+1,hltConfig_.triggerName(j));
      if (hltConfig_.triggerName(j)=="HLT_PhysicsDeclared") 
	if (physdecME) physdecME->setBinLabel(2,"PhysicsDeclared");
    }
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
  if(theTcMETAnalyzerFlag) theTcMETAnalyzer->beginRun(iRun, iSetup);
  if(theMuCorrMETAnalyzerFlag) theMuCorrMETAnalyzer->beginRun(iRun, iSetup);
  if(thePfMETAnalyzerFlag) thePfMETAnalyzer->beginRun(iRun, iSetup);
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
  if(theTcMETAnalyzerFlag)     theTcMETAnalyzer->endRun(iRun, iSetup, dbe);
  if(theMuCorrMETAnalyzerFlag) theMuCorrMETAnalyzer->endRun(iRun, iSetup, dbe);
  if(thePfMETAnalyzerFlag)     thePfMETAnalyzer->endRun(iRun, iSetup, dbe);
  //if(theHTMHTAnalyzerFlag)    theHTMHTAnalyzer->endRun(iRun, iSetup, dbe);

}

// ***********************************************************
void JetMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  LogTrace(metname)<<"[JetMETAnalyzer] Analysis of event # ";

  // *** Fill lumisection ME
  int myLuminosityBlock;
  myLuminosityBlock = iEvent.luminosityBlock();
  lumisecME->Fill(myLuminosityBlock);

  if (myLuminosityBlock<_LSBegin) return;
  if (myLuminosityBlock>_LSEnd && _LSEnd>0) return;

  // **** Get the TriggerResults container
  edm::Handle<TriggerResults> triggerResults;
  iEvent.getByLabel(theTriggerResultsLabel, triggerResults);

  // *** Fill trigger results ME
  //if (&triggerResults){
  if (physdecME) physdecME->Fill(0.5);

  bool bPhysicsDeclared = false;
  if(!_doHLTPhysicsOn) bPhysicsDeclared = true;

  if (triggerResults.isValid()){
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
    
    if( triggerNames.triggerIndex("HLT_PhysicsDeclared") != triggerNames.size() )
      if (triggerResults->accept(triggerNames.triggerIndex("HLT_PhysicsDeclared"))) {
	if (physdecME) physdecME->Fill(1.5);
	if(_doHLTPhysicsOn) bPhysicsDeclared = true;
      }

    //sanity check
    assert(triggerResults->size()==hltConfig_.size());
    //check the trigger results
    for (unsigned int j=0; j!=hltConfig_.size(); ++j) {
      if (triggerResults->accept(j)){
        if (hltpathME) hltpathME->Fill(j);
      }
    }
  }

  if (DEBUG)  std::cout << "trigger label " << theTriggerResultsLabel << std::endl;

  Int_t JetLoPass = 0;
  Int_t JetHiPass = 0;

  if (triggerResults.isValid()) {

    if (DEBUG) std::cout << "trigger valid " << std::endl;
    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
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

  // ==========================================================
  //Vertex information
  
  bool bPrimaryVertex = true;
  if(_doPVCheck){
    bPrimaryVertex = false;
    Handle<VertexCollection> vertexHandle;

    iEvent.getByLabel(_theVertexLabel, vertexHandle);

    if (!vertexHandle.isValid()) {
      LogInfo("JetMETAnalyzer") << "JetMETAnalyzer: Could not find vertex collection" << std::endl;
      if (DEBUG) std::cout << "JetMETAnalyzer: Could not find vertex collection" << std::endl;
    }
    
    if ( vertexHandle.isValid() ){
      VertexCollection vertexCollection = *(vertexHandle.product());
      int vertex_number     = vertexCollection.size();
      VertexCollection::const_iterator v = vertexCollection.begin();
      double vertex_chi2    = v->normalizedChi2();
      //double vertex_d0      = sqrt(v->x()*v->x()+v->y()*v->y());
      //double vertex_numTrks = v->tracksSize();
      double vertex_ndof    = v->ndof();
      bool   fakeVtx        = v->isFake();
      double vertex_sumTrks = 0.0;
      double vertex_Z       = v->z();
      for (Vertex::trackRef_iterator vertex_curTrack = v->tracks_begin(); vertex_curTrack!=v->tracks_end(); vertex_curTrack++) {
	vertex_sumTrks += (*vertex_curTrack)->pt();
      }
      
      if (  !fakeVtx
	    && vertex_number>=_nvtx_min
	    //&& vertex_numTrks>_nvtxtrks_min
	    && vertex_ndof   >_vtxndof_min
	    && vertex_chi2   <_vtxchi2_max
	    && fabs(vertex_Z)<_vtxz_max ) bPrimaryVertex = true;
    }
  }
  // ==========================================================

  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( _theGTLabel, gtReadoutRecord);

  if (!gtReadoutRecord.isValid()) {
    LogInfo("JetMETAnalyzer") << "JetMETAnalyzer: Could not find GT readout record" << std::endl;
    if (DEBUG) std::cout << "JetMETAnalyzer: Could not find GT readout record product" << std::endl;
  }
  
  bool bTechTriggers    = true;
  bool bTechTriggersAND = true;
  bool bTechTriggersOR  = false;
  bool bTechTriggersNOT = false;

  if (gtReadoutRecord.isValid()) {
    const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
    
    for (unsigned ttr = 0; ttr != _techTrigsAND.size(); ttr++) {
      bTechTriggersAND = bTechTriggersAND && technicalTriggerWordBeforeMask.at(_techTrigsAND.at(ttr));
    }
    
    for (unsigned ttr = 0; ttr != _techTrigsOR.size(); ttr++) {
      bTechTriggersOR = bTechTriggersOR || technicalTriggerWordBeforeMask.at(_techTrigsOR.at(ttr));
    }
    
    for (unsigned ttr = 0; ttr != _techTrigsNOT.size(); ttr++) {
      bTechTriggersNOT = bTechTriggersNOT || technicalTriggerWordBeforeMask.at(_techTrigsNOT.at(ttr));
    }
  }
    
  bTechTriggers = bTechTriggersAND && bTechTriggersOR && !bTechTriggersNOT;
    
  bool bJetCleanup = bTechTriggers && bPrimaryVertex && bPhysicsDeclared;

  // **** Get the Calo Jet container
  edm::Handle<reco::CaloJetCollection> caloJets;
  

  // **** Get the AntiKt Jet container
  iEvent.getByLabel(theAKJetCollectionLabel, caloJets);    
  if(caloJets.isValid()){
  if(theJetAnalyzerFlag){
    theAKJetAnalyzer->setJetHiPass(JetHiPass);
    theAKJetAnalyzer->setJetLoPass(JetLoPass);
    theAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
  }
  if(theJetPtAnalyzerFlag){
    LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt anti-Kt analyzer";
    thePtAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
  }
  }  

  if(caloJets.isValid() && bJetCleanup){
  if(theJetCleaningFlag){
    theCleanedAKJetAnalyzer->setJetHiPass(JetHiPass);
    theCleanedAKJetAnalyzer->setJetLoPass(JetLoPass);
    theCleanedAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
  }
  if(theJetPtCleaningFlag){
    LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned Jet Pt anti-Kt analyzer";
    theCleanedPtAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
  }
  }


  if(caloJets.isValid() && bJetCleanup){
    if(theDiJetSelectionFlag){
    theDiJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
    }
  }


  if(caloJets.isValid()){
    if(theJetPtAnalyzerFlag){
      LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt anti-Kt analyzer";
      thePtAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
    }
  }

  if(caloJets.isValid() && bJetCleanup){
    if(theJetPtCleaningFlag){
      LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned Jet Pt anti-Kt analyzer";
      theCleanedPtAKJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
    }
  }
  
  
  // **** Get the SISCone Jet container
  iEvent.getByLabel(theSCJetCollectionLabel, caloJets);    
  if(theSConeJetAnalyzerFlag){
    if(caloJets.isValid()){
      if(theJetAnalyzerFlag){
	theSCJetAnalyzer->setJetHiPass(JetHiPass);
	theSCJetAnalyzer->setJetLoPass(JetLoPass);
	theSCJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
      if(theJetPtAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt SisCone analyzer";
	thePtSCJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
    }
    
    if(caloJets.isValid() && bJetCleanup){
      if(theJetCleaningFlag){
	theCleanedSCJetAnalyzer->setJetHiPass(JetHiPass);
	theCleanedSCJetAnalyzer->setJetLoPass(JetLoPass);
	theCleanedSCJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
      if(theJetPtCleaningFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned Jet Pt SisCone analyzer";
	theCleanedPtSCJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
    }
  }
  // **** Get the Iterative Cone Jet container  
  iEvent.getByLabel(theICJetCollectionLabel, caloJets);
  if(theIConeJetAnalyzerFlag) {
    if(caloJets.isValid()){
      if(theJetAnalyzerFlag){
	theICJetAnalyzer->setJetHiPass(JetHiPass);
	theICJetAnalyzer->setJetLoPass(JetLoPass);
	theICJetAnalyzer->analyze(iEvent, iSetup, *caloJets);	
      }
      if(theJetPtAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet Pt ICone analyzer";
	thePtICJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
    }
    
    if(caloJets.isValid() && bJetCleanup){
      if(theJetCleaningFlag){
	theCleanedICJetAnalyzer->setJetHiPass(JetHiPass);
	theCleanedICJetAnalyzer->setJetLoPass(JetLoPass);
	theCleanedICJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
      if(theJetPtCleaningFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned Jet Pt ICone analyzer";
	theCleanedPtICJetAnalyzer->analyze(iEvent, iSetup, *caloJets);
      }
    }
  }

  // **** Get the JPT Jet container
  edm::Handle<reco::JPTJetCollection> jptJets;
  iEvent.getByLabel(theJPTJetCollectionLabel, jptJets);
  if(jptJets.isValid() && theJPTJetAnalyzerFlag){
    //theJPTJetAnalyzer->setJetHiPass(JetHiPass);
    //theJPTJetAnalyzer->setJetLoPass(JetLoPass);
    theJPTJetAnalyzer->analyze(iEvent, iSetup, *jptJets);
  }
  
  if(jptJets.isValid() && bJetCleanup && theJPTJetCleaningFlag){
    //theCleanedJPTJetAnalyzer->setJetHiPass(JetHiPass);
    //theCleanedJPTJetAnalyzer->setJetLoPass(JetLoPass);
    theCleanedJPTJetAnalyzer->analyze(iEvent, iSetup, *jptJets);
  }
  
  // **** Get the PFlow Jet container
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel(thePFJetCollectionLabel, pfJets);

  if(pfJets.isValid()){
    if(thePFJetAnalyzerFlag){
      thePFJetAnalyzer->setJetHiPass(JetHiPass);
      thePFJetAnalyzer->setJetLoPass(JetLoPass);
      LogTrace(metname)<<"[JetMETAnalyzer] Call to the PFJet analyzer";
      thePFJetAnalyzer->analyze(iEvent, iSetup, *pfJets);
    }
    if(thePFJetCleaningFlag){
      theCleanedPFJetAnalyzer->setJetHiPass(JetHiPass);
      theCleanedPFJetAnalyzer->setJetLoPass(JetLoPass);
      LogTrace(metname)<<"[JetMETAnalyzer] Call to the Cleaned PFJet analyzer";
      theCleanedPFJetAnalyzer->analyze(iEvent, iSetup, *pfJets);
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
  // **** MuCorrMETAnalyzer **** //
  //
  if(theMuCorrMETAnalyzerFlag){
     
    theMuCorrMETAnalyzer->analyze(iEvent, iSetup, *triggerResults);
     
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


  //--- MET
  if(theCaloMETAnalyzerFlag){
    theCaloMETAnalyzer->endJob();
    theCaloMETNoHFAnalyzer->endJob();
    theCaloMETHOAnalyzer->endJob();
    theCaloMETNoHFHOAnalyzer->endJob();
  }
  if(theTcMETAnalyzerFlag) theTcMETAnalyzer->endJob();
  if(theMuCorrMETAnalyzerFlag) theMuCorrMETAnalyzer->endJob();
  if(thePfMETAnalyzerFlag)   thePfMETAnalyzer->endJob();
  //if(theHTMHTAnalyzerFlag) theHTMHTAnalyzer->endJob();
  
  //--- Jet
  if(theJetAnalyzerFlag) { 
    theAKJetAnalyzer->endJob();
    if(theSConeJetAnalyzerFlag) theSCJetAnalyzer->endJob(); 
    if(theIConeJetAnalyzerFlag) theICJetAnalyzer->endJob();
  }

  if(theJetCleaningFlag) {
    theCleanedAKJetAnalyzer->endJob();
    if(theSConeJetAnalyzerFlag) theCleanedSCJetAnalyzer->endJob(); 
    if(theIConeJetAnalyzerFlag) theCleanedICJetAnalyzer->endJob();
  }

  if(theDiJetSelectionFlag) theDiJetAnalyzer->endJob();

  if(theJetPtAnalyzerFlag ) {
    thePtAKJetAnalyzer->endJob();
    if(theSConeJetAnalyzerFlag) thePtSCJetAnalyzer->endJob();
    if(theIConeJetAnalyzerFlag) thePtICJetAnalyzer->endJob();
  }

  if(theJetPtCleaningFlag ) {
    theCleanedPtAKJetAnalyzer->endJob();
    if(theSConeJetAnalyzerFlag) theCleanedPtSCJetAnalyzer->endJob();
    if(theIConeJetAnalyzerFlag) theCleanedPtICJetAnalyzer->endJob();
  }

  if(theJPTJetAnalyzerFlag) theJPTJetAnalyzer->endJob();
  if(theJPTJetCleaningFlag) theCleanedJPTJetAnalyzer->endJob();

  if(outputMEsInRootFile){
    dbe->save(outputFileName);
  }
  
}

