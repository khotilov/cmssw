/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/09/12 19:13:26 $
 *  $Revision: 1.8 $
 *  \author F. Chlebana - Fermilab
 */

#include "DQMOffline/JetMET/src/JetMETAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
using namespace std;
using namespace edm;


JetMETAnalyzer::JetMETAnalyzer(const edm::ParameterSet& pSet) {

  cout<<"[JetMETAnalyzer] Constructor called!"<<endl;

  parameters = pSet;
  
  // Calo Jet Collection Label
  theSCJetCollectionLabel   = parameters.getParameter<edm::InputTag>("SCJetsCollectionLabel");
  theICJetCollectionLabel   = parameters.getParameter<edm::InputTag>("ICJetsCollectionLabel");
//theCaloJetCollectionLabel = parameters.getParameter<edm::InputTag>("CaloJetsCollectionLabel");

  thePFJetCollectionLabel   = parameters.getParameter<edm::InputTag>("PFJetsCollectionLabel");

  theCaloMETCollectionLabel     = parameters.getParameter<edm::InputTag>("CaloMETCollectionLabel");
  theCaloMETNoHFCollectionLabel = parameters.getParameter<edm::InputTag>("CaloMETNoHFCollectionLabel");
  
//theSCJetAnalyzerFlag      = parameters.getUntrackedParameter<bool>("DoSCJetAnalysis",true);
//theICJetAnalyzerFlag      = parameters.getUntrackedParameter<bool>("DoICJetAnalysis",true);
  theJetAnalyzerFlag        = parameters.getUntrackedParameter<bool>("DoJetAnalysis",true);

  thePFJetAnalyzerFlag      = parameters.getUntrackedParameter<bool>("DoPFJetAnalysis",true);
  theCaloMETAnalyzerFlag    = parameters.getUntrackedParameter<bool>("DoCaloMETAnalysis",true);

  // --- do the analysis on the Jets
  if(theJetAnalyzerFlag) {
    //    theJetAnalyzer    = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    //    theJetAnalyzer->setSource("CaloJets");
    theSCJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    theSCJetAnalyzer->setSource("SISConeJets");
    theICJetAnalyzer  = new JetAnalyzer(parameters.getParameter<ParameterSet>("jetAnalysis"));
    theICJetAnalyzer->setSource("IterativeConeJets");

  }

  // --- do the analysis on the PFJets
  if(thePFJetAnalyzerFlag)
    thePFJetAnalyzer = new PFJetAnalyzer(parameters.getParameter<ParameterSet>("pfJetAnalysis"));

  // --- do the analysis on the MET
  if(theCaloMETAnalyzerFlag)
    theCaloMETAnalyzer = new CaloMETAnalyzer(parameters.getParameter<ParameterSet>("caloMETAnalysis"));

}

JetMETAnalyzer::~JetMETAnalyzer() {   
  if(theJetAnalyzerFlag) {
    //    delete theJetAnalyzer;
    delete theSCJetAnalyzer;
    delete theICJetAnalyzer;
  }
  if(thePFJetAnalyzerFlag)   delete thePFJetAnalyzer;
  if(theCaloMETAnalyzerFlag) delete theCaloMETAnalyzer;
}


void JetMETAnalyzer::beginJob(edm::EventSetup const& iSetup) {

  metname = "JetMETAnalyzer";

  LogTrace(metname)<<"[JetMETAnalyzer] Parameters initialization";
  dbe = edm::Service<DQMStore>().operator->();

  if(theJetAnalyzerFlag) { 
    //    theJetAnalyzer->beginJob(iSetup, dbe);
    theSCJetAnalyzer->beginJob(iSetup, dbe);
    theICJetAnalyzer->beginJob(iSetup, dbe);
  }

  if(thePFJetAnalyzerFlag)    thePFJetAnalyzer->beginJob(iSetup, dbe);
  if(theCaloMETAnalyzerFlag)  theCaloMETAnalyzer->beginJob(iSetup, dbe);

}


void JetMETAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  LogTrace(metname)<<"[JetMETAnalyzer] Analysis of event # ";
  
  // **** Get the TriggerResults container
  edm::Handle<TriggerResults> triggerResults;
  iEvent.getByLabel(theTriggerResultsLabel,triggerResults);

  // **** Get the Calo Jet container
  edm::Handle<reco::CaloJetCollection> caloJets;

  /****
  iEvent.getByLabel(theCaloJetCollectionLabel, caloJets);
  if(caloJets.isValid()){
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the Jet analyzer";
	theJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }
  ****/

  // **** Get the SISCone Jet container
  iEvent.getByLabel(theSCJetCollectionLabel, caloJets);

  if(caloJets.isValid()){
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the SC Jet analyzer";
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

  // **** Get the Iterative Cone  Jet container
  iEvent.getByLabel(theICJetCollectionLabel, caloJets);

  if(caloJets.isValid()){
    for (reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal!=caloJets->end(); ++cal){
      if(theJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the IC Jet analyzer";
	if (cal == caloJets->begin()) {
	  theICJetAnalyzer->setNJets(caloJets->size());
	  theICJetAnalyzer->setLeadJetFlag(1);
	} else {
	  theICJetAnalyzer->setLeadJetFlag(0);
	}
	theICJetAnalyzer->analyze(iEvent, iSetup, *cal);	
      }
    }
  }
  
  // **** Get the PFlow Jet container
  edm::Handle<reco::PFJetCollection> pfJets;
  iEvent.getByLabel(thePFJetCollectionLabel, pfJets);

  if(pfJets.isValid()){
    for (reco::PFJetCollection::const_iterator cal = pfJets->begin(); cal!=pfJets->end(); ++cal){
      if(thePFJetAnalyzerFlag){
	LogTrace(metname)<<"[JetMETAnalyzer] Call to the PFJet analyzer";
	thePFJetAnalyzer->analyze(iEvent, iSetup, *cal);
      }
    }
  }


  // **** Get the MET container  
  edm::Handle<reco::CaloMETCollection> calometcoll;
  iEvent.getByLabel(theCaloMETCollectionLabel, calometcoll);
  edm::Handle<reco::CaloMETCollection> calometNoHFcoll;
  iEvent.getByLabel(theCaloMETNoHFCollectionLabel, calometNoHFcoll);

  if(calometcoll.isValid() && calometNoHFcoll.isValid()){
    const CaloMETCollection *calometcol = calometcoll.product();
    const CaloMET *calomet;
    calomet = &(calometcol->front());
    const CaloMETCollection *calometNoHFcol = calometNoHFcoll.product();
    const CaloMET *calometNoHF;
    calometNoHF = &(calometNoHFcol->front());

    if(theCaloMETAnalyzerFlag){
      LogTrace(metname)<<"[JetMETAnalyzer] Call to the CaloMET analyzer";
      theCaloMETAnalyzer->analyze(iEvent, iSetup,
				  *triggerResults,
				  *calomet, *calometNoHF);
    }
  }


}


void JetMETAnalyzer::endJob(void) {
  LogTrace(metname)<<"[JetMETAnalyzer] Saving the histos";
  bool outputMEsInRootFile   = parameters.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = parameters.getParameter<std::string>("OutputFileName");
  if(outputMEsInRootFile){
    dbe->showDirStructure();
    dbe->save(outputFileName);
  }
}

