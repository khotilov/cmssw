#include <iostream>

#include "HeavyFlavorAnalysis/InclB/interface/HFTree.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAna00Event.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaTrack.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TGenCand.hh"
#include "AnalysisDataFormats/HeavyFlavorObjects/rootio/TAnaVertex.hh"

// -- Yikes!
TAna00Event  *gHFEvent;

using namespace::std;

// ----------------------------------------------------------------------
HFTree::HFTree(const edm::ParameterSet& iConfig) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTree constructor" << endl;
  cout << "----------------------------------------------------------------------" << endl;
  fFile = new TFile(iConfig.getParameter<string>("fileName").c_str(), "RECREATE");
  fTree = new TTree("T1","CMSSW HF tree");
  fEvent = new TAna00Event(0);
  fTree->Branch("TAna00Event", "TAna00Event", &fEvent, 256000/8, 1);

  gHFEvent = fEvent;
  nevt=0;

}


// ----------------------------------------------------------------------
HFTree::~HFTree() {
  
  // -- Save output
  fFile->cd();
  fTree->Write();
  fFile->Write();
  fFile->Close();
  delete fFile;
}


// ----------------------------------------------------------------------
void HFTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  nevt++;
  gHFEvent->fRunNumber   = iEvent.id().run();
  gHFEvent->fEventNumber = iEvent.id().event();
  cout << "HFTree> filling tree for run: " << gHFEvent->fRunNumber
       << " event: " << gHFEvent->fEventNumber 
       << endl;
  // cout << "HFTree> number of events: " << nevt << endl;
      
  fTree->Fill();
  gHFEvent->Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void  HFTree::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFTree::endJob() { 
  cout << "HFTree> Summary: Events processed: " << nevt << endl;
}

//define this as a plug-in
//DEFINE_FWK_MODULE(HFTree);
