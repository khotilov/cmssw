#include "L1Trigger/L1GctAnalyzer/interface/DumpGctDigis.h"

// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"

using std::string;
using std::ios;
using std::endl;

//
// constructors and destructor
//
DumpGctDigis::DumpGctDigis(const edm::ParameterSet& iConfig) :
  rawLabel_( iConfig.getUntrackedParameter<edm::InputTag>("rawInput", edm::InputTag("L1GctRawDigis") ) ),
  emuRctLabel_( iConfig.getUntrackedParameter<edm::InputTag>("emuRctInput", edm::InputTag("L1RctEmuDigis") ) ),
  emuGctLabel_( iConfig.getUntrackedParameter<edm::InputTag>("emuGctInput", edm::InputTag("L1GctEmuDigis") ) ),
  outFilename_( iConfig.getUntrackedParameter<string>("outFile", "gctAnalyzer.txt") ),
  doHW_( iConfig.getUntrackedParameter<bool>("doHardware", true) ),
  doEmu_( iConfig.getUntrackedParameter<bool>("doEmulated", true) ),
  doRctEM_( iConfig.getUntrackedParameter<bool>("doRctEm", true) ),
  doEM_( iConfig.getUntrackedParameter<bool>("doEm", true) ),
  doRegions_( iConfig.getUntrackedParameter<bool>("doRegions", 0) ),
  doJets_( iConfig.getUntrackedParameter<bool>("doJets", 0) ),
  rctEmMinRank_( iConfig.getUntrackedParameter<unsigned>("rctEmMinRank", 0) ),
  doInternEM_( iConfig.getUntrackedParameter<bool>("doInternEm", true) ),
  doFibres_( iConfig.getUntrackedParameter<bool>("doFibres", false) )
{
  //now do what ever initialization is needed

  outFile_.open(outFilename_.c_str(), ios::out);

}


DumpGctDigis::~DumpGctDigis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  outFile_.close();

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DumpGctDigis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  outFile_ << "Run :" << iEvent.id().run() << "  Event :" << iEvent.id().event() << endl;
  
  // EM
  if (doRctEM_ && doHW_) { doRctEM(iEvent, rawLabel_); }
  if (doRctEM_ && doEmu_) { doRctEM(iEvent, emuRctLabel_); }
  if (doEM_ && doHW_) { doEM(iEvent, rawLabel_); }
  if (doEM_ && doEmu_){ doEM(iEvent, emuGctLabel_); }

  // Jets
  if (doRegions_ && doHW_) { doRegions(iEvent, rawLabel_); }
  if (doRegions_ && doEmu_) { doRegions(iEvent, emuRctLabel_); }
  if (doJets_ && doHW_) { doJets(iEvent, rawLabel_); }
  if (doJets_ && doEmu_) { doJets(iEvent, emuGctLabel_); }

  // debugging
  if (doInternEM_ && doHW_) { doInternEM(iEvent, rawLabel_); }
  if (doFibres_ && doHW_) { doFibres(iEvent, rawLabel_); }

}

void DumpGctDigis::doEM(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1GctEmCandCollection> isoEm;
  Handle<L1GctEmCandCollection> nonIsoEm;

  L1GctEmCandCollection::const_iterator ie;
  L1GctEmCandCollection::const_iterator ne;
  
  iEvent.getByLabel(label.label(),"isoEm",isoEm);
  iEvent.getByLabel(label.label(),"nonIsoEm",nonIsoEm);

  outFile_ << "Iso EM :" << " from : " << label.label() << endl;
  for (ie=isoEm->begin(); ie!=isoEm->end(); ie++) {
    outFile_ << (*ie) << endl;
  } 
  outFile_ << endl;
  
  outFile_ << "Non-iso EM :" << " from : " << label.label() << endl;
  for (ne=nonIsoEm->begin(); ne!=nonIsoEm->end(); ne++) {
    outFile_ << (*ne) << endl;
  } 
  outFile_ << endl;

}

void DumpGctDigis::doRctEM(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1CaloEmCollection> em;

  L1CaloEmCollection::const_iterator e;
 
  iEvent.getByLabel(label, em);

  outFile_ << "RCT EM :" << " from : " << label.label() << endl;
  for (e=em->begin(); e!=em->end(); e++) {
    if (e->rank() >= rctEmMinRank_) {
      outFile_ << (*e) << endl;
    }
  } 
  outFile_ << endl;
  
}


void DumpGctDigis::doRegions(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1CaloRegionCollection> rgns;

  L1CaloRegionCollection::const_iterator r;
  
  iEvent.getByLabel(label, rgns);

  outFile_ << "Regions :" << " from : " << label.label() << endl;
  for (r=rgns->begin(); r!=rgns->end(); r++) {
    outFile_ << (*r) << endl;
  } 
  outFile_ << endl;

}


void DumpGctDigis::doJets(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1GctJetCandCollection> cenJets;
  Handle<L1GctJetCandCollection> forJets;
  Handle<L1GctJetCandCollection> tauJets;
  
  L1GctJetCandCollection::const_iterator cj;
  L1GctJetCandCollection::const_iterator fj;
  L1GctJetCandCollection::const_iterator tj;
  
  
  iEvent.getByLabel(label.label(),"cenJets",cenJets);
  iEvent.getByLabel(label.label(),"forJets",forJets);
  iEvent.getByLabel(label.label(),"tauJets",tauJets);
  
  outFile_ << "Central jets :" << " from : " << label.label() << endl;
  for (cj=cenJets->begin(); cj!=cenJets->end(); cj++) {
    outFile_ << (*cj) << endl;
  } 
  outFile_ << endl;
  
  outFile_ << "Forward jets : " << " from : " << label.label() << endl;
  for (fj=forJets->begin(); fj!=forJets->end(); fj++) {
    outFile_ << (*fj) << endl;
  } 
  outFile_ << endl;
  
  outFile_ << "Tau jets :" << " from : " << label.label() << endl;
  for (tj=tauJets->begin(); tj!=tauJets->end(); tj++) {
    outFile_ << (*tj) << endl;
  } 
  outFile_ << endl; 

}


void DumpGctDigis::doInternEM(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1GctInternEmCandCollection> em;

  L1GctInternEmCandCollection::const_iterator e;
  
  iEvent.getByLabel(label, em);

  outFile_ << "Internal EM :" << " from : " << label.label() << endl;
  for (e=em->begin(); e!=em->end(); e++) {
    outFile_ << (*e) << endl;
  } 
  outFile_ << endl;
  
}


void DumpGctDigis::doFibres(const edm::Event& iEvent, edm::InputTag label) {

  using namespace edm;

  Handle<L1GctFibreCollection> fibres;

  L1GctFibreCollection::const_iterator f;
  
  iEvent.getByLabel(label, fibres);

  outFile_ << "Fibres :" << " from : " << label.label() << endl;
  for (f=fibres->begin(); f!=fibres->end(); f++) {
    outFile_ << (*f) << endl;
  } 
  outFile_ << endl;
  
}

