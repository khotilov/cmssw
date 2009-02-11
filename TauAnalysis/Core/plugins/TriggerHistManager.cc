#include "TauAnalysis/Core/plugins/TriggerHistManager.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/TriggerNames.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TauAnalysis/Core/interface/histManagerAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

TriggerHistManager::TriggerHistManager(const edm::ParameterSet& cfg)
{
  //std::cout << "<TriggerHistManager::TriggerHistManager>:" << std::endl;

  l1GtReadoutRecordSrc_ = cfg.getParameter<edm::InputTag>("l1GtReadoutRecordSource");
  //std::cout << " l1GtReadoutRecordSrc = " << l1GtReadoutRecordSrc_ << std::endl;

  l1GtObjectMapRecordSrc_ = cfg.getParameter<edm::InputTag>("l1GtObjectMapRecordSource");
  //std::cout << " l1GtObjectMapRecordSrc = " << l1GtObjectMapRecordSrc_ << std::endl;

  l1Bits_ = cfg.getParameter<vstring>("l1Bits");
  //std::cout << " l1Bits = " << format_vstring(l1Bits_) << std::endl;

  hltResultsSrc_ = cfg.getParameter<edm::InputTag>("hltResultsSource");
  //std::cout << " hltResultsSrc = " << hltResultsSrc_ << std::endl;

  hltPaths_ = cfg.getParameter<vstring>("hltPaths");
  //std::cout << " hltPaths = " << format_vstring(hltPaths_) << std::endl;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;
}

TriggerHistManager::~TriggerHistManager()
{
//--- nothing to be done yet...
}

void TriggerHistManager::bookHistograms(const edm::EventSetup& setup)
{
  //std::cout << "<TriggerHistManager::bookHistograms>:" << std::endl;

  if ( edm::Service<DQMStore>().isAvailable() ) {
    DQMStore& dqmStore = (*edm::Service<DQMStore>());

    dqmStore.setCurrentFolder(dqmDirectory_store_);

//--- book histograms for L1 trigger bits
    for ( vstring::const_iterator l1Bit = l1Bits_.begin();
	  l1Bit != l1Bits_.end(); ++l1Bit ) {
      std::string meName = std::string("Trigger").append(*l1Bit);
      MonitorElement* me = dqmStore.book1D(meName, meName, 2, -0.5, 1.5);
      hL1triggerBits_.insert(std::pair<std::string, MonitorElement*>(*l1Bit, me));
    }

//--- book histograms for HLT results
    for ( vstring::const_iterator hltPath = hltPaths_.begin();
	  hltPath != hltPaths_.end(); ++hltPath ) {
      std::string meName = std::string("Trigger").append(*hltPath);
      MonitorElement* me = dqmStore.book1D(meName, meName, 2, -0.5, 1.5);
      hHLTresults_.insert(std::pair<std::string, MonitorElement*>(*hltPath, me));
    }
  }
}

void TriggerHistManager::fillHistograms(const edm::Event& iEvent, const edm::EventSetup& iSetup)

{  
  //std::cout << "<TriggerHistManager::fillHistograms>:" << std::endl; 

//--- fill histograms for L1 trigger bits
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(l1GtReadoutRecordSrc_, l1GtReadoutRecord);
  edm::Handle<L1GlobalTriggerObjectMapRecord> l1GtObjectMapRecord;
  iEvent.getByLabel(l1GtObjectMapRecordSrc_, l1GtObjectMapRecord);

  DecisionWord l1GtDecision = l1GtReadoutRecord->decisionWord();
  const std::vector<L1GlobalTriggerObjectMap>& l1GtObjectMaps = l1GtObjectMapRecord->gtObjectMap();

  for ( vstring::const_iterator l1Bit = l1Bits_.begin();
	l1Bit != l1Bits_.end(); ++l1Bit ) {
    bool isMatch = false;
    for ( std::vector<L1GlobalTriggerObjectMap>::const_iterator l1GtObjectMap = l1GtObjectMaps.begin();
	  l1GtObjectMap != l1GtObjectMaps.end(); ++l1GtObjectMap ) {
      std::string l1Bit_i = (*l1GtObjectMap).algoName();
      if ( l1Bit_i == (*l1Bit) ) {
	int index = (*l1GtObjectMap).algoBitNumber();
	bool isTriggered = l1GtDecision[index];	
	hL1triggerBits_[*l1Bit]->Fill(isTriggered);
	isMatch = true;
      }
    }
    
    if ( !isMatch ) {
      edm::LogError ("TriggerHistManager::fillHistograms") << " Undefined L1 bit = " << (*l1Bit) << " --> skipping !!";
      continue;
    }
  }
   
//--- fill histograms for HLT results 
  edm::Handle<edm::TriggerResults> hltResults;
  iEvent.getByLabel(hltResultsSrc_, hltResults);

  edm::TriggerNames triggerNames;
  triggerNames.init(*hltResults);

  for ( vstring::const_iterator hltPath = hltPaths_.begin();
	hltPath != hltPaths_.end(); ++hltPath ) {
    unsigned int index = triggerNames.triggerIndex(*hltPath);
    if ( index < triggerNames.size() ) {
      bool isTriggered = ( hltResults->accept(index) ) ? true : false;
      hHLTresults_[*hltPath]->Fill(isTriggered);
    } else {
      edm::LogError ("TriggerResultEventSelector::operator") << " Undefined HLT path = " << (*hltPath) << " --> skipping !!";
      continue;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(HistManagerPluginFactory, TriggerHistManager, "TriggerHistManager");

#include "TauAnalysis/Core/interface/HistManagerAdapter.h"

typedef HistManagerAdapter<TriggerHistManager> TriggerAnalyzer;

DEFINE_ANOTHER_FWK_MODULE(TriggerAnalyzer);
