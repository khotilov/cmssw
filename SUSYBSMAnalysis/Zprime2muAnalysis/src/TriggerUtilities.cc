#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/Event.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TriggerUtilities.h"

Zprime2muTriggerPathsAndFilters::Zprime2muTriggerPathsAndFilters(const edm::Event& event) {
  // Here we explicitly specify the HLT path/final filter name as a
  // function of run number, rather than relying on HLTConfigProvider
  // to give it to us, or using globbing on different path version
  // numbers, etc.
    
  // During running, online HLT menus change often so maintaining this
  // list is burdensome. Dealing with subtle bugs/inconsistencies from
  // assumptions made is worse.
    
  // If none of the data/MC/run number conditions below are
  // satisfied, path and filter will remain empty strings, and the
  // valid bit will be set false.

  valid = true; // assume until proven otherwise
  const unsigned run = event.id().run();

  if (!event.isRealData())                 { path = "HLT_Mu15_v2",        filter = "hltL3Muon15",                               prescaled_path = "HLT_Mu15_v2",  prescaled_filter = "hltL3Muon15"; }
  else if (run >= 160404 && run <= 163261) { path = "HLT_Mu30_v1",        filter = "hltSingleMu30L3Filtered30",                 prescaled_path = "HLT_Mu15_v2",  prescaled_filter = "hltL3Muon15"; }
  else if (run >= 163269 && run <= 163869) { path = "HLT_Mu30_v2",        filter = "hltSingleMu30L3Filtered30",                 prescaled_path = "HLT_Mu15_v3",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 165088 && run <= 167043) { path = "HLT_Mu40_v1",        filter = "hltSingleMu40L3Filtered40",                 prescaled_path = "HLT_Mu15_v4",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 166346 && run <= 166346) { path = "HLT_Mu40_v2",        filter = "hltSingleMu40L3Filtered40",                 prescaled_path = "HLT_Mu15_v5",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 167078 && run <= 167913) { path = "HLT_Mu40_v3",        filter = "hltSingleMu40L3Filtered40",                 prescaled_path = "HLT_Mu15_v6",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 170249 && run <= 173198) { path = "HLT_Mu40_v5",        filter = "hltSingleMu40L2QualL3Filtered40",           prescaled_path = "HLT_Mu15_v8",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 173236 && run <= 178380) { path = "HLT_Mu40_eta2p1_v1", filter = "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40", prescaled_path = "HLT_Mu15_v9",  prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 178420 && run <= 179889) { path = "HLT_Mu40_eta2p1_v4", filter = "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40", prescaled_path = "HLT_Mu15_v12", prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else if (run >= 179959 && run <= 180252) { path = "HLT_Mu40_eta2p1_v5", filter = "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40", prescaled_path = "HLT_Mu15_v13", prescaled_filter = "hltSingleMu15L3Filtered15"; }
  else
    valid = false;
}

trigger::TriggerObjectCollection get_L3_muons(const edm::Event& event, const std::string& filter, const edm::InputTag& trigger_summary_src, const std::string& collection) {
  trigger::TriggerObjectCollection L3_muons;

  edm::Handle<trigger::TriggerEvent> trigEvent;
  event.getByLabel(trigger_summary_src, trigEvent);
  if (!trigEvent.isValid())
    throw cms::Exception("get_L3_muons") << "couldn't get hltTriggerSummaryAOD " << trigger_summary_src.encode() << " from event\n";

  // The TriggerEvent object keeps a list of all trigger-firing
  // objects, associated to the original collection module name
  // (usually "hltL3MuonCandidates"). Determine which entries are
  // muons, with "keys" being indices into the TriggerObjectCollection
  // below. If the collection is not found, the final loop will not
  // run thanks to the sentinel values.
  std::pair<int, int> collection_keys(-1,-2);
  int key_prev = 0;
  for (trigger::size_type iC = 0; iC < trigEvent->sizeCollections(); ++iC) {
    int key = trigEvent->collectionKey(iC);
    if (trigEvent->collectionTag(iC).label() == collection) {
      collection_keys = std::make_pair(key_prev, key - 1);
      break;
    }
    key_prev = key;
  }

  // Same idea for the filter name, but a slightly different way to
  // get the keys. If the filter is not found, the keys_passing_filter
  // vector will be empty and no muons will be kept in the loop below.
  trigger::Keys keys_passing_filter;
  for (trigger::size_type iF = 0; iF < trigEvent->sizeFilters(); ++iF) {
    if (trigEvent->filterTag(iF).label() == filter) {
      keys_passing_filter = trigEvent->filterKeys(iF);
      break;
    }
  }

  // Finally, run over the pre-selected range of keys to grab the
  // entries that also passed the filter.
  const trigger::TriggerObjectCollection& TOC(trigEvent->getObjects());
  for (trigger::size_type iO = collection_keys.first; iO <= collection_keys.second; ++iO)
    if (std::find(keys_passing_filter.begin(), keys_passing_filter.end(), iO) != keys_passing_filter.end())
      L3_muons.push_back(TOC.at(iO));

  return L3_muons;
}
