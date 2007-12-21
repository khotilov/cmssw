
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoTracker/FinalTrackSelectors/interface/SimpleTrackListMerger.h"
#include "RecoTracker/FinalTrackSelectors/interface/TrackMultiSelector.h"

using cms::SimpleTrackListMerger;
using reco::modules::TrackMultiSelector;
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(SimpleTrackListMerger);
DEFINE_ANOTHER_FWK_MODULE(TrackMultiSelector);
