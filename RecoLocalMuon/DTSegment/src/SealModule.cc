#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();

// ----------- Modules for 2D-segments reco -----------
#include "RecoLocalMuon/DTSegment/src/DTRecSegment2DProducer.h"
#include "RecoLocalMuon/DTSegment/src/DTRecSegment2DAlgoFactory.h"
#include "RecoLocalMuon/DTSegment/src/DTCombinatorialPatternReco.h"

DEFINE_ANOTHER_FWK_MODULE(DTRecSegment2DProducer);
DEFINE_SEAL_PLUGIN (DTRecSegment2DAlgoFactory, DTCombinatorialPatternReco, "DTCombinatorialPatternReco");
//-------------------------------------------------------------------------------------------------------

// ----------- Modules for 4D-segments reco -----------
#include "RecoLocalMuon/DTSegment/src/DTRecSegment4DProducer.h"
#include "RecoLocalMuon/DTSegment/src/DTRecSegment4DAlgoFactory.h"
#include "RecoLocalMuon/DTSegment/src/DTCombinatorialPatternReco4D.h"
#include "RecoLocalMuon/DTSegment/src/DTRefitAndCombineReco4D.h"

DEFINE_ANOTHER_FWK_MODULE(DTRecSegment4DProducer);
DEFINE_SEAL_PLUGIN (DTRecSegment4DAlgoFactory, DTCombinatorialPatternReco4D, "DTCombinatorialPatternReco4D");
DEFINE_SEAL_PLUGIN (DTRecSegment4DAlgoFactory, DTRefitAndCombineReco4D, "DTRefitAndCombineReco4D");
//-------------------------------------------------------------------------------------------------------
