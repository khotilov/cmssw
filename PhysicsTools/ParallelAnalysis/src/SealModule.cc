#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/ParallelAnalysis/interface/TSelectorAnalyzer.h"
#include "PhysicsTools/ParallelAnalysis/interface/TrackAnalysisAlgorithm.h"

typedef TSelectorAnalyzer<examples::TrackAnalysisAlgorithm> TrackTSelectorAnalyzer;

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE( TrackTSelectorAnalyzer );
