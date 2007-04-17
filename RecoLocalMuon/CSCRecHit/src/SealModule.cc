#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DProducer.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DBuilderPluginFactory.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DFromStripsAndWires.h>
#include <RecoLocalMuon/CSCRecHit/src/CSCRecHit2DFromORedStrips.h>

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CSCRecHit2DProducer);
DEFINE_EDM_PLUGIN(CSCRecHit2DBuilderPluginFactory, CSCRecHit2DFromStripsAndWires, "CSCRecHit2DFromStripsAndWires");
DEFINE_EDM_PLUGIN(CSCRecHit2DBuilderPluginFactory, CSCRecHit2DFromORedStrips, "CSCRecHit2DFromORedStrips");
