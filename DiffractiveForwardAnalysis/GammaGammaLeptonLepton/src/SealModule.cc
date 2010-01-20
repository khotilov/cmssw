#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CosmicsMuMu.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/CollisionsMuMu.h"
#include "DiffractiveForwardAnalysis/GammaGammaLeptonLepton/interface/ExclusiveTrackTrack.h"

DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(CosmicsMuMu);
DEFINE_ANOTHER_FWK_MODULE(CollisionsMuMu);
DEFINE_ANOTHER_FWK_MODULE(ExclusiveTrackTrack);
