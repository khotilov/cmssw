#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoMuon/MuonSeedGenerator/plugins/CosmicMuonSeedGenerator.h"
#include "RecoMuon/MuonSeedGenerator/plugins/MuonSeedGenerator.h"
#include "RecoMuon/MuonSeedGenerator/plugins/MuonSeedProducer.h"
#include "RecoMuon/MuonSeedGenerator/plugins/RPCSeedGenerator.h"
#include "RecoMuon/MuonSeedGenerator/plugins/MuonSeedMerger.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CosmicMuonSeedGenerator);
DEFINE_ANOTHER_FWK_MODULE(MuonSeedGenerator);
DEFINE_ANOTHER_FWK_MODULE(MuonSeedProducer);
DEFINE_ANOTHER_FWK_MODULE(RPCSeedGenerator);
DEFINE_ANOTHER_FWK_MODULE(MuonSeedMerger);

