#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTElectronTrackIsolationProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTCombinedIsolationProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTElectronCombinedIsolationProducer.h"

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTElectronDetaDphiProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTRecoEcalCandidateProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTPixelMatchElectronProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTRegionalPixelSeedGeneratorProducers.h"

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTClusterShapeProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTHybridClusterProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTIslandClusterProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTMulti5x5ClusterProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTRemoveDuplicatedSC.h"
// #include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTRemoveSpikesSC.h"

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTR9Producer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTHcalIsolationProducersRegional.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTEcalIsolationProducersRegional.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTEcalRecIsolationProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTHcalIsolationDoubleConeProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTPhotonTrackIsolationProducersRegional.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EcalListOfFEDSProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EcalRecHitsMerger.h"
#include "RecoEgamma/EgammaHLTProducers/interface/ESListOfFEDSProducer.h"
#include "RecoEgamma/EgammaHLTProducers/interface/ESRecHitsMerger.h"

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTNxNClusterProducer.h"


DEFINE_FWK_MODULE(EgammaHLTElectronTrackIsolationProducers);
DEFINE_FWK_MODULE(EgammaHLTElectronDetaDphiProducer);
DEFINE_FWK_MODULE(EgammaHLTRecoEcalCandidateProducers);
DEFINE_FWK_MODULE(EgammaHLTPixelMatchElectronProducers);
DEFINE_FWK_MODULE(EgammaHLTRegionalPixelSeedGeneratorProducers);
DEFINE_FWK_MODULE(EgammaHLTClusterShapeProducer);
DEFINE_FWK_MODULE(EgammaHLTR9Producer);
DEFINE_FWK_MODULE(EgammaHLTHybridClusterProducer);
DEFINE_FWK_MODULE(EgammaHLTIslandClusterProducer);
DEFINE_FWK_MODULE(EgammaHLTMulti5x5ClusterProducer);
DEFINE_FWK_MODULE(EgammaHLTEcalIsolationProducersRegional);
DEFINE_FWK_MODULE(EgammaHLTEcalRecIsolationProducer);
DEFINE_FWK_MODULE(EgammaHLTHcalIsolationProducersRegional);
DEFINE_FWK_MODULE(EgammaHLTHcalIsolationDoubleConeProducers);
DEFINE_FWK_MODULE(EgammaHLTPhotonTrackIsolationProducersRegional);
DEFINE_FWK_MODULE(EgammaHLTRemoveDuplicatedSC);
// DEFINE_FWK_MODULE(EgammaHLTRemoveSpikesSC);
DEFINE_FWK_MODULE(EcalListOfFEDSProducer);
DEFINE_FWK_MODULE(EcalRecHitsMerger);
DEFINE_FWK_MODULE(ESListOfFEDSProducer);
DEFINE_FWK_MODULE(ESRecHitsMerger);
DEFINE_FWK_MODULE(EgammaHLTNxNClusterProducer);
DEFINE_FWK_MODULE(EgammaHLTCombinedIsolationProducer);
DEFINE_FWK_MODULE(EgammaHLTElectronCombinedIsolationProducer);
