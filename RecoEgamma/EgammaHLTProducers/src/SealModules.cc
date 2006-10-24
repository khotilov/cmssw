#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTEcalIsolationProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTHcalIsolationProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTPhotonTrackIsolationProducers.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTRecoEcalCandidateProducers.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(EgammaHLTEcalIsolationProducers);
DEFINE_ANOTHER_FWK_MODULE(EgammaHLTHcalIsolationProducers);
DEFINE_ANOTHER_FWK_MODULE(EgammaHLTPhotonTrackIsolationProducers);
DEFINE_ANOTHER_FWK_MODULE(EgammaHLTRecoEcalCandidateProducers);
