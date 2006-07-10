#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/PhotonProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/SimplePhotonAnalyzer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/PhotonCorrectionProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConvertedPhotonProducer.h"

DEFINE_SEAL_MODULE();


DEFINE_ANOTHER_FWK_MODULE(PhotonProducer);
DEFINE_ANOTHER_FWK_MODULE(SimplePhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(PhotonCorrectionProducer);
DEFINE_ANOTHER_FWK_MODULE(ConvertedPhotonProducer);
