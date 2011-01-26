#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/PhotonCoreProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/PhotonProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConvertedPhotonProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionTrackCandidateProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/TrackProducerWithSCAssociation.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionTrackProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionTrackMerger.h"


DEFINE_FWK_MODULE(PhotonCoreProducer);
DEFINE_FWK_MODULE(PhotonProducer);
DEFINE_FWK_MODULE(ConvertedPhotonProducer);
DEFINE_FWK_MODULE(ConversionTrackCandidateProducer);
DEFINE_FWK_MODULE(TrackProducerWithSCAssociation);
DEFINE_FWK_MODULE(ConversionProducer);
DEFINE_FWK_MODULE(ConversionTrackProducer);
DEFINE_FWK_MODULE(ConversionTrackMerger);
