#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/PhotonProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConvertedPhotonProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionTrackCandidateProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/TrackProducerWithSCAssociation.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/SoftConversionTrackCandidateProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/SoftConversionProducer.h"
#include "RecoEgamma/EgammaPhotonProducers/interface/TrackerOnlyConversionProducer.h"

DEFINE_SEAL_MODULE();

DEFINE_ANOTHER_FWK_MODULE(PhotonProducer);
DEFINE_ANOTHER_FWK_MODULE(ConvertedPhotonProducer);
DEFINE_ANOTHER_FWK_MODULE(ConversionTrackCandidateProducer);
DEFINE_ANOTHER_FWK_MODULE(TrackProducerWithSCAssociation);
DEFINE_ANOTHER_FWK_MODULE(TrackerOnlyConversionProducer);
DEFINE_ANOTHER_FWK_MODULE(SoftConversionTrackCandidateProducer);
DEFINE_ANOTHER_FWK_MODULE(SoftConversionProducer);
