#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoParticleFlow/PFClusterProducer/plugins/PFClusterProducer.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFRecHitProducerECAL.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFRecHitProducerHCAL.h"
#include "RecoParticleFlow/PFClusterProducer/plugins/PFRecHitProducerPS.h"



DEFINE_FWK_MODULE(PFClusterProducer);
DEFINE_FWK_MODULE(PFRecHitProducerECAL);
DEFINE_FWK_MODULE(PFRecHitProducerHCAL);
DEFINE_FWK_MODULE(PFRecHitProducerPS);
