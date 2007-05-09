#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

#include "PixelLayerPairsESProducer.h"
#include "PixelLayerTripletsESProducer.h"
#include "MixedLayerTripletsESProducer.h"

DEFINE_FWK_EVENTSETUP_MODULE(PixelLayerPairsESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(PixelLayerTripletsESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(MixedLayerTripletsESProducer);
