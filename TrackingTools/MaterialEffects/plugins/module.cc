#include "TrackingTools/MaterialEffects/interface/MaterialEffectsUpdator.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterialESProducer.h"


#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

DEFINE_FWK_EVENTSETUP_MODULE(PropagatorWithMaterialESProducer);
