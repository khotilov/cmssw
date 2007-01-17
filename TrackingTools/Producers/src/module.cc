#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagatorESProducer.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePropagatorESProducer.h"
#include "TrackingTools/GeomPropagators/interface/SmartPropagatorESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"


#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

DEFINE_FWK_EVENTSETUP_MODULE(StraightLinePropagatorESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(AnalyticalPropagatorESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(SmartPropagatorESProducer);
