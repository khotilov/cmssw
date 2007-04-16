#include "TrackingTools/KalmanUpdators/interface/KFUpdatorESProducer.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

DEFINE_FWK_EVENTSETUP_MODULE(KFUpdatorESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(Chi2MeasurementEstimatorESProducer);

