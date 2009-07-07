#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPEESProducer.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPEfromTrackAngleESProducer.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPEfromTrackAngle2ESProducer.h"

#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitMatcherESProducer.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitConverter.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(StripCPEESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(StripCPEfromTrackAngleESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(StripCPEfromTrackAngle2ESProducer);
DEFINE_ANOTHER_FWK_EVENTSETUP_MODULE(SiStripRecHitMatcherESProducer);
using cms::SiStripRecHitConverter;
DEFINE_ANOTHER_FWK_MODULE(SiStripRecHitConverter);

