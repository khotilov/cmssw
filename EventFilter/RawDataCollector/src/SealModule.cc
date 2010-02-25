#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "EventFilter/RawDataCollector/interface/RawDataCollectorModule.h"

using namespace edm::serviceregistry;


DEFINE_FWK_MODULE(RawDataCollectorModule);
