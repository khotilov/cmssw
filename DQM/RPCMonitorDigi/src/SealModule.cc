#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();
#include "DQM/RPCMonitorDigi/interface/RPCMonitorDigi.h"
DEFINE_ANOTHER_FWK_MODULE(RPCMonitorDigi);
#include "DQM/RPCMonitorDigi/interface/RPCTTUMonitor.h"
DEFINE_ANOTHER_FWK_MODULE(RPCTTUMonitor);
