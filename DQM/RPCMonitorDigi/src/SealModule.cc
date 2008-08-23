#include "DQM/RPCMonitorDigi/interface/RPCMonitorDigi.h"
#include "DQM/RPCMonitorDigi/interface/RPCMonitorSync.h"
#include "DQM/RPCMonitorDigi/interface/RPCMonitorEfficiency.h"
#include "DQM/RPCMonitorDigi/interface/MuonSegmentEff.h"
//#include "DQM/RPCMonitorDigi/interface/RPCReadoutError.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(RPCMonitorDigi);
DEFINE_ANOTHER_FWK_MODULE(RPCMonitorSync);
DEFINE_ANOTHER_FWK_MODULE(RPCMonitorEfficiency);
DEFINE_ANOTHER_FWK_MODULE(MuonSegmentEff);
//DEFINE_ANOTHER_FWK_MODULE(RPCReadoutError);
