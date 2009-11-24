#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();


//Used to read ME from ROOT files
#include <DQM/RPCMonitorClient/interface/ReadMeFromFile.h>
DEFINE_ANOTHER_FWK_MODULE(ReadMeFromFile);

//General Client
#include <DQM/RPCMonitorClient/interface/RPCDqmClient.h>
DEFINE_ANOTHER_FWK_MODULE(RPCDqmClient);



