#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "EventFilter/DTTFRawToDigi/interface/DTTFFEDReader.h"
#include "EventFilter/DTTFRawToDigi/interface/DTTFFEDSim.h"

DEFINE_FWK_MODULE(DTTFFEDReader);
DEFINE_ANOTHER_FWK_MODULE(DTTFFEDSim);
