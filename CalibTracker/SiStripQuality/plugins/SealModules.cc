#include "PluginManager/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "CalibTracker/SiStripQuality/plugins/SiStripBadModuleByHandBuilder.h"
#include "CalibTracker/SiStripQuality/plugins/SiStripBadStripFromConstructionDB.h"
#include "CalibTracker/SiStripQuality/plugins/SiStripQualityStatistics.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(SiStripBadModuleByHandBuilder);
DEFINE_ANOTHER_FWK_MODULE(SiStripBadStripFromConstructionDB);
DEFINE_ANOTHER_FWK_MODULE(SiStripQualityStatistics);
