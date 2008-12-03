#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_SEAL_MODULE();

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripFEDRawDataAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripFEDRawDataAnalyzer);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripDigiAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripDigiAnalyzer);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripTrivialClusterSource.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripTrivialClusterSource);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripTrivialDigiSource.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripTrivialDigiSource);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripRawToClustersDummyUnpacker.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripRawToClustersDummyUnpacker);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripClustersDSVBuilder.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripClustersDSVBuilder);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripDigiValidator.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripDigiValidator);

#include "EventFilter/SiStripRawToDigi/test/plugins/SiStripModuleTimer.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripModuleTimer);
