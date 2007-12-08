//#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

#include "ElectronIDAnalyzer.h"
#include "ElectronAnalyzer.h"
#include "ElectronPixelSeedAnalyzer.h"
#include "MCElectronAnalyzer.h"
#include "MCPhotonAnalyzer.h"
#include "MCPizeroAnalyzer.h"
#include "GsfElectronAnalyzer.h"
#include "SimpleConvertedPhotonAnalyzer.h"
#include "SimplePhotonAnalyzer.h"
#include "SiStripElectronAnalyzer.h"

DEFINE_SEAL_MODULE();

DEFINE_FWK_MODULE(ElectronIDAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(ElectronAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(ElectronPixelSeedAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCElectronAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCPhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCPizeroAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(GsfElectronAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SimpleConvertedPhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SimplePhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SiStripElectronAnalyzer);
