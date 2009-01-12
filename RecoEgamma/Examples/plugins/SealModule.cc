//#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

#include "ElectronIDAnalyzer.h"
#include "ElectronSeedAnalyzer.h"
#include "MCElectronAnalyzer.h"
#include "MCPhotonAnalyzer.h"
#include "MCPizeroAnalyzer.h"
#include "SimpleConvertedPhotonAnalyzer.h"
#include "SimplePhotonAnalyzer.h"
#include "SiStripElectronAnalyzer.h"
#include "PhotonsWithConversionsAnalyzer.h"
#include "GsfElectronMCAnalyzer.h"
#include "GsfElectronDataAnalyzer.h"
#include "GsfElectronFakeAnalyzer.h"
#include "PatPhotonSimpleAnalyzer.h"
DEFINE_SEAL_MODULE();

#include "PhysicsTools/UtilAlgos/interface/Merger.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "PhysicsTools/UtilAlgos/interface/EventSetupInitTrait.h"

typedef Merger<reco::SuperClusterCollection> EgammaSuperClusterMerger;
DEFINE_FWK_MODULE( EgammaSuperClusterMerger );
DEFINE_FWK_MODULE(ElectronIDAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(ElectronSeedAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCElectronAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCPhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(MCPizeroAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(GsfElectronMCAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(GsfElectronDataAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(GsfElectronFakeAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SimpleConvertedPhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SimplePhotonAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(SiStripElectronAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(PhotonsWithConversionsAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(PatPhotonSimpleAnalyzer);
