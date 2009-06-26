#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "Validation/RecoMET/interface/CaloTowerAnalyzer.h"
//#include "Validation/RecoMET/interface/ECALRecHitAnalyzer.h"
//#include "Validation/RecoMET/interface/HCALRecHitAnalyzer.h"
#include "Validation/RecoMET/interface/METTester.h"
#include "Validation/RecoMET/interface/METFileSaver.h"

DEFINE_SEAL_MODULE();

DEFINE_ANOTHER_FWK_MODULE (METFileSaver) ;
DEFINE_ANOTHER_FWK_MODULE (METTester) ;
//DEFINE_ANOTHER_FWK_MODULE (CaloTowerAnalyzer) ;
//DEFINE_ANOTHER_FWK_MODULE (ECALRecHitAnalyzer) ;
//DEFINE_ANOTHER_FWK_MODULE (HCALRecHitAnalyzer) ;

