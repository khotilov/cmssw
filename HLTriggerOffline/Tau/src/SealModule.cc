#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "HLTriggerOffline/Tau/interface/HLTTauRefCombiner.h"
#include "HLTriggerOffline/Tau/interface/HLTTauAnalyzer.h"
#include "HLTriggerOffline/Tau/interface/L25TauAnalyzer.h"
#include "HLTriggerOffline/Tau/interface/L2TauAnalyzer.h"
#include "HLTriggerOffline/Tau/interface/L1TauAnalyzer.h"
#include "HLTriggerOffline/Tau/interface/HLTTauValidation.h"
#include "HLTriggerOffline/Tau/interface/TauJetMCFilter.h"
#include "HLTriggerOffline/Tau/interface/ElectronOfETauDQM.h"
#include "HLTriggerOffline/Tau/interface/HLTTauMCProducer.h"
#include "HLTriggerOffline/Tau/interface/HLTMuonTauAnalyzer.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HLTTauRefCombiner);
DEFINE_ANOTHER_FWK_MODULE(HLTTauAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(L2TauAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(L1TauAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(L25TauAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(HLTTauValidation);
DEFINE_ANOTHER_FWK_MODULE(TauJetMCFilter);
DEFINE_ANOTHER_FWK_MODULE(ElectronOfETauDQM);
DEFINE_ANOTHER_FWK_MODULE(HLTTauMCProducer);
DEFINE_ANOTHER_FWK_MODULE(HLTMuonTauAnalyzer);
