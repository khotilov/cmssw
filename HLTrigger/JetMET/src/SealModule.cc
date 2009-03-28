#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HLTrigger/JetMET/interface/HLTJetVBFFilter.h"
#include "HLTrigger/JetMET/interface/HLTNVFilter.h"
#include "HLTrigger/JetMET/interface/HLTPhi2METFilter.h"
#include "HLTrigger/JetMET/interface/HLTAcoFilter.h"
#include "HLTrigger/JetMET/interface/HLTDiJetAveFilter.h"
#include "HLTrigger/JetMET/interface/HLTRapGapFilter.h"
#include "HLTrigger/JetMET/interface/HLT2jetGapFilter.h"
#include "HLTrigger/JetMET/interface/HLTHPDFilter.h"
#include "HLTrigger/JetMET/interface/HLTMhtHtFilter.h"
#include "HLTrigger/JetMET/interface/HLTHcalMETNoiseFilter.h"


DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HLTJetVBFFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTNVFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTPhi2METFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTAcoFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTDiJetAveFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTRapGapFilter);
DEFINE_ANOTHER_FWK_MODULE(HLT2jetGapFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTHPDFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTMhtHtFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTHcalMETNoiseFilter);
