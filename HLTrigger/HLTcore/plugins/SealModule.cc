#include "FWCore/Framework/interface/MakerMacros.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaler.h"
#include "HLTrigger/HLTcore/interface/HLTMakeSummaryObjects.h"
#include "HLTrigger/HLTcore/interface/TriggerSummaryProducerAOD.h"
#include "HLTrigger/HLTcore/interface/TriggerSummaryProducerRAW.h"

DEFINE_FWK_MODULE(HLTPrescaler);
DEFINE_FWK_MODULE(HLTMakeSummaryObjects);
DEFINE_FWK_MODULE(TriggerSummaryProducerAOD);
DEFINE_FWK_MODULE(TriggerSummaryProducerRAW);
