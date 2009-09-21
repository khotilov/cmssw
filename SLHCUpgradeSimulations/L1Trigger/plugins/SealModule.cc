#include "SLHCUpgradeSimulations/L1Trigger/interface/L1CaloTowerProducer.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1CaloTriggerSetupProducer.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1CaloClusterProducer.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1CaloClusterFilter.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1JetMaker.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1ExtraMaker.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/SLHCCaloTriggerAccessor.h"
#include "SLHCUpgradeSimulations/L1Trigger/interface/L1DTslhcProd.h"

DEFINE_SEAL_MODULE();
DEFINE_FWK_MODULE(L1CaloTowerProducer);
DEFINE_ANOTHER_FWK_MODULE(L1CaloClusterProducer);
DEFINE_ANOTHER_FWK_MODULE(L1CaloClusterFilter);
DEFINE_ANOTHER_FWK_MODULE(L1JetMaker);
DEFINE_ANOTHER_FWK_MODULE(L1ExtraMaker);
DEFINE_ANOTHER_FWK_MODULE(SLHCCaloTriggerAccessor);
DEFINE_FWK_EVENTSETUP_MODULE(L1CaloTriggerSetupProducer);
DEFINE_FWK_MODULE(L1DTslhcProd);

