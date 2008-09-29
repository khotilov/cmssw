#include "ElectroWeakAnalysis/EWKTau/interface/MuonFilter.h"
#include "ElectroWeakAnalysis/EWKTau/interface/PFTauFilter.h"
#include "ElectroWeakAnalysis/EWKTau/interface/CaloTauFilter.h"
#include "ElectroWeakAnalysis/EWKTau/interface/MuTauFilter.h"
#include "ElectroWeakAnalysis/EWKTau/interface/ETauAnalyzer.h"
//#include "ElectroWeakAnalysis/EWKTau/interface/WToTauNuFilter.h"
#include "ElectroWeakAnalysis/EWKTau/interface/DiscriminationByLdgTrackProd.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(MuonFilter);
DEFINE_ANOTHER_FWK_MODULE(PFTauFilter);
DEFINE_ANOTHER_FWK_MODULE(CaloTauFilter);
DEFINE_ANOTHER_FWK_MODULE(MuTauFilter);
DEFINE_ANOTHER_FWK_MODULE(ETauAnalyzer);
DEFINE_ANOTHER_FWK_MODULE(DiscriminationByLdgTrackProd);
//DEFINE_ANOTHER_FWK_MODULE(WToTauNuFilter);

