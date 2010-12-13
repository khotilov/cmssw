#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

namespace {
 namespace {
  susybsm::HSCParticle pa;
/*  susybsm::DriftTubeTOF dtitof;

  susybsm::TimeMeasurement tm;
  std::vector<susybsm::TimeMeasurement> tmv;
  susybsm::MuonTOFCollection mtc; 
  susybsm::MuonTOF mt;
  susybsm::MuonTOFRef mtr;
  susybsm::MuonTOFRefProd mtrp;
  susybsm::MuonTOFRefVector mtrv;
  edm::Wrapper<susybsm::MuonTOFCollection> wr;
  std::vector<std::pair<edm::Ref<std::vector<reco::Muon>,reco::Muon,edm::refhelper::FindUsingAdvance<std::vector<reco::Muon>,reco::Muon> >,susybsm::DriftTubeTOF> > a;
  std::vector<susybsm::DriftTubeTOF> b;
  */
  susybsm::RPCHit4D rpc4h;
  std::vector<susybsm::RPCHit4D> rpc4hv;
  susybsm::CaloBetaMeasurement calobeta;
  susybsm::RPCBetaMeasurement rpcbeta;
//susybsm::DeDxBeta dedxbeta;
  susybsm::HSCParticleCollection hc;
  susybsm::HSCParticleRef hr;
  susybsm::HSCParticleRefProd hp;
  susybsm::HSCParticleRefVector hv;
  edm::Wrapper<susybsm::HSCParticleCollection> wr1;
  
  susybsm::TracksEcalRecHitsMap terhm;
  edm::Wrapper<susybsm::TracksEcalRecHitsMap> wr2;
  edm::helpers::KeyVal<edm::RefProd<std::vector<reco::Track> >,edm::RefProd<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > > hlpr1;

  susybsm::HSCPIsolation hscpI;
  susybsm::HSCPIsolationCollection hscpIc;
  susybsm::HSCPIsolationValueMap hscpIvm; 
  edm::Wrapper<susybsm::HSCPIsolation> hscpIW; 
  edm::Wrapper<susybsm::HSCPIsolationCollection> hscpIcW;
  edm::Wrapper<susybsm::HSCPIsolationValueMap> hscpIvmW;
  
 }
}
