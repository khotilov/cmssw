#include "SUSYBSMAnalysis/HSCP/interface/BetaCalculatorTK.h"

Beta_Calculator_TK::Beta_Calculator_TK(const edm::ParameterSet& iConfig){
  m_trackDeDxEstimatorTag = iConfig.getParameter<edm::InputTag>("dEdXEstimator");
}


void Beta_Calculator_TK::addInfoToCandidate(HSCParticle& candidate, edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   if(!candidate.isTrack())return;

   edm::Handle<DeDxDataValueMap> dedxH;
   iEvent.getByLabel(m_trackDeDxEstimatorTag,dedxH);
   const ValueMap<DeDxData> dEdxTrack = *dedxH.product();

   reco::TrackRef track = candidate.getTrack();
   float k=0.4;
   DeDxBeta result = DeDxBeta(track,dEdxTrack[track],k);
   candidate.setTk(result);
}

