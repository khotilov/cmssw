#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <algorithm>

typedef std::vector<reco::PFCandidatePtr> PFCandPtrs; 
typedef PFCandPtrs::iterator PFCandIter;

namespace reco { namespace tau {

class SortPFCandsDescendingPt {
  public:
    bool operator()(const PFCandidatePtr& a, const PFCandidatePtr& b) const {
      return a->pt() > b->pt();
    }
};

std::vector<PFCandidatePtr> flattenPiZeros(const std::vector<RecoTauPiZero>& piZeros) {
  std::vector<PFCandidatePtr> output;

  for(std::vector<RecoTauPiZero>::const_iterator piZero = piZeros.begin();
      piZero != piZeros.end(); ++piZero)
  {
    for(size_t iDaughter = 0; iDaughter < piZero->numberOfDaughters(); ++iDaughter)
    {
      output.push_back(PFCandidatePtr(piZero->daughterPtr(iDaughter)));
    }
  }
  return output;
}

std::vector<reco::PFCandidatePtr> pfCandidates(const reco::PFJet& jet, int particleId, bool sort)
{
  PFCandPtrs pfCands = jet.getPFConstituents();
  PFCandPtrs selectedPFCands;

  for(PFCandIter cand = pfCands.begin(); cand != pfCands.end(); ++cand)
  {
    if( (**cand).particleId() == particleId )
      selectedPFCands.push_back(*cand);
  }

  if ( sort ) std::sort(selectedPFCands.begin(), selectedPFCands.end(), SortPFCandsDescendingPt());

  return selectedPFCands;
}

std::vector<reco::PFCandidatePtr> pfCandidates(const reco::PFJet& jet, const std::vector<int>& particleIds, bool sort)
{
  PFCandPtrs output;
  // Get each desired candidate type, unsorted for now
  for(std::vector<int>::const_iterator particleId = particleIds.begin(); 
      particleId != particleIds.end(); ++particleId) {
    PFCandPtrs selectedPFCands = pfCandidates(jet, *particleId, false);
    output.insert(output.end(), selectedPFCands.begin(), selectedPFCands.end());
  } 
  if (sort) std::sort(output.begin(), output.end(), SortPFCandsDescendingPt());
  return output;
}

std::vector<reco::PFCandidatePtr> pfGammas(const reco::PFJet& jet, bool sort)
{
  return pfCandidates(jet, reco::PFCandidate::gamma, sort);
}

} }
