#ifndef RecoParticleFlow_Benchmark_PFBenchmarkAlgo_h
#define RecoParticleFlow_Benchmark_PFBenchmarkAlgo_h

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include <vector>

class PFBenchmarkAlgo {
public:

  // optional c'tor -- all methods can be used without an instance
  PFBenchmarkAlgo();
  virtual ~PFBenchmarkAlgo();

  // calculate base quantities for the given pair of candidates
  static double deltaEt(const reco::Candidate *, const reco::Candidate *);
  static double deltaEta(const reco::Candidate *, const reco::Candidate *);
  static double deltaPhi(const reco::Candidate *, const reco::Candidate *);
  static double deltaR(const reco::Candidate *, const reco::Candidate *);

  // simple candidate matching
  static const reco::Candidate *matchByDeltaEt(const reco::Candidate *, const reco::CandidateCollection *);
  static const reco::Candidate *matchByDeltaR(const reco::Candidate *, const reco::CandidateCollection *);

  // find a duplicate of the candidate in the collection and return a pointer to that element
  static const reco::Candidate *recoverCandidate(const reco::Candidate *, const reco::CandidateCollection *);

  // sorting functions - returns a sorted copy of the input
  static reco::CandidateCollection sortByDeltaR(const reco::Candidate *, const reco::CandidateCollection *);
  static reco::CandidateCollection sortByDeltaEt(const reco::Candidate *, const reco::CandidateCollection *);

  // multi-match function - returns a sorted, constrained collection
  static reco::CandidateCollection findAllInCone(const reco::Candidate *, const reco::CandidateCollection *, double ConeSize);
  static reco::CandidateCollection findAllInEtWindow(const reco::Candidate *, const reco::CandidateCollection *, double EtWindow);

  // make CandidateCollection out of vector<T> (like a PFCandidateCollection)
  template <typename CandidateDerived>
  static const reco::CandidateCollection *newCandidateCollection(const std::vector<CandidateDerived> *);
  static void deleteCandidateCollection(const reco::CandidateCollection *);

};

// template implementation (required in header)
template <typename CandidateDerived>
const reco::CandidateCollection *PFBenchmarkAlgo::newCandidateCollection(const std::vector<CandidateDerived> *InputCollection) {

  reco::CandidateCollection *copy_candidates = new reco::CandidateCollection();

  for (unsigned int i = 0; i < InputCollection->size(); i++) {
    CandidateDerived *c = (*InputCollection)[i].clone();
    copy_candidates->push_back((CandidateDerived* const)c);
  }

  return copy_candidates;

}

#endif // RecoParticleFlow_Benchmark_PFBenchmarkAlgo_h
