#include <vector>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"

class TtDecayChannelSelector {
 public:
  enum Leaf
    {
      Elec=0,
      Muon=1,
      Tau =2
    };

  enum TauType
    {
      Leptonic=0,
      OneProng=1,
      ThreeProng=2
    };
    
  typedef std::vector<int> Decay;

  TtDecayChannelSelector(const edm::ParameterSet&);
  ~TtDecayChannelSelector();
  bool operator()(const reco::CandidateCollection&) const;

 private:

  void parseDecayInput(Decay&, Decay&) const;
  void parseTauDecayInput(Decay&) const;
  unsigned int countChargedParticles(const reco::Candidate& part) const;
  bool checkTauDecay(const reco::Candidate&) const;

 private:
  bool  invert_;  //inversion flag
  int   channel_; //top decay channel
  int   summed_;
  Decay decay_;
  Decay chn1_;
  Decay chn2_;
  Decay tauDecay_;
};
