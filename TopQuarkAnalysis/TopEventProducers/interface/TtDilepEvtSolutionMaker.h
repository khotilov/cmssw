#include <string>
#include <vector>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include <DataFormats/Candidate/interface/Candidate.h>

class TtDilepEvtSolutionMaker : public edm::EDProducer {

  public:

    explicit TtDilepEvtSolutionMaker(const edm::ParameterSet & iConfig);
    ~TtDilepEvtSolutionMaker();
  
    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

  private:  

    // next methods are avoidable but they make the code legible
    inline bool PTComp(const reco::Candidate*, const reco::Candidate*) const;
    inline bool LepDiffCharge(const reco::Candidate* , const reco::Candidate*) const;
    inline bool HasPositiveCharge(const reco::Candidate*) const;

  private:

    edm::InputTag electronSource_;
    edm::InputTag muonSource_;
    edm::InputTag tauSource_;
    edm::InputTag metSource_;
    edm::InputTag jetSource_;
    edm::InputTag evtSource_;
    unsigned int nrCombJets_;
    bool matchToGenEvt_, calcTopMass_, useMCforBest_;
    bool eeChannel_, emuChannel_, mumuChannel_, etauChannel_, mutauChannel_, tautauChannel_;
    double tmassbegin_, tmassend_, tmassstep_;

};

inline bool TtDilepEvtSolutionMaker::PTComp(const reco::Candidate* l1, const reco::Candidate* l2) const {
  return (l1->pt() > l2->pt());
}

inline bool TtDilepEvtSolutionMaker::LepDiffCharge(const reco::Candidate* l1, const reco::Candidate* l2) const {
  return (l1->charge() != l2->charge());
}

inline bool TtDilepEvtSolutionMaker::HasPositiveCharge(const reco::Candidate* l) const {
  return (l->charge() > 0);
}

