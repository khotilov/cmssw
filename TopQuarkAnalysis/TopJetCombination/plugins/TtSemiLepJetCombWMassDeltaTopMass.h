#ifndef TtSemiLepJetCombWMassDeltaTopMass_h
#define TtSemiLepJetCombWMassDeltaTopMass_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

class TtSemiLepJetCombWMassDeltaTopMass : public edm::EDProducer {

 public:
  
  explicit TtSemiLepJetCombWMassDeltaTopMass(const edm::ParameterSet&);
  ~TtSemiLepJetCombWMassDeltaTopMass();
  
 private:

  virtual void beginJob() {};
  virtual void produce(edm::Event& evt, const edm::EventSetup& setup);
  virtual void endJob() {};

  bool isValid(const int& idx, const edm::Handle<std::vector<pat::Jet> >& jets){ return (0<=idx && idx<(int)jets->size()); };

  edm::InputTag jets_;
  edm::InputTag leps_;
  edm::InputTag mets_;
  int maxNJets_;
  double wMass_;
  bool useBTagging_;
  std::string bTagAlgorithm_;
  double minBDiscBJets_;
  double maxBDiscLightJets_;
  int neutrinoSolutionType_;
};

#endif
