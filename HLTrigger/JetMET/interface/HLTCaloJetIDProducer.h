#ifndef HLTJetIDProducer_h
#define HLTJetIDProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoJets/JetProducers/interface/JetIDHelper.h"

class HLTJetIDProducer : public edm::EDProducer {
 public:
  explicit HLTJetIDProducer(const edm::ParameterSet&);
  ~HLTJetIDProducer();
  virtual void beginJob() ; 
  virtual void produce(edm::Event &, const edm::EventSetup&);
 private:
  edm::InputTag jetsInput_;
  double min_EMF_;         // minimum EMF
  double max_EMF_;         // maximum EMF
  int min_N90_;            // mininum N90
  int min_N90hits_;        // mininum Nhit90

  reco::helper::JetIDHelper jetID_;

};

#endif
