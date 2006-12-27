#ifndef EMIsolatedTauJetSelector_H
#define EMIsolatedTauJetSelector_H

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Handle.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/BTauReco/interface/EMIsolatedTauTagInfo.h"
#include "DataFormats/BTauReco/interface/JetCrystalsAssociation.h"



class EMIsolatedTauJetsSelector: public edm::EDProducer {
 public:
  explicit EMIsolatedTauJetsSelector(const edm::ParameterSet&);
  ~EMIsolatedTauJetsSelector();
  virtual void produce(edm::Event&, const edm::EventSetup&);

 private:
  edm::InputTag singleTau;
  edm::InputTag doubleTau;
  edm::InputTag l1Code;


};
#endif
