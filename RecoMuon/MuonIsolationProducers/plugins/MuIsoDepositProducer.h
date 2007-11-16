#ifndef MuonIsolationProducers_MuIsoDepositProducer_H
#define MuonIsolationProducers_MuIsoDepositProducer_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/MuonIsolation/interface/MuIsoExtractor.h"
#include <string>

namespace edm { class Event; }
namespace edm { class EventSetup; }

class MuIsoDepositProducer : public edm::EDProducer {

public:

  MuIsoDepositProducer(const edm::ParameterSet&);

  virtual ~MuIsoDepositProducer();

  virtual void produce(edm::Event&, const edm::EventSetup&);
  
private:
  template <typename T>
    void callProduces(const std::vector<std::string>& instLabels);


  edm::ParameterSet theConfig;
  std::string theInputType;

  std::string theOutputType;
  bool theExtractForCandidate;

  std::string theMuonTrackRefType;
  edm::InputTag theMuonCollectionTag;
  std::vector<std::string> theDepositNames;
  bool theMultipleDepositsFlag;
  muonisolation::MuIsoExtractor * theExtractor;

};
#endif
