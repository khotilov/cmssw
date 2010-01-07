#ifndef L1RCTProducer_h
#define L1RCTProducer_h 

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

class L1RCT;
class L1RCTLookupTables;

class L1RCTProducer : public edm::EDProducer
{
 public:
  explicit L1RCTProducer(const edm::ParameterSet& ps);
  virtual ~L1RCTProducer();
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
 private:
  L1RCTLookupTables* rctLookupTables;
  L1RCT* rct;
  bool useEcal;
  bool useHcal;
  std::vector<edm::InputTag> ecalDigis;
  std::vector<edm::InputTag> hcalDigis;
  std::vector<int> bunchCrossings; 
  bool useDebugTpgScales;
};
#endif
