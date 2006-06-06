#ifndef CSCDigitizer_CSCDigiDump_h
#define CSCDigitizer_CSCDigiDump_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <string>

class CSCDigiDump : public edm::EDAnalyzer {
public:
  explicit CSCDigiDump(edm::ParameterSet const& conf);
  virtual void analyze(edm::Event const& e, edm::EventSetup const& c);

private:
  edm::ProductTag wireDigiTag_;
  edm::ProductTag stripDigiTag_;
  edm::ProductTag comparatorDigiTag_;
};

#endif

