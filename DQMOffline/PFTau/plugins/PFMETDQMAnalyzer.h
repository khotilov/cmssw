#ifndef __DQMOffline_PFTau_PFMETDQMAnalyzer__
#define __DQMOffline_PFTau_PFMETDQMAnalyzer__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMOffline/PFTau/interface/PFMETMonitor.h"


class PFMETDQMAnalyzer: public edm::EDAnalyzer {
 public:
  
  PFMETDQMAnalyzer(const edm::ParameterSet& parameterSet);
  
 private:
  void analyze(edm::Event const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob();

  void storeBadEvents(edm::Event const&, float& val);

  edm::InputTag matchLabel_;
  edm::InputTag inputLabel_;
  std::string benchmarkLabel_;
  
  PFMETMonitor pfMETMonitor_;

  edm::ParameterSet pSet_;


};

#endif 
