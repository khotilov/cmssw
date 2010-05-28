#ifndef HLTriggerOffline_Egamma_EmDQMPosProcessor_H
#define HLTriggerOffline_Egamma_EmDQMPosProcessor_H


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "TGraphAsymmErrors.h"
#include "DQMServices/Core/interface/DQMStore.h"

class EmDQMPostProcessor : public edm::EDAnalyzer, public TGraphAsymmErrors{
 public:
  EmDQMPostProcessor(const edm::ParameterSet& pset);
  ~EmDQMPostProcessor() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) {};
  void endRun(edm::Run const&, edm::EventSetup const&);
  TProfile* dividehistos(DQMStore * dqm, const std::string& num, const std::string& denom, const std::string& out,const std::string& label, const std::string& titel= "");

 private:
  std::string subDir_;
  
  /** dataset with which these histograms were produced.
      This is set by a user parameter in the configuration
      file.

      It is just used for writing it to the DQM output file.
      Useful to remember with which dataset a histogram file
      was produced. This code does not do much with this information
      (apart from copying it to the output file) but 
      it can be used when generating reports.
  */
  std::string dataSet_;
};

#endif
