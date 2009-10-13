#ifndef EEDataCertificationTask_h
#define EEDataCertificationTask_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TH2F.h"

class EEDataCertificationTask: public edm::EDAnalyzer{

public:

/// Constructor
EEDataCertificationTask(const edm::ParameterSet& ps);

/// Destructor
virtual ~EEDataCertificationTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

/// BeginJob
void beginJob(const edm::EventSetup& c);

/// EndJob
void endJob(void);

/// BeginLuminosityBlock
void beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const  edm::EventSetup& iSetup);

/// EndLuminosityBlock
void endLuminosityBlock(const edm::LuminosityBlock&  lumiBlock, const  edm::EventSetup& iSetup);

/// Reset
void reset(void);

/// Cleanup
void cleanup(void);
  
private:

bool cloneME_;
  
DQMStore* dqmStore_;

std::string prefixME_;

bool enableCleanup_;

bool mergeRuns_;

TH2F *hDQM_;
TH2F *hDAQ_;
TH2F *hDCS_;

MonitorElement* meEEDataCertificationSummary_;
MonitorElement* meEEDataCertification_[18];
MonitorElement* meEEDataCertificationSummaryMap_;

};

#endif
