#ifndef EEStatusFlagsTask_H
#define EEStatusFlagsTask_H

/*
 * \file EEStatusFlagsTask.h
 *
 * $Date: 2009/06/18 14:47:10 $
 * $Revision: 1.8 $
 * \author G. Della Ricca
 *
*/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class MonitorElement;
class DQMStore;

class EEStatusFlagsTask: public edm::EDAnalyzer{

public:

/// Constructor
EEStatusFlagsTask(const edm::ParameterSet& ps);

/// Destructor
virtual ~EEStatusFlagsTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

/// BeginJob
void beginJob(void);

/// EndJob
void endJob(void);

/// BeginRun
void beginRun(const edm::Run & r, const edm::EventSetup & c);

/// EndRun
void endRun(const edm::Run & r, const edm::EventSetup & c);

/// Reset
void reset(void);

/// Setup
void setup(void);

/// Cleanup
void cleanup(void);

private:

int ievt_;

DQMStore* dqmStore_;

std::string prefixME_;

bool enableCleanup_;

bool mergeRuns_;

edm::InputTag EcalRawDataCollection_;

MonitorElement* meEvtType_[18];

MonitorElement* meFEchErrors_[18][3];

bool init_;

};

#endif
