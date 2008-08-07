#ifndef EERawDataTask_H
#define EERawDataTask_H

/*
 * \file EERawDataTask.h
 *
 * $Date: 2008/07/12 08:26:48 $
 * $Revision: 1.2 $
 * \author E. Di Marco
 *
*/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class MonitorElement;
class DQMStore;

class EERawDataTask : public edm::EDAnalyzer {

public:

/// Constructor
EERawDataTask(const edm::ParameterSet& ps);

/// Destructor
virtual ~EERawDataTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

/// BeginJob
void beginJob(const edm::EventSetup& c);

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

edm::InputTag FEDRawDataCollection_;
edm::InputTag EcalRawDataCollection_;
edm::InputTag GTEvmSource_;

MonitorElement* meEECRCErrors_;

MonitorElement* meEERunNumberErrors_;
MonitorElement* meEEL1AErrors_;
MonitorElement* meEEOrbitNumberErrors_;
MonitorElement* meEEBunchCrossingErrors_;
MonitorElement* meEETriggerTypeErrors_;

bool init_;

enum activeEVM { TCS, FDLEVM };

};

#endif
