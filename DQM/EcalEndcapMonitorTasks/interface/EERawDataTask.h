#ifndef EERawDataTask_H
#define EERawDataTask_H

/*
 * \file EERawDataTask.h
 *
 * $Date: 2009/04/29 07:46:18 $
 * $Revision: 1.6 $
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

edm::InputTag FEDRawDataCollection_;
edm::InputTag EcalRawDataCollection_;

MonitorElement* meEECRCErrors_;
MonitorElement* meEEEventTypePreCalibrationBX_;
MonitorElement* meEEEventTypeCalibrationBX_;
MonitorElement* meEEEventTypePostCalibrationBX_;
MonitorElement* meEERunNumberErrors_;
MonitorElement* meEEOrbitNumberErrors_;
MonitorElement* meEETriggerTypeErrors_;
MonitorElement* meEECalibrationEventErrors_;
MonitorElement* meEEL1ADCCErrors_;
MonitorElement* meEEBunchCrossingDCCErrors_;
MonitorElement* meEEL1AFEErrors_;
MonitorElement* meEEBunchCrossingFEErrors_;
MonitorElement* meEEL1ATCCErrors_;
MonitorElement* meEEBunchCrossingTCCErrors_;
MonitorElement* meEEL1ASRPErrors_;
MonitorElement* meEEBunchCrossingSRPErrors_;

bool init_;

float calibrationBX_;

enum activeEVM { TCS, FDLEVM };

};

#endif
