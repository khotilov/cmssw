#ifndef EcalEndcapMonitorModule_H
#define EcalEndcapMonitorModule_H

/*
 * \file EcalEndcapMonitorModule.h
 *
 * $Date: 2007/11/27 10:43:26 $
 * $Revision: 1.8 $
 * \author G. Della Ricca
 *
*/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

class MonitorElement;
class DaqMonitorBEInterface;

class EcalEndcapMonitorModule: public edm::EDAnalyzer{

public:

/// Constructor
EcalEndcapMonitorModule(const edm::ParameterSet& ps);

/// Destructor
virtual ~EcalEndcapMonitorModule();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

/// Setup
void setup(void);

/// Cleanup
void cleanup(void);

private:

int runType_;
int evtType_;

int runNumber_;
int evtNumber_;

bool fixedRunNumber_;

int ievt_;

edm::InputTag EcalTBEventHeader_;
edm::InputTag EcalRawDataCollection_;
edm::InputTag EEDigiCollection_;
edm::InputTag EcalUncalibratedRecHitCollection_;
edm::InputTag EcalTrigPrimDigiCollection_;

bool verbose_;

bool enableMonitorDaemon_;

bool enableEventDisplay_;

DaqMonitorBEInterface* dbe_;

bool enableCleanup_;

MonitorElement* meStatus_;

MonitorElement* meRun_;
MonitorElement* meEvt_;

MonitorElement* meRunType_;
MonitorElement* meEvtType_;

MonitorElement* meEEDCC_;

MonitorElement* meEEdigis_[2];
MonitorElement* meEEhits_[2];
MonitorElement* meEEtpdigis_[2];

MonitorElement* meEvent_[18];

bool init_;

};

#endif
