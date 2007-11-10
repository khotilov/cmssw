#ifndef EETestPulseTask_H
#define EETestPulseTask_H

/*
 * \file EETestPulseTask.h
 *
 * $Date: 2007/08/14 17:44:47 $
 * $Revision: 1.6 $
 * \author G. Della Ricca
 *
*/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/MonitorElement.h"
class DaqMonitorBEInterface;

class EETestPulseTask: public edm::EDAnalyzer{

public:

/// Constructor
EETestPulseTask(const edm::ParameterSet& ps);

/// Destructor
virtual ~EETestPulseTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

/// BeginJob
void beginJob(const edm::EventSetup& c);

/// EndJob
void endJob(void);

/// Setup
void setup(void);

/// Cleanup
void cleanup(void);

private:

int ievt_;

DaqMonitorBEInterface* dbe_;

bool enableCleanup_;

edm::InputTag EcalRawDataCollection_;
edm::InputTag EEDigiCollection_;
edm::InputTag EcalPnDiodeDigiCollection_;
edm::InputTag EcalUncalibratedRecHitCollection_;

MonitorElement* meShapeMapG01_[18];
MonitorElement* meShapeMapG06_[18];
MonitorElement* meShapeMapG12_[18];

MonitorElement* meAmplMapG01_[18];
MonitorElement* meAmplMapG06_[18];
MonitorElement* meAmplMapG12_[18];

MonitorElement* mePnAmplMapG01_[18];
MonitorElement* mePnAmplMapG16_[18];

MonitorElement* mePnPedMapG01_[18];
MonitorElement* mePnPedMapG16_[18];

// Quality check on crystals, one per each gain

float amplitudeThreshold_;

bool init_;

};

#endif
