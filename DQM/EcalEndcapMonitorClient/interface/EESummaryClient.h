#ifndef EESummaryClient_H
#define EESummaryClient_H

/*
 * \file EESummaryClient.h
 *
 * $Date: 2009/03/05 10:41:43 $
 * $Revision: 1.36 $
 * \author G. Della Ricca
 *
*/

#include <vector>
#include <string>
#include <fstream>

#include "TROOT.h"
#include "TProfile2D.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQM/EcalEndcapMonitorClient/interface/EEClient.h"

class MonitorElement;
class DQMStore;
class EcalCondDBInterface;
class RunIOV;
class MonRunIOV;

class EESummaryClient : public EEClient {

public:

/// Constructor
EESummaryClient(const edm::ParameterSet& ps);

/// Destructor
virtual ~EESummaryClient();

/// Analyze
void analyze(void);

/// BeginJob
void beginJob(DQMStore* dqmStore);

/// EndJob
void endJob(void);

/// BeginRun
void beginRun(void);

/// EndRun
void endRun(void);

/// Setup
void setup(void);

/// Cleanup
void cleanup(void);

/// SoftReset
void softReset(bool flag);

/// WriteDB
bool writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov, bool& status, bool flag);

/// Get Functions
inline int getEvtPerJob() { return ievt_; }
inline int getEvtPerRun() { return jevt_; }

/// Set Clients
inline void setFriends(std::vector<EEClient*> clients) { clients_ = clients; }

private:

int ievt_;
int jevt_;

bool cloneME_;

bool verbose_;
bool debug_;

std::string prefixME_;

bool enableCleanup_;

std::vector<int> superModules_;

std::vector<EEClient*> clients_;

DQMStore* dqmStore_;

MonitorElement* meIntegrity_[2];
MonitorElement* meIntegrityErr_;
MonitorElement* meStatusFlags_[2];
MonitorElement* meStatusFlagsErr_;
MonitorElement* meOccupancy_[2];
MonitorElement* meOccupancy1D_;
MonitorElement* mePedestalOnline_[2];
MonitorElement* mePedestalOnlineErr_;
MonitorElement* mePedestalOnlineMean_;
MonitorElement* mePedestalOnlineRMS_;
MonitorElement* mePedestalOnlineRMSMap_[2];
MonitorElement* meLaserL1_[2];
MonitorElement* meLaserL1Err_;
MonitorElement* meLaserL1PN_[2];
MonitorElement* meLaserL1PNErr_;
MonitorElement* meLaserL1Ampl_; 
MonitorElement* meLaserL1Timing_;
MonitorElement* meLaserL1AmplOverPN_;
MonitorElement* meLedL1_[2];
MonitorElement* meLedL1Err_;
MonitorElement* meLedL1PN_[2];
MonitorElement* meLedL1PNErr_;
MonitorElement* meLedL1Ampl_;
MonitorElement* meLedL1Timing_;
MonitorElement* meLedL1AmplOverPN_;
MonitorElement* mePedestal_[2];
MonitorElement* mePedestalG01_[2];
MonitorElement* mePedestalG06_[2];
MonitorElement* mePedestalG12_[2];
MonitorElement* mePedestalErr_;
MonitorElement* mePedestalPN_[2];
MonitorElement* mePedestalPNG01_[2];
MonitorElement* mePedestalPNG16_[2];
MonitorElement* mePedestalPNErr_;
MonitorElement* meTestPulse_[2];
MonitorElement* meTestPulseG01_[2];
MonitorElement* meTestPulseG06_[2];
MonitorElement* meTestPulseG12_[2];
MonitorElement* meTestPulseErr_;
MonitorElement* meTestPulsePN_[2];
MonitorElement* meTestPulsePNG01_[2];
MonitorElement* meTestPulsePNG16_[2];
MonitorElement* meTestPulsePNErr_;
MonitorElement* meTestPulseAmplG01_;
MonitorElement* meTestPulseAmplG06_;
MonitorElement* meTestPulseAmplG12_;

MonitorElement* meCosmic_[2];
MonitorElement* meTiming_[2];
MonitorElement* meTriggerTowerEt_[2];
MonitorElement* meTriggerTowerEtSpectrum_[2];
MonitorElement* meTriggerTowerEmulError_[2];
MonitorElement* meTriggerTowerTiming_[2];

MonitorElement* meGlobalSummary_[2];

TProfile2D* hpot01_[18];
TProfile2D* httt01_[18];

};

#endif
