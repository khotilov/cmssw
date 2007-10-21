#ifndef EESummaryClient_H
#define EESummaryClient_H

/*
 * \file EESummaryClient.h
 *
 * $Date: 2007/10/18 08:17:30 $
 * $Revision: 1.7 $
 * \author G. Della Ricca
 *
*/

#include <vector>
#include <string>
#include <fstream>

#include "TROOT.h"
#include "TProfile2D.h"
#include "TH1F.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "OnlineDB/EcalCondDB/interface/EcalCondDBInterface.h"
#include "OnlineDB/EcalCondDB/interface/MonRunIOV.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/MonitorUserInterface.h"
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"

#include "DQM/EcalEndcapMonitorClient/interface/EEClient.h"

class EESummaryClient : public EEClient {

public:

/// Constructor
EESummaryClient(const edm::ParameterSet& ps);

/// Destructor
virtual ~EESummaryClient();

/// Subscribe/Unsubscribe to Monitoring Elements
void subscribe(void);
void subscribeNew(void);
void unsubscribe(void);

/// softReset
void softReset(void);

/// Analyze
void analyze(void);

/// BeginJob
void beginJob(MonitorUserInterface* mui);

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

/// HtmlOutput
void htmlOutput(int run, string htmlDir, string htmlName);

/// WriteDB
bool writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov);

/// Get Functions
inline int getEvtPerJob() { return ievt_; }
inline int getEvtPerRun() { return jevt_; }

/// Set Clients
inline void setFriends(vector<EEClient*> clients) { clients_ = clients; }

private:

void writeMap( std::ofstream& hf, std::string mapname );

int ievt_;
int jevt_;

bool cloneME_;
bool enableQT_;

bool verbose_;

bool enableMonitorDaemon_;

string prefixME_;

vector<int> superModules_;

vector<EEClient*> clients_;

MonitorUserInterface* mui_;
DaqMonitorBEInterface* dbe_;

MonitorElement* meIntegrity_[2];
MonitorElement* meOccupancy_[2];
MonitorElement* mePedestalOnline_[2];
MonitorElement* meLaserL1_[2];
MonitorElement* meLaserL1PN_[2];
MonitorElement* meLed_[2];
MonitorElement* meLedPN_[2];
MonitorElement* mePedestal_[2];
MonitorElement* mePedestalPN_[2];
MonitorElement* meTestPulse_[2];
MonitorElement* meTestPulsePN_[2];

MonitorElement* meCosmic_[2];
MonitorElement* meTiming_[2];
MonitorElement* meTriggerTowerEt_[2];
MonitorElement* meTriggerTowerEmulError_[2];

MEContentsTH2FWithinRangeROOT* qtg01_[2];
MEContentsTH2FWithinRangeROOT* qtg02_[2];
MEContentsTH2FWithinRangeROOT* qtg03_[2];
MEContentsTH2FWithinRangeROOT* qtg04_[2];
MEContentsTH2FWithinRangeROOT* qtg04PN_[2];
MEContentsTH2FWithinRangeROOT* qtg05_[2];
MEContentsTH2FWithinRangeROOT* qtg05PN_[2];
MEContentsTH2FWithinRangeROOT* qtg06_[2];
MEContentsTH2FWithinRangeROOT* qtg06PN_[2];
MEContentsTH2FWithinRangeROOT* qtg07_[2];
MEContentsTH2FWithinRangeROOT* qtg07PN_[2];

MEContentsTH2FWithinRangeROOT* qtg08_[2];
MEContentsTH2FWithinRangeROOT* qtg09_[2];
MEContentsTH2FWithinRangeROOT* qtg10_[2];
MEContentsTH2FWithinRangeROOT* qtg11_[2];

MonitorElement* meGlobalSummary_[2];

MEContentsTH2FWithinRangeROOT* qtg99_[2];

};

#endif
