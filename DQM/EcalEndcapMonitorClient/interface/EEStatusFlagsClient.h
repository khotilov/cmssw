#ifndef EEStatusFlagsClient_H
#define EEStatusFlagsClient_H

/*
 * \file EEStatusFlagsClient.h
 *
 * $Date: 2009/08/27 15:41:03 $
 * $Revision: 1.16 $
 * \author G. Della Ricca
 *
*/

#include <vector>
#include <string>

#include "TROOT.h"
#include "TProfile2D.h"
#include "TH1F.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQM/EcalEndcapMonitorClient/interface/EEClient.h"

class MonitorElement;
class DQMStore;
class EcalCondDBInterface;
class RunIOV;
class MonRunIOV;

class EEStatusFlagsClient : public EEClient {

friend class EESummaryClient;

public:

/// Constructor
EEStatusFlagsClient(const edm::ParameterSet& ps);

/// Destructor
virtual ~EEStatusFlagsClient();

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

/// WriteDB
bool writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov, bool& status);

/// Get Functions
inline int getEvtPerJob() { return ievt_; }
inline int getEvtPerRun() { return jevt_; }

private:

int ievt_;
int jevt_;

bool cloneME_;

bool verbose_;
bool debug_;

std::string prefixME_;

bool enableCleanup_;

std::vector<int> superModules_;

DQMStore* dqmStore_;

MonitorElement* meh01_[18];

TH2F* h01_[18];

MonitorElement* meh02_[18];

TH1F* h02_[18];

MonitorElement* meh03_[18];

TH2F* h03_[18];

};

#endif
