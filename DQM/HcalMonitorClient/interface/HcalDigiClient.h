#ifndef HcalDigiClient_H
#define HcalDigiClient_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/UI/interface/MonitorUIRoot.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace cms;
using namespace edm;
using namespace std;

class HcalDigiClient{

public:

/// Constructor
HcalDigiClient(const ParameterSet& ps, MonitorUserInterface* mui);
HcalDigiClient();

/// Destructor
virtual ~HcalDigiClient();

/// Subscribe/Unsubscribe to Monitoring Elements
void subscribe(void);
void subscribeNew(void);
void unsubscribe(void);

/// Analyze
void analyze(void);

/// BeginJob
void beginJob(void);

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


 ///process report
  void report();
  
  /// WriteDB
  void htmlOutput(int run, string htmlDir, string htmlName);
  void getHistograms();
  void loadHistograms(TFile* f);

  void errorOutput();
  void getErrors(map<string, vector<QReport*> > out1, map<string, vector<QReport*> > out2, map<string, vector<QReport*> > out3);
  bool hasErrors() const { return dqmReportMapErr_.size(); }
  bool hasWarnings() const { return dqmReportMapWarn_.size(); }
  bool hasOther() const { return dqmReportMapOther_.size(); }

  void resetME();
  void createTests();

private:

  int ievt_;
  int jevt_;
  int kevt_;
  
  bool collateSources_;
  bool cloneME_;
  bool verbose_;
  string process_;

  MonitorUserInterface* mui_;
  
  TH2F* gl_occ_geo[4];
  TH2F* gl_occ_elec[3];
  TH1F* gl_occ_eta;
  TH1F* gl_occ_phi;
  TH2F* gl_err_geo;
  TH2F* gl_err_elec[3];

  TH2F* sub_occ_geo[4][4];
  TH2F* sub_occ_elec[4][3];
  TH1F* sub_occ_eta[4];
  TH1F* sub_occ_phi[4];

  TH2F* sub_err_geo[4];
  TH2F* sub_err_elec[4][3];

  TH2F* geoRef;

  TH1F* qie_adc[4];
  TH1F* num_digi[4];
  TH1F* qie_capid[4];

  // Quality criteria for data integrity
  map<string, vector<QReport*> > dqmReportMapErr_;
  map<string, vector<QReport*> > dqmReportMapWarn_;
  map<string, vector<QReport*> > dqmReportMapOther_;
  map<string, string> dqmQtests_;

};

#endif
