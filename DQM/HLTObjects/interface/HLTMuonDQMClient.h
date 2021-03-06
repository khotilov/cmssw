#ifndef HLTMuonDQMClient_H
#define HLTMuonDQMClient_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Core/interface/MonitorElementBaseT.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

using namespace edm;
using namespace std;

class HLTMuonDQMClient: public EDAnalyzer{

 public:
  
  HLTMuonDQMClient(const ParameterSet& ps);
  virtual ~HLTMuonDQMClient();
  
 protected:
  
  void beginJob(const EventSetup& c);
  void analyze(const Event& e, const EventSetup& c);
  void endJob();
  void setup();
  void cleanup();
  void doQT();


  void htmlOutput(int run, string htmlDir, string htmlName);

  void CRtd(void);          // initialize the track display
  void PlotTrackDisplay(void);

 private:
  
  bool writeHisto_;
  bool writeHTML_;
  int dumpRate_;
  string outputFile_;
  string outputFileName_;
  string rootFolder_;
  string htmlDir_;
  string htmlName_;
  int count_;
  int run_;
  bool sta_;
  bool init_;

  DaqMonitorBEInterface* dbe_;

  TH1F *hNMu_[4];
  TH1F *hQ_[4];
  TH1F *hPt_[4];
  TH1F *hPtlx_[2];
  TH1F *hEta_[4];
  TH1F *hPhi_[4];
  TH1F *hPtres_[3];
  TH1F *hEtares_[3];
  TH1F *hPhires_[3];
  TProfile *hEtareseta_[3], *hPhiresphi_[3],*hDrphi_[3],*hD0phi_[3],*hDzeta_[3] ;
  TH2F *hEtaphi_[4];
  TH2F *hPtphi_[4];
  TH2F *hPteta_[4];
  TH1F *hD0_[2];
  TH1F *hDr_[2];
  TH1F *hDz_[2];
  TH1F *hErr0_[2];
  TH1F *hNhit_[2];
  TH1F *hIso_[2];
  TH1F *hDimumass_[2];
  TH1F *hL1Qual_;
  bool x11_;                 // enable the ROOT application's x11 interface?
  TApplication *theApp_;     // the ROOT application
  TPad *CRtd_pad_[3];        // Cosmic Run track display pads
};

#endif
