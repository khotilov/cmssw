#ifndef FOURVECTORHLTCLIENT_H
#define FOURVECTORHLTCLIENT_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile2D.h>

using namespace std;

class FourVectorHLTClient: public edm::EDAnalyzer {

public:

  /// Constructor
  FourVectorHLTClient(const edm::ParameterSet& ps);
  
  /// Destructor
  virtual ~FourVectorHLTClient();
 
protected:

  /// BeginJob
  void beginJob(const edm::EventSetup& c);

  /// BeginRun
  void beginRun(const edm::Run& r, const edm::EventSetup& c);

  /// Fake Analyze
  void analyze(const edm::Event& e, const edm::EventSetup& c) ;

  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                            const edm::EventSetup& context) ;

  /// DQM Client Diagnostic
  void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c);

  /// EndRun
  void endRun(const edm::Run& r, const edm::EventSetup& c);

  /// Endjob
  void endJob();

private:

  void initialize();
  TH1F * get1DHisto(string meName, DQMStore * dbi);
  TH2F * get2DHisto(string meName, DQMStore * dbi);
  TProfile2D * get2DProfile(string meName, DQMStore * dbi);
  TProfile * get1DProfile(string meName, DQMStore * dbi);
  edm::ParameterSet parameters_;

  DQMStore* dbe_;  
  string monitorDir_;
  std::vector<std::string> fourVectorMEName;
  int counterLS_;      ///counter
  int counterEvt_;     ///counter
  int prescaleLS_;     ///units of lumi sections
  int prescaleEvt_;    ///prescale on number of events
  int nChannels;
  Float_t reportSummary;
  Float_t summarySum;
  Float_t summaryContent[20];
  // -------- member data --------

  MonitorElement * reportSummary_;
  MonitorElement * reportSummaryContent_[20];
  MonitorElement * reportSummaryMap_;


};

#endif
