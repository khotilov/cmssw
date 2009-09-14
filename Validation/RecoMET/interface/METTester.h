#ifndef METTESTER_H
#define METTESTER_H

// author: Mike Schmitt (The University of Florida)
// date: 8/24/2006
// modification: Bobby Scurlock 
// date: 03.11.2006
// note: added RMS(METx) vs SumET capability 
// modification: Rick Cavanaugh
// date: 05.11.2006 
// note: added configuration parameters 
// modification: Mike Schmitt
// date: 02.28.2007
// note: code rewrite

#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include "DQMServices/Core/interface/MonitorElement.h"

class METTester: public edm::EDAnalyzer {
public:

  explicit METTester(const edm::ParameterSet&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
  //virtual void beginJob() ;
  virtual void endJob() ;

 private:

  // DAQ Tools
  DQMStore* dbe_;
  std::map<std::string, MonitorElement*> me;

  // Inputs from Configuration File
  std::string METType_;
  std::string FolderName_;
  edm::InputTag inputMETLabel_;
  edm::InputTag inputTrackLabel_;
  edm::InputTag inputMuonLabel_;
  edm::InputTag inputElectronLabel_;
  edm::InputTag inputBeamSpotLabel_;
  bool finebinning_;

};

#endif // METTESTER_H
