#ifndef CALOTOWERANALYZER_H
#define CALOTOWERANALYZER_H

// author: Bobby Scurlock (The University of Florida)
// date: 8/24/2006
// modification: Mike Schmitt
// date: 02.28.2007
// note: code rewrite

#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include "DQMServices/Core/interface/MonitorElement.h"

class CaloTowerAnalyzer: public edm::EDAnalyzer {
public:

  explicit CaloTowerAnalyzer(const edm::ParameterSet&);

  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run& ,const edm::EventSetup&);
  //virtual void beginJob();
  virtual void endJob();

private:

  // DAQ Tools
  DQMStore* dbe_;
  std::map<std::string, MonitorElement*> me;

  // Inputs from Configuration
  edm::InputTag caloTowersLabel_;
  std::vector< edm::InputTag >  HLTBitLabel_ ;
  edm::InputTag HLTResultsLabel_;
  bool debug_;
  double energyThreshold_;
  bool finebinning_;
  bool hltselection_;
  std::string FolderName_;
  int Nevents;
};

#endif
