#ifndef L1TRPCTF_H
#define L1TRPCTF_H

/*
 * \file L1TRPCTF.h
 *
 * $Date: 2007/08/27 16:39:18 $
 * $Revision: 1.4 $
 * \author J. Berryhill
 *
*/

// system include files
#include <memory>
#include <unistd.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

#include <iostream>
#include <fstream>
#include <vector>

//
// class decleration
//

class L1TRPCTF : public edm::EDAnalyzer {

public:

// Constructor
L1TRPCTF(const edm::ParameterSet& ps);

// Destructor
virtual ~L1TRPCTF();

protected:
// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

private:
  // ----------member data ---------------------------
  DaqMonitorBEInterface * dbe;

  MonitorElement* rpctfetavalue[3];
  MonitorElement* rpctfphivalue[3];
  MonitorElement* rpctfptvalue[3];
  MonitorElement* rpctfchargevalue[3];
  MonitorElement* rpctfquality[3];
  MonitorElement* rpctfntrack;
  MonitorElement* rpctfbx;

  int nev_; // Number of events processed
  std::string outputFile_; //file name for ROOT ouput
  bool verbose_;
  bool monitorDaemon_;
  ofstream logFile_;
  edm::InputTag rpctfSource_ ;
};

#endif
