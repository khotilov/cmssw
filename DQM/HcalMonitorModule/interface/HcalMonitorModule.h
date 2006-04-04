#ifndef HcalMonitorModule_H
#define HcalMonitorModule_H

/*
 * \file HcalMonitorModule.h
 *
 * $Date: 2006/01/03 19:49:45 $
 * $Revision: 1.5 $
 * \author W. Fisher
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQM/HcalMonitorModule/interface/HcalMonitorSelector.h"
#include "DQM/HcalMonitorTasks/interface/HcalDigiMonitor.h"
#include "DQM/HcalMonitorTasks/interface/HcalDataFormatMonitor.h"
#include "DQM/HcalMonitorTasks/interface/HcalRecHitMonitor.h"
#include "DQM/HcalMonitorTasks/interface/HcalPedestalMonitor.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace cms;
using namespace std;

class HcalMonitorModule: public edm::EDAnalyzer{

public:

/// Constructor
HcalMonitorModule(const edm::ParameterSet& ps);

/// Destructor
~HcalMonitorModule();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

private:

  int m_ievt;
  int m_runNum;
  DaqMonitorBEInterface* m_dbe;
  
  MonitorElement* m_meStatus;
  MonitorElement* m_meRunNum;
  MonitorElement* m_meRunType;
  MonitorElement* m_meEvtNum;
  MonitorElement* m_meEvtMask;
  
  HcalMonitorSelector*    m_evtSel;
  HcalDigiMonitor*        m_digiMon;
  HcalDataFormatMonitor*  m_dfMon;
  HcalRecHitMonitor*      m_rhMon;
  HcalPedestalMonitor*    m_pedMon;
  
  bool m_monitorDaemon;

  string m_outputFile;
  ofstream m_logFile;
  
};

#endif
