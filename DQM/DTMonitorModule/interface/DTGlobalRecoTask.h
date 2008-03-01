#ifndef DTGlobalRecoTask_H
#define DTGlobalRecoTask_H

/*
 * \file DTGlobalRecoTask.h
 *
 * $Date: 2006/10/18 18:04:26 $
 * $Revision: 1.2 $
 * \author M. Zanetti - INFN Padova
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include <FWCore/Framework/interface/EDAnalyzer.h>

#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <iostream>
#include <fstream>
#include <vector>

class DTGlobalRecoTask: public edm::EDAnalyzer{

friend class DTMonitorModule;

public:

/// Constructor
DTGlobalRecoTask(const edm::ParameterSet& ps, DQMStore* dbe,
		 const edm::EventSetup& context);

/// Destructor
virtual ~DTGlobalRecoTask();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

private:

  int nevents;
  
  // My monitor elements
  
  std::ofstream logFile;
  
};

#endif
