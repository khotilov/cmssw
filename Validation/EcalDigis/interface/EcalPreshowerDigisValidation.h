#ifndef EcalPreshowerDigisValidation_H
#define EcalPreshowerDigisValidation_H

/*
 * \file EcalPreshowerDigisValidation.h
 *
 * $Date: 2006/06/20 16:25:59 $
 * $Revision: 1.2 $
 * \author F. Cossutti
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace cms;
using namespace edm;
using namespace std;

class EcalPreshowerDigisValidation: public EDAnalyzer{

    typedef map<uint32_t,float,less<uint32_t> >  MapType;

public:

/// Constructor
EcalPreshowerDigisValidation(const ParameterSet& ps);

/// Destructor
~EcalPreshowerDigisValidation();

protected:

/// Analyze
void analyze(const Event& e, const EventSetup& c);

// BeginJob
void beginJob(const EventSetup& c);

// EndJob
void endJob(void);

private:

 bool verbose_;
 
 DaqMonitorBEInterface* dbe_;
 
 string outputFile_;

 MonitorElement* meESDigiMultiplicity_;
 
 MonitorElement* meESDigiADC_[3];

};

#endif
