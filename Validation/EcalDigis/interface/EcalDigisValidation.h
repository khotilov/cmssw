#ifndef EcalDigisValidation_H
#define EcalDigisValidation_H

/*
 * \file EcalDigisValidation.h
 *
 * $Date: 2006/04/10 08:50:50 $
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

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/ESDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

using namespace cms;
using namespace edm;
using namespace std;

class EcalDigisValidation: public EDAnalyzer{

    typedef map<uint32_t,float,less<uint32_t> >  MapType;

public:

/// Constructor
EcalDigisValidation(const ParameterSet& ps);

/// Destructor
~EcalDigisValidation();

protected:

/// Analyze
void analyze(const Event& e, const EventSetup& c);

// BeginJob
void beginJob(const EventSetup& c);

// EndJob
void endJob(void);

private:

 string HepMCLabel;
 string SimTkLabel;
 string SimVtxLabel;
 
 bool verbose_;
 
 DaqMonitorBEInterface* dbe_;
 
 string outputFile_;
 
 map<int, double, less<int> > gainConv_;

 double barrelADCtoGeV_;
 double endcapADCtoGeV_;
 
 MonitorElement* meGunEnergy_;
 MonitorElement* meGunEta_;
 MonitorElement* meGunPhi_;   

 MonitorElement* meEBDigiSimRatio_;
 MonitorElement* meEEDigiSimRatio_;

 MonitorElement* meEBDigiSimRatiogt10ADC_;
 MonitorElement* meEEDigiSimRatiogt10ADC_;

 MonitorElement* meEBDigiSimRatiogt100ADC_;
 MonitorElement* meEEDigiSimRatiogt100ADC_;

};

#endif
