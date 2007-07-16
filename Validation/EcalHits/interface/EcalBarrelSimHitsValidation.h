#ifndef EcalBarrelSimHitsValidation_H
#define EcalBarrelSimHitsValidation_H

/*
 * \file EcalBarrelSimHitsValidation.h
 *
 * \author C.Rovelli
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

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/EcalValidation/interface/PEcalValidInfo.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


class EcalBarrelSimHitsValidation: public edm::EDAnalyzer{

    typedef std::map<uint32_t,float,std::less<uint32_t> >  MapType;

public:

/// Constructor
EcalBarrelSimHitsValidation(const edm::ParameterSet& ps);

/// Destructor
~EcalBarrelSimHitsValidation();

protected:

/// Analyze
void analyze(const edm::Event& e, const edm::EventSetup& c);

// BeginJob
void beginJob(const edm::EventSetup& c);

// EndJob
void endJob(void);

private:

 uint32_t getUnitWithMaxEnergy(MapType& themap);

 virtual float energyInMatrixEB(int nCellInEta, int nCellInPhi, 
                                int centralEta, int centralPhi, int centralZ,
                                MapType& themap); 
 
 bool  fillEBMatrix(int nCellInEta, int nCellInPhi,
                    int CentralEta, int CentralPhi,int CentralZ,
                    MapType& fillmap, MapType&  themap);
 
 float eCluster2x2( MapType& themap);
 float eCluster4x4(float e33,MapType& themap);

 std::string g4InfoLabel;
 std::string EBHitsCollection;
 std::string ValidationCollection;
 
 bool verbose_;
 
 DaqMonitorBEInterface* dbe_;
 
 std::string outputFile_;

 MonitorElement* menEBHits_;

 MonitorElement* menEBCrystals_;

 MonitorElement* meEBOccupancy_;

 MonitorElement* meEBLongitudinalShower_;

 MonitorElement* meEBhitEnergy_;

 MonitorElement* meEBe1_; 
 MonitorElement* meEBe4_; 
 MonitorElement* meEBe9_; 
 MonitorElement* meEBe16_; 
 MonitorElement* meEBe25_; 

 MonitorElement* meEBe1oe4_;
 MonitorElement* meEBe4oe9_;
 MonitorElement* meEBe9oe16_;
 MonitorElement* meEBe1oe25_;
 MonitorElement* meEBe9oe25_; 
 MonitorElement* meEBe16oe25_;
};

#endif
