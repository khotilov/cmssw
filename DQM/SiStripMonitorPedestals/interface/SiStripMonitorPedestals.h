#ifndef SiStripMonitorPedestals_SiStripMonitorPedestalsTemp_h
#define SiStripMonitorPedestals_SiStripMonitorPedestalsTemp_h
// -*- C++ -*-
//
// Package:     SiStripMonitorPedestals
// Class  :     SiStripMonitorPedestals
// 
/**\class SiStripMonitorPedestals SiStripMonitorPedestals.h 

 Description: <one line class summary>

 Usage:
    <usage>

*/
//
// Original Author:  gennai, dutta
//         Created:  Sat Feb  4 20:49:51 CET 2006
// $Id: SiStripMonitorPedestals.h,v 1.3 2006/07/06 11:08:18 gennai Exp $
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

// data formats
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
//#include "DataFormats/SiStripDigi/interface/SiStripEventSummary.h"
// cabling
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include "CondFormats/SiStripObjects/interface/FedChannelConnection.h"
#include "DataFormats/SiStripDetId/interface/SiStripReadoutKey.h"
//#include "DQM/SiStripCommon/interface/SiStripGenerateKey.h"
//
#include "CalibTracker/SiStripAPVAnalysis/interface/ApvAnalysisFactory.h"

#include "Fed9UUtils.hh"
#include "boost/cstdint.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>


class FEDRawDataCollection;
class FEDRawData;
class SiStripDigi;
class SiStripRawDigi;
class SiStripEventSummary;
class SiStripFedCabling;


class SiStripMonitorPedestals : public edm::EDAnalyzer {
   public:
      explicit SiStripMonitorPedestals(const edm::ParameterSet&);
      ~SiStripMonitorPedestals();

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginJob(edm::EventSetup const&) ;
      virtual void endJob() ;
      
   
 private:
       struct ModMEs{
	 MonitorElement* PedsPerStrip;
	 MonitorElement* PedsDistribution;
	 MonitorElement* PedsEvolution;
	 MonitorElement* CMSubNoisePerStrip;
	 MonitorElement* RawNoisePerStrip;
	 MonitorElement* NoisyStrips;
	 MonitorElement* CMDistribution;
	 
       };
       DaqMonitorBEInterface* dbe_;
       edm::ParameterSet conf_;
       // uint32_t me_type: 1=#digis/module; 2=adcs of hottest strip/module; 3= adcs of coolest strips/module.
       std::map<uint32_t, ModMEs> DigiMEs;
       SiStripFedCabling* fedCabling_;

       bool analyzed;
       
       //The following to be put inside the parametersets
       int16_t nEvUpdate_;
       int16_t signalCutPeds_;
       int16_t nEvTot_;
       int16_t nEvInit_;
       int nIteration_;
       ApvAnalysisFactory* apvFactory_;
       edm::ParameterSet pedsPSet_;
       int  theEventInitNumber_; 
       int theEventIterNumber_;
       int NumCMstripsInGroup_;
       string outPutFileName;
};

#endif
