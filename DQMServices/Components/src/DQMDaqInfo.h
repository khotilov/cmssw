#ifndef DQMSERVICES_COMPONENTS_DQMDaqInfo_H
# define DQMSERVICES_COMPONENTS_DQMDaqInfo_H
// -*- C++ -*-
//
// Package:    DQMDaqInfo
// Class:      DQMDaqInfo
// 
/**\class DQMDaqInfo DQMDaqInfo.cc CondCore/DQMDaqInfo/src/DQMDaqInfo.cc
   
 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ilaria SEGONI
//         Created:  Thu Sep 25 11:17:43 CEST 2008
// $Id: DQMDaqInfo.h,v 1.1 2008/10/21 16:00:07 segoni Exp $
//
//

// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

//Run Info
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
#include "CondFormats/RunInfo/interface/RunSummary.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"

//DQM
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"



class DQMDaqInfo : public edm::EDAnalyzer {
public:
  explicit DQMDaqInfo(const edm::ParameterSet&);
  ~DQMDaqInfo();
  

private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& , const  edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& , const  edm::EventSetup&);
  virtual void endJob() ;
  
  bool saveDCFile_;
  std::string outputFile_;
  std::ofstream  dataCertificationFile;
  
  DQMStore *dbe_;  
  std::string outputFileName;
  bool saveData;
  bool FedGranularity;

  enum subDetList { Pixel , SiStrip , ECAL , HCAL , DT , CSC , RPC };  
  
  MonitorElement*  DaqFraction[7];
  
};

#endif
