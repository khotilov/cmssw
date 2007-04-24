// Last commit: $Id: SiStripRawToDigiModule.h,v 1.14 2007/03/21 16:38:13 bainbrid Exp $

#ifndef EventFilter_SiStripRawToDigi_SiStripRawToDigiModule_H
#define EventFilter_SiStripRawToDigi_SiStripRawToDigiModule_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>

class SiStripRawToDigiUnpacker;

/**
   @file EventFilter/SiStripRawToDigi/interface/SiStripRawToDigiModule.h
   @class SiStripRawToDigiModule 
   
   @brief A plug-in module that takes a FEDRawDataCollection as input
   from the Event and creates EDProducts containing StripDigis.
*/
class SiStripRawToDigiModule : public edm::EDProducer {
  
 public:
  
  SiStripRawToDigiModule( const edm::ParameterSet& );
  ~SiStripRawToDigiModule();

  virtual void beginJob( const edm::EventSetup& ) {;}
  virtual void endJob() {;}
  
  virtual void produce( edm::Event&, const edm::EventSetup& );
  
 private: 
  
  SiStripRawToDigiUnpacker* rawToDigi_;

  std::string label_;
  std::string instance_;
  
};

#endif // EventFilter_SiStripRawToDigi_SiStripRawToDigiModule_H

