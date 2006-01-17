// -*- C++ -*-
//
// Package:    HcalDbProducer
// Class:      HcalDbProducer
// 
/**\class HcalDbProducer HcalDbProducer.h CalibFormats/HcalDbProducer/interface/HcalDbProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Fedor Ratnikov
//         Created:  Tue Aug  9 19:10:10 CDT 2005
// $Id: HcalDbProducer.cc,v 1.10 2006/01/10 19:29:40 fedor Exp $
//
//


// system include files
#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"


#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"
#include "CondFormats/HcalObjects/interface/HcalGainWidths.h"
#include "CondFormats/HcalObjects/interface/HcalGains.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidths.h"
#include "CondFormats/HcalObjects/interface/HcalPedestals.h"
#include "CondFormats/HcalObjects/interface/HcalQIEData.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h" 
#include "CondFormats/DataRecord/interface/HcalElectronicsMapRcd.h"  	 
#include "CondFormats/DataRecord/interface/HcalGainWidthsRcd.h"
#include "CondFormats/DataRecord/interface/HcalGainsRcd.h"
#include "CondFormats/DataRecord/interface/HcalPedestalWidthsRcd.h"
#include "CondFormats/DataRecord/interface/HcalPedestalsRcd.h" 
#include "CondFormats/DataRecord/interface/HcalQIEDataRcd.h"

#include "HcalDbProducer.h"

HcalDbProducer::HcalDbProducer( const edm::ParameterSet&)
  : mService (new HcalDbService ())
{
  //the following line is needed to tell the framework what
  // data is being produced
  setWhatProduced (this, (dependsOn (&HcalDbProducer::pedestalsCallback) &
			  &HcalDbProducer::pedestalWidthsCallback &
			  &HcalDbProducer::gainsCallback &
			  &HcalDbProducer::gainWidthsCallback &
			  &HcalDbProducer::QIEDataCallback &
			  &HcalDbProducer::channelQualityCallback &
			  &HcalDbProducer::electronicsMapCallback
			  )
		   );
  
  //now do what ever other initialization is needed
}


HcalDbProducer::~HcalDbProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
boost::shared_ptr<HcalDbService> HcalDbProducer::produce( const HcalDbRecord&)
{
  std::cout << "HcalDbProducer::produce..." << std::endl;
  return mService;
}

void HcalDbProducer::pedestalsCallback (const HcalPedestalsRcd& fRecord) {
  std::cout << "HcalDbProducer::poolPedestalsCallback->..." << std::endl;
  edm::ESHandle <HcalPedestals> item;
  fRecord.get (item);
  mService->setData (item.product ());
}

  void HcalDbProducer::pedestalWidthsCallback (const HcalPedestalWidthsRcd& fRecord) {
  std::cout << "HcalDbProducer::poolPedestalWidthsCallback->..." << std::endl;
  edm::ESHandle <HcalPedestalWidths> item;
  fRecord.get (item);
  mService->setData (item.product ());
}


  void HcalDbProducer::gainsCallback (const HcalGainsRcd& fRecord) {
  std::cout << "HcalDbProducer::poolGainsCallback->..." << std::endl;
  edm::ESHandle <HcalGains> item;
  fRecord.get (item);
  mService->setData (item.product ());
}


  void HcalDbProducer::gainWidthsCallback (const HcalGainWidthsRcd& fRecord) {
  std::cout << "HcalDbProducer::poolGainWidthsCallback->..." << std::endl;
  edm::ESHandle <HcalGainWidths> item;
  fRecord.get (item);
  mService->setData (item.product ());
}

void HcalDbProducer::QIEDataCallback (const HcalQIEDataRcd& fRecord) {
  std::cout << "HcalDbProducer::QIEDataCallback->..." << std::endl;
  edm::ESHandle <HcalQIEData> item;
  fRecord.get (item);
  mService->setData (item.product ());
}

void HcalDbProducer::channelQualityCallback (const HcalChannelQualityRcd& fRecord) {
  std::cout << "HcalDbProducer::channelQualityCallback->..." << std::endl;
  edm::ESHandle <HcalChannelQuality> item;
  fRecord.get (item);
  mService->setData (item.product ());
}

void HcalDbProducer::electronicsMapCallback (const HcalElectronicsMapRcd& fRecord) {
  std::cout << "HcalDbProducer::electronicsMapCallback->..." << std::endl;
  edm::ESHandle <HcalElectronicsMap> item;
  fRecord.get (item);
  mService->setData (item.product ());
}



