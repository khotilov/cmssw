// -*- C++ -*-
// Original Author:  Fedor Ratnikov
// $Id: HcalTextCalibrations.cc,v 1.11 2009/05/06 22:24:12 mansj Exp $
//
//

#include <memory>
#include <iostream>
#include <fstream>

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "FWCore/Framework/interface/ValidityInterval.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalDbASCIIIO.h"

#include "CondFormats/DataRecord/interface/HcalAllRcds.h"

#include "HcalTextCalibrations.h"
//
// class decleration
//

using namespace cms;

HcalTextCalibrations::HcalTextCalibrations ( const edm::ParameterSet& iConfig ) 
  
{
  //parsing parameters
  std::vector<edm::ParameterSet> data = iConfig.getParameter<std::vector<edm::ParameterSet> >("input");
  std::vector<edm::ParameterSet>::iterator request = data.begin ();
  for (; request != data.end (); request++) {
    std::string objectName = request->getParameter<std::string> ("object");
    edm::FileInPath fp = request->getParameter<edm::FileInPath>("file");
    mInputs [objectName] = fp.fullPath();
    if (objectName == "Pedestals") {
      setWhatProduced (this, &HcalTextCalibrations::producePedestals);
      findingRecord <HcalPedestalsRcd> ();
    }
    else if (objectName == "PedestalWidths") {
      setWhatProduced (this, &HcalTextCalibrations::producePedestalWidths);
      findingRecord <HcalPedestalWidthsRcd> ();
    }
    else if (objectName == "Gains") {
      setWhatProduced (this, &HcalTextCalibrations::produceGains);
      findingRecord <HcalGainsRcd> ();
    }
    else if (objectName == "GainWidths") {
      setWhatProduced (this, &HcalTextCalibrations::produceGainWidths);
      findingRecord <HcalGainWidthsRcd> ();
    }
    else if (objectName == "QIEData") {
      setWhatProduced (this, &HcalTextCalibrations::produceQIEData);
      findingRecord <HcalQIEDataRcd> ();
    }
    else if (objectName == "ChannelQuality") {
      setWhatProduced (this, &HcalTextCalibrations::produceChannelQuality);
      findingRecord <HcalChannelQualityRcd> ();
    }
    else if (objectName == "ZSThresholds") {
      setWhatProduced (this, &HcalTextCalibrations::produceZSThresholds);
      findingRecord <HcalZSThresholdsRcd> ();
    }
    else if (objectName == "RespCorrs") {
      setWhatProduced (this, &HcalTextCalibrations::produceRespCorrs);
      findingRecord <HcalRespCorrsRcd> ();
    }
    else if (objectName == "LUTCorrs") {
      setWhatProduced (this, &HcalTextCalibrations::produceLUTCorrs);
      findingRecord <HcalLUTCorrsRcd> ();
    }
    else if (objectName == "TimeCorrs") {
      setWhatProduced (this, &HcalTextCalibrations::produceTimeCorrs);
      findingRecord <HcalTimeCorrsRcd> ();
    }
    else if (objectName == "L1TriggerObjects") {
      setWhatProduced (this, &HcalTextCalibrations::produceL1TriggerObjects);
      findingRecord <HcalL1TriggerObjectsRcd> ();
    }
    else if (objectName == "ElectronicsMap") {
      setWhatProduced (this, &HcalTextCalibrations::produceElectronicsMap);
      findingRecord <HcalElectronicsMapRcd> ();
    }
    else {
      std::cerr << "HcalTextCalibrations-> Unknown object name '" << objectName 
		<< "', known names are: "
		<< "Pedestals PedestalWidths Gains GainWidths QIEData ChannelQuality ElectronicsMap "
		<< "ZSThresholds RespCorrs LUTCorrs TimeCorrs L1TriggerObjects"
		<< std::endl;
    }
  }
  //  setWhatProduced(this);
}


HcalTextCalibrations::~HcalTextCalibrations()
{
}


//
// member functions
//
void 
HcalTextCalibrations::setIntervalFor( const edm::eventsetup::EventSetupRecordKey& iKey, const edm::IOVSyncValue& iTime, edm::ValidityInterval& oInterval ) {
  std::string record = iKey.name ();
  oInterval = edm::ValidityInterval (edm::IOVSyncValue::beginOfTime(), edm::IOVSyncValue::endOfTime()); //infinite
}

template <class T>
std::auto_ptr<T> produce_impl (const std::string& fFile) {
  std::auto_ptr<T> result (new T ());
  //  std::auto_ptr<T> result;
  std::ifstream inStream (fFile.c_str ());
  if (!inStream.good ()) {
    std::cerr << "HcalTextCalibrations-> Unable to open file '" << fFile << "'" << std::endl;
    throw cms::Exception("FileNotFound") << "Unable to open '" << fFile << "'" << std::endl;
  }
  if (!HcalDbASCIIIO::getObject (inStream, &*result)) {
    std::cerr << "HcalTextCalibrations-> Can not read object from file '" << fFile << "'" << std::endl;
    throw cms::Exception("ReadError") << "Can not read object from file '" << fFile << "'" << std::endl;
  }
  return result;
}



std::auto_ptr<HcalPedestals> HcalTextCalibrations::producePedestals (const HcalPedestalsRcd&) {
  return produce_impl<HcalPedestals> (mInputs ["Pedestals"]);
}

std::auto_ptr<HcalPedestalWidths> HcalTextCalibrations::producePedestalWidths (const HcalPedestalWidthsRcd&) {
  return produce_impl<HcalPedestalWidths> (mInputs ["PedestalWidths"]);
}

std::auto_ptr<HcalGains> HcalTextCalibrations::produceGains (const HcalGainsRcd&) {
  return produce_impl<HcalGains> (mInputs ["Gains"]);
}

std::auto_ptr<HcalGainWidths> HcalTextCalibrations::produceGainWidths (const HcalGainWidthsRcd&) {
  return produce_impl<HcalGainWidths> (mInputs ["GainWidths"]);
}

std::auto_ptr<HcalQIEData> HcalTextCalibrations::produceQIEData (const HcalQIEDataRcd& rcd) {
  return produce_impl<HcalQIEData> (mInputs ["QIEData"]);
}

std::auto_ptr<HcalChannelQuality> HcalTextCalibrations::produceChannelQuality (const HcalChannelQualityRcd& rcd) {
  return produce_impl<HcalChannelQuality> (mInputs ["ChannelQuality"]);
}

std::auto_ptr<HcalZSThresholds> HcalTextCalibrations::produceZSThresholds (const HcalZSThresholdsRcd& rcd) {
  return produce_impl<HcalZSThresholds> (mInputs ["ZSThresholds"]);
}

std::auto_ptr<HcalRespCorrs> HcalTextCalibrations::produceRespCorrs (const HcalRespCorrsRcd& rcd) {
  return produce_impl<HcalRespCorrs> (mInputs ["RespCorrs"]);
}

std::auto_ptr<HcalLUTCorrs> HcalTextCalibrations::produceLUTCorrs (const HcalLUTCorrsRcd& rcd) {
  return produce_impl<HcalLUTCorrs> (mInputs ["LUTCorrs"]);
}

std::auto_ptr<HcalTimeCorrs> HcalTextCalibrations::produceTimeCorrs (const HcalTimeCorrsRcd& rcd) {
  return produce_impl<HcalTimeCorrs> (mInputs ["TimeCorrs"]);
}

std::auto_ptr<HcalL1TriggerObjects> HcalTextCalibrations::produceL1TriggerObjects (const HcalL1TriggerObjectsRcd& rcd) {
  return produce_impl<HcalL1TriggerObjects> (mInputs ["L1TriggerObjects"]);
}

std::auto_ptr<HcalElectronicsMap> HcalTextCalibrations::produceElectronicsMap (const HcalElectronicsMapRcd& rcd) {
  return produce_impl<HcalElectronicsMap> (mInputs ["ElectronicsMap"]);
}

