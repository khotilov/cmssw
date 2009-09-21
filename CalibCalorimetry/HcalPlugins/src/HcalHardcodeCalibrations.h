//
// Original Author:  Fedor Ratnikov Oct 21, 2005
// $Id: HcalHardcodeCalibrations.h,v 1.12 2009/07/30 20:09:03 kukartse Exp $
//
// ESSource to generate default HCAL calibration objects 
//
#include <map>
#include <string>

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/HcalObjects/interface/AllObjects.h"
class ParameterSet;

class HcalPedestalsRcd;
class HcalPedestalWidthsRcd;
class HcalGainsRcd;
class HcalGainWidthsRcd;
class HcalQIEDataRcd;
class HcalChannelQualityRcd;
class HcalElectronicsMapRcd;
class HcalRespCorrsRcd;
class HcalZSThresholdsRcd;
class HcalL1TriggerObjectsRcd;
class HcalTimeCorrsRcd;
class HcalLUTCorrsRcd;
class HcalPFCorrsRcd;
class HcalValidationCorrsRcd;
class HcalLutMetadataRcd;

class HcalHardcodeCalibrations : public edm::ESProducer,
		       public edm::EventSetupRecordIntervalFinder
{
public:
  HcalHardcodeCalibrations (const edm::ParameterSet& );
  ~HcalHardcodeCalibrations ();

  void produce () {};
  
protected:
  virtual void setIntervalFor(const edm::eventsetup::EventSetupRecordKey&,
			      const edm::IOVSyncValue& , 
			      edm::ValidityInterval&) ;

  std::auto_ptr<HcalPedestals> producePedestals (const HcalPedestalsRcd& rcd);
  std::auto_ptr<HcalPedestalWidths> producePedestalWidths (const HcalPedestalWidthsRcd& rcd);
  std::auto_ptr<HcalGains> produceGains (const HcalGainsRcd& rcd);
  std::auto_ptr<HcalGainWidths> produceGainWidths (const HcalGainWidthsRcd& rcd);
  std::auto_ptr<HcalQIEData> produceQIEData (const HcalQIEDataRcd& rcd);
  std::auto_ptr<HcalChannelQuality> produceChannelQuality (const HcalChannelQualityRcd& rcd);
  std::auto_ptr<HcalElectronicsMap> produceElectronicsMap (const HcalElectronicsMapRcd& rcd);

  std::auto_ptr<HcalRespCorrs> produceRespCorrs (const HcalRespCorrsRcd& rcd);
  std::auto_ptr<HcalZSThresholds> produceZSThresholds (const HcalZSThresholdsRcd& rcd);
  std::auto_ptr<HcalL1TriggerObjects> produceL1TriggerObjects (const HcalL1TriggerObjectsRcd& rcd);
  std::auto_ptr<HcalTimeCorrs> produceTimeCorrs (const HcalTimeCorrsRcd& rcd);
  std::auto_ptr<HcalLUTCorrs> produceLUTCorrs (const HcalLUTCorrsRcd& rcd);
  std::auto_ptr<HcalPFCorrs> producePFCorrs (const HcalPFCorrsRcd& rcd);

  std::auto_ptr<HcalValidationCorrs> produceValidationCorrs (const HcalValidationCorrsRcd& rcd);
  std::auto_ptr<HcalLutMetadata> produceLutMetadata (const HcalLutMetadataRcd& rcd);

  bool h2mode_;
};

