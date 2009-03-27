#ifndef CondFormats_SiStripCondDataRecords_h
#define CondFormats_SiStripCondDataRecords_h

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"

/*Recod associated to SiStripApvGain Object: the SimRcd is used in simulation only*/
class SiStripApvGainRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripApvGainRcd> {};
class SiStripApvGainSimRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripApvGainSimRcd> {};

/*Record associated to SiStripBadStrip Object*/
class SiStripBadChannelRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripBadChannelRcd> {};
class SiStripBadFiberRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripBadFiberRcd> {};
class SiStripBadModuleRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripBadModuleRcd> {};
class SiStripBadStripRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripBadStripRcd> {};
class SiStripDCSStatusRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripDCSStatusRcd> {};

class SiStripFedCablingRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripFedCablingRcd> {};

/*Recod associated to SiStripLorenzaAngle Object: the SimRcd is used in simulation only*/
class SiStripLorentzAngleRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripLorentzAngleRcd> {};
class SiStripLorentzAngleSimRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripLorentzAngleSimRcd> {};

/*Record associated to SiStripModuleHV object*/
class SiStripModuleHVRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripModuleHVRcd> {};
class SiStripModuleLVRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripModuleLVRcd> {};

/*Record associated to SiStripDetVOff object*/
class SiStripDetVOffRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripDetVOffRcd> {};

class SiStripNoisesRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripNoisesRcd> {};

class SiStripPedestalsRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripPedestalsRcd> {};

class SiStripPerformanceSummaryRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripPerformanceSummaryRcd> {};

class SiStripRunSummaryRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripRunSummaryRcd> {};

class SiStripSummaryRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripSummaryRcd> {};

/*Record Associated to SiStripThreshold Object*/
class SiStripThresholdRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripThresholdRcd> {};
class SiStripClusterThresholdRcd : public edm::eventsetup::EventSetupRecordImplementation<SiStripClusterThresholdRcd> {};


#endif
