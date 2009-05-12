#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/SourceFactory.h"

DEFINE_SEAL_MODULE();

#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"

//-----------------------------------------------------------------------------------------

#include "CalibTracker/SiStripESProducers/plugins/DBWriter/SiStripFedCablingManipulator.h"
DEFINE_ANOTHER_FWK_MODULE(SiStripFedCablingManipulator);
//-----------------------------------------------------------------------------------------

#include "CalibTracker/SiStripESProducers/plugins/DBWriter/DummyCondDBWriter.h"

#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
struct CabRcdName{ static const char* name(){return "SiStripFedCablingRcd";}};
typedef DummyCondDBWriter<SiStripFedCabling,SiStripFedCabling,SiStripFedCablingRcd,CabRcdName> SiStripFedCablingDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripFedCablingDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
struct PedRcdName{ static const char* name(){return "SiStripPedestalsRcd";}};
typedef DummyCondDBWriter<SiStripPedestals,SiStripPedestals,SiStripPedestalsRcd,PedRcdName> SiStripPedestalsDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripPedestalsDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
struct NoiseRcdName{ static const char* name(){return "SiStripNoisesRcd";}};
typedef DummyCondDBWriter<SiStripNoises,SiStripNoises,SiStripNoisesRcd,NoiseRcdName> SiStripNoisesDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripNoisesDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripApvGain.h"
struct GainRcdName{ static const char* name(){return "SiStripApvGainRcd";}};
typedef DummyCondDBWriter<SiStripApvGain,SiStripApvGain,SiStripApvGainRcd,GainRcdName> SiStripApvGainDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripApvGainDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
struct LARcdName{ static const char* name(){return "SiStripLorentzAngleRcd";}};
typedef DummyCondDBWriter<SiStripLorentzAngle,SiStripLorentzAngle,SiStripLorentzAngleRcd,LARcdName> SiStripLorentzAngleDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripLorentzAngleDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripThreshold.h"
struct ThRcdName{ static const char* name(){return "SiStripThresholdRcd";}};
typedef DummyCondDBWriter<SiStripThreshold,SiStripThreshold,SiStripThresholdRcd,ThRcdName> SiStripThresholdDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripThresholdDummyDBWriter);


#include "CondFormats/SiStripObjects/interface/SiStripBadStrip.h"
struct BadStripRcdName{ static const char* name(){return "SiStripBadStrip";}};
typedef DummyCondDBWriter<SiStripBadStrip,SiStripBadStrip,SiStripBadStripRcd,BadStripRcdName> SiStripBadStripDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripBadStripDummyDBWriter);

typedef DummyCondDBWriter<SiStripBadStrip,SiStripBadStrip,SiStripBadModuleRcd,BadStripRcdName> SiStripBadModuleDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripBadModuleDummyDBWriter);

typedef DummyCondDBWriter<SiStripBadStrip,SiStripBadStrip,SiStripBadFiberRcd,BadStripRcdName> SiStripBadFiberDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripBadFiberDummyDBWriter);

typedef DummyCondDBWriter<SiStripBadStrip,SiStripBadStrip,SiStripBadChannelRcd,BadStripRcdName> SiStripBadChannelDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripBadChannelDummyDBWriter);

#include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
#include "CalibTracker/Records/interface/SiStripQualityRcd.h"
typedef DummyCondDBWriter<SiStripQuality,SiStripBadStrip,SiStripQualityRcd,BadStripRcdName> SiStripBadStripFromQualityDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripBadStripFromQualityDummyDBWriter);

#include "CondFormats/SiStripObjects/interface/SiStripModuleHV.h"
struct ModuleHVRcdName{ static const char* name(){return "SiStripModuleHVRcd";}};
typedef DummyCondDBWriter<SiStripModuleHV,SiStripModuleHV,SiStripModuleHVRcd,ModuleHVRcdName> SiStripModuleHVDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripModuleHVDummyDBWriter);

#include "CondFormats/SiStripObjects/interface/SiStripDetVOff.h"
struct DetVOffRcdName{ static const char* name(){return "SiStripDetVOffRcd";}};
typedef DummyCondDBWriter<SiStripDetVOff,SiStripDetVOff,SiStripDetVOffRcd,DetVOffRcdName> SiStripDetVOffDummyDBWriter;
DEFINE_ANOTHER_FWK_MODULE(SiStripDetVOffDummyDBWriter);


//---------------------------------------------------------------------------------------------------------------
// Dummy printers

#include "CalibTracker/SiStripESProducers/plugins/DBWriter/DummyCondObjPrinter.h"

typedef DummyCondObjPrinter<SiStripNoises,SiStripNoisesRcd> SiStripNoisesDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripNoisesDummyPrinter);

typedef DummyCondObjPrinter<SiStripThreshold,SiStripThresholdRcd> SiStripThresholdDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripThresholdDummyPrinter);

typedef DummyCondObjPrinter<SiStripLorentzAngle,SiStripLorentzAngleRcd> SiStripLorentzAngleDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripLorentzAngleDummyPrinter);

typedef DummyCondObjPrinter<SiStripLorentzAngle,SiStripLorentzAngleSimRcd> SiStripLorentzAngleSimDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripLorentzAngleSimDummyPrinter);

typedef DummyCondObjPrinter<SiStripPedestals,SiStripPedestalsRcd> SiStripPedestalsDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripPedestalsDummyPrinter);

#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
typedef DummyCondObjPrinter<SiStripGain,SiStripGainRcd> SiStripGainDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripGainDummyPrinter);

#include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
typedef DummyCondObjPrinter<SiStripGain,SiStripGainSimRcd> SiStripGainSimDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripGainSimDummyPrinter);

typedef DummyCondObjPrinter<SiStripDetVOff,SiStripDetVOffRcd> SiStripDetVOffDummyPrinter;
DEFINE_ANOTHER_FWK_MODULE(SiStripDetVOffDummyPrinter);
