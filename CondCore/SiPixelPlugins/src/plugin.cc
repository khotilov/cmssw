#include "CondCore/PluginSystem/interface/registration_macros.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelFedCablingMap.h"
#include "CondFormats/DataRecord/interface/SiPixelFedCablingMapRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelGainCalibration.h"
#include "CondFormats/DataRecord/interface/SiPixelGainCalibrationRcd.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleRcd.h"
#include "CondFormats/SiPixelObjects/interface/PixelCPEParmErrors.h"
#include "CondFormats/DataRecord/interface/PixelCPEParmErrorsRcd.h"
#include "CalibFormats/SiPixelObjects/interface/SiPixelCalibConfiguration.h"
#include "CondFormats/DataRecord/interface/SiPixelCalibConfigurationRcd.h"

DEFINE_SEAL_MODULE();
REGISTER_PLUGIN(SiPixelFedCablingMapRcd,SiPixelFedCablingMap);
REGISTER_PLUGIN(SiPixelGainCalibrationRcd,SiPixelGainCalibration);
REGISTER_PLUGIN(SiPixelLorentzAngleRcd,SiPixelLorentzAngle);
REGISTER_PLUGIN(PixelCPEParmErrorsRcd,PixelCPEParmErrors);
REGISTER_PLUGIN(SiPixelCalibConfigurationRcd,SiPixelCalibConfiguration);

