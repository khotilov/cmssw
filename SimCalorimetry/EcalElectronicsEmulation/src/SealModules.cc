#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimCalorimetry/EcalElectronicsEmulation/interface/EcalSimpleSource.h"
#include "SimCalorimetry/EcalElectronicsEmulation/interface/EcalSimRawData.h"
#include "SimCalorimetry/EcalElectronicsEmulation/interface/EcalFEtoDigi.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_INPUT_SOURCE(EcalSimpleSource);
DEFINE_ANOTHER_FWK_MODULE(EcalSimRawData);
DEFINE_ANOTHER_FWK_MODULE(EcalFEtoDigi);
