/*
 *  plugin.cc
 *  CMSSW
 *
 *  Created by Chris Jones on 7/24/05.
 *
 */

#include "CondCore/ESSources/interface/registration_macros.h"
#include "CondFormats/CSCObjects/interface/CSCPedestals.h"
#include "CondFormats/CSCObjects/interface/CSCDBPedestals.h"
#include "CondFormats/DataRecord/interface/CSCPedestalsRcd.h"
#include "CondFormats/DataRecord/interface/CSCDBPedestalsRcd.h"
#include "CondFormats/CSCObjects/interface/CSCGains.h"
#include "CondFormats/CSCObjects/interface/CSCDBGains.h"
#include "CondFormats/DataRecord/interface/CSCGainsRcd.h"
#include "CondFormats/DataRecord/interface/CSCDBGainsRcd.h"
#include "CondFormats/CSCObjects/interface/CSCcrosstalk.h"
#include "CondFormats/CSCObjects/interface/CSCDBCrosstalk.h"
#include "CondFormats/DataRecord/interface/CSCcrosstalkRcd.h"
#include "CondFormats/DataRecord/interface/CSCDBCrosstalkRcd.h"
#include "CondFormats/CSCObjects/interface/CSCIdentifier.h"
#include "CondFormats/DataRecord/interface/CSCIdentifierRcd.h"
#include "CondFormats/CSCObjects/interface/CSCNoiseMatrix.h"
#include "CondFormats/CSCObjects/interface/CSCDBNoiseMatrix.h"
#include "CondFormats/DataRecord/interface/CSCNoiseMatrixRcd.h"
#include "CondFormats/DataRecord/interface/CSCDBNoiseMatrixRcd.h"
#include "CondFormats/CSCObjects/interface/CSCChamberMap.h"
#include "CondFormats/DataRecord/interface/CSCChamberMapRcd.h"
#include "CondFormats/CSCObjects/interface/CSCChamberIndex.h"
#include "CondFormats/DataRecord/interface/CSCChamberIndexRcd.h"
#include "CondFormats/CSCObjects/interface/CSCCrateMap.h"
#include "CondFormats/DataRecord/interface/CSCCrateMapRcd.h"
#include "CondFormats/CSCObjects/interface/CSCDDUMap.h"
#include "CondFormats/DataRecord/interface/CSCDDUMapRcd.h"
#include "CondFormats/CSCObjects/interface/CSCBadChambers.h"
#include "CondFormats/DataRecord/interface/CSCBadChambersRcd.h"
#include "CondFormats/CSCObjects/interface/CSCBadStrips.h"
#include "CondFormats/DataRecord/interface/CSCBadStripsRcd.h"
#include "CondFormats/CSCObjects/interface/CSCBadWires.h"
#include "CondFormats/DataRecord/interface/CSCBadWiresRcd.h"
#include "CondFormats/CSCObjects/interface/CSCReadoutMapping.h"
#include "CondFormats/DataRecord/interface/CSCReadoutMappingRcd.h"
#include "CondFormats/CSCObjects/interface/CSCL1TPParameters.h"
#include "CondFormats/DataRecord/interface/CSCL1TPParametersRcd.h"
#include "CondFormats/CSCObjects/interface/CSCDQM_DCSData.h"
#include "CondFormats/DataRecord/interface/CSCDCSDataRcd.h"


REGISTER_PLUGIN(CSCPedestalsRcd,CSCPedestals);
REGISTER_PLUGIN(CSCDBPedestalsRcd,CSCDBPedestals);
REGISTER_PLUGIN(CSCGainsRcd,CSCGains);
REGISTER_PLUGIN(CSCDBGainsRcd,CSCDBGains);
REGISTER_PLUGIN(CSCcrosstalkRcd,CSCcrosstalk);
REGISTER_PLUGIN(CSCDBCrosstalkRcd,CSCDBCrosstalk);
REGISTER_PLUGIN(CSCNoiseMatrixRcd,CSCNoiseMatrix);
REGISTER_PLUGIN(CSCDBNoiseMatrixRcd,CSCDBNoiseMatrix);
REGISTER_PLUGIN(CSCChamberMapRcd,CSCChamberMap);
REGISTER_PLUGIN(CSCChamberIndexRcd,CSCChamberIndex);
REGISTER_PLUGIN(CSCCrateMapRcd,CSCCrateMap);
REGISTER_PLUGIN(CSCDDUMapRcd,CSCDDUMap);
REGISTER_PLUGIN(CSCBadChambersRcd,CSCBadChambers);
REGISTER_PLUGIN(CSCBadStripsRcd,CSCBadStrips);
REGISTER_PLUGIN(CSCBadWiresRcd,CSCBadWires);
REGISTER_PLUGIN(CSCIdentifierRcd,CSCIdentifier);
REGISTER_PLUGIN(CSCReadoutMappingRcd,CSCReadoutMapping);
REGISTER_PLUGIN(CSCL1TPParametersRcd,CSCL1TPParameters);
REGISTER_PLUGIN(CSCDCSDataRcd,  cscdqm::DCSData);
