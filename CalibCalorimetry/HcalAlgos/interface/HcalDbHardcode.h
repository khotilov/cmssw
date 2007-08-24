//
// F.Ratnikov (UMd), Dec. 14, 2005
//
#ifndef HcalDbHardcodeIn_h
#define HcalDbHardcodeIn_h

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"
#include "CondFormats/HcalObjects/interface/HcalGain.h"
#include "CondFormats/HcalObjects/interface/HcalGainWidth.h"
#include "CondFormats/HcalObjects/interface/HcalQIECoder.h"
#include "CondFormats/HcalObjects/interface/HcalQIEShape.h"
#include "CondFormats/HcalObjects/interface/HcalCalibrationQIECoder.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"


/**

   \class HcalDbHardcode
   \brief Hardcode implementation of some conditions data
   \author Fedor Ratnikov
   
*/
namespace HcalDbHardcode {
  HcalPedestal makePedestal (HcalGenericDetId fId, bool fSmear = false);
  HcalPedestalWidth makePedestalWidth (HcalGenericDetId fId);
  HcalGain makeGain (HcalGenericDetId fId, bool fSmear = false);
  HcalGainWidth makeGainWidth (HcalGenericDetId fId);
  HcalQIECoder makeQIECoder (HcalGenericDetId fId);
  HcalCalibrationQIECoder makeCalibrationQIECoder (HcalGenericDetId fId);
  HcalQIEShape makeQIEShape ();
  void makeHardcodeMap(HcalElectronicsMap& emap);
}
#endif
