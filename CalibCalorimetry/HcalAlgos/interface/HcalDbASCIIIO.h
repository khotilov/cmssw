//
// F.Ratnikov (UMd), Jul. 19, 2005
//
#ifndef HcalDbASCIIIO_h
#define HcalDbASCIIIO_h

#include <iostream>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CondFormats/HcalObjects/interface/AllClasses.h"

/**
   \class HcalDbASCIIIO
   \brief IO for ASCII instances of Hcal Calibrations
   \author Fedor Ratnikov Oct. 28, 2005
   $Id: HcalDbASCIIIO.h,v 1.3 2005/12/15 23:37:56 fedor Exp $
   
Text file formats for different data types is as following:
- # in first column comments the line
- HcalPedestals, HcalPedestalWidths, HcalGains, HcalGainWidths have identical formats:
  eta(int)  phi(int) depth(int) det(HB,HE,HF) cap1_value(float) cap2_value(float) cap3_value(float) cap4_value(float)  HcalDetId(int,optional)
- HcalQIEShape:
  33 x floats - low edges for first 33 channels of ideal QIE
- HcalQIEData:
  eta phi depth det 4x offsets_cap1 4x offsets_cap2 4x offsets_cap3 4x offsets_cap4 4x slopes_cap1 4x slopes_cap2 4x slopes_cap3 4x slopes_cap4
- HcalChannelQuality:
  eta phi depth det status(GOOD/BAD/HOT/DEAD)
- HcalElectronicsMap:
  line#  crate HTR_slot top_bottom(t/b) dcc# dcc_spigot fiber fiberchan subdet(HB/HE/HF/HO/HT) eta phi depth
  if electronics channel is known to be unconnected, either "subdet" or "eta" should be NA
*/
namespace HcalDbASCIIIO {
  bool getObject (std::istream& fInput, HcalPedestals* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalPedestals& fObject);
  bool getObject (std::istream& fInput, HcalPedestalWidths* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalPedestalWidths& fObject);
  bool getObject (std::istream& fInput, HcalGains* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalGains& fObject);
  bool getObject (std::istream& fInput, HcalGainWidths* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalGainWidths& fObject);
  bool getObject (std::istream& fInput, HcalQIEData* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalQIEData& fObject);
  bool getObject (std::istream& fInput, HcalChannelQuality* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalChannelQuality& fObject);
  bool getObject (std::istream& fInput, HcalElectronicsMap* fObject);
  bool dumpObject (std::ostream& fOutput, const HcalElectronicsMap& fObject);
} 
#endif
