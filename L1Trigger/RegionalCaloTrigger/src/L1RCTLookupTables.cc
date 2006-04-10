#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTLookupTables.h"

unsigned short L1RCTLookupTables::lookup(unsigned short hfenergy){
  float energy = (float)hfenergy*0.5;
  return convertTo10Bits(energy);
}

unsigned long L1RCTLookupTables::lookup(unsigned short ecal,unsigned short hcal,
					 unsigned short fgbit){
  float ecalLinear = convertEcal(ecal);
  float hcalLinear = convertHcal(hcal);
  unsigned long HE_FGBit = (calcHEBit(ecalLinear,hcalLinear) || fgbit);
  unsigned long etIn7Bits = convertTo7Bits(ecalLinear+hcalLinear);
  unsigned long etIn9Bits = convertTo9Bits(ecalLinear+hcalLinear);
  unsigned long activityBit = calcActivityBit(ecalLinear,hcalLinear);

  unsigned long shiftEtIn9Bits = etIn9Bits<<7;
  unsigned long shiftHE_FGBit = HE_FGBit<<16;
  unsigned long shiftActivityBit = activityBit<<17;
  unsigned long output=etIn7Bits+shiftEtIn9Bits+shiftHE_FGBit+shiftActivityBit;
  return output;
}

float L1RCTLookupTables::convertEcal(unsigned short ecal){
  return (float)ecal*0.5;
}

float L1RCTLookupTables::convertHcal(unsigned short hcal){
  return (float)hcal*0.5;
}

unsigned short L1RCTLookupTables::calcActivityBit(float ecal, float hcal){
  return ((ecal > 2) || (hcal > 4));
}

unsigned short L1RCTLookupTables::calcHEBit(float ecal, float hcal){
  return (hcal > ecal);
}

unsigned long L1RCTLookupTables::convertTo7Bits(float et){
  unsigned long etBits = (unsigned long)(et/0.5);
  unsigned long sevenBits = pow(2,7)-1;
  if(etBits > sevenBits)
    return sevenBits;
  else
    return etBits;
}

unsigned long L1RCTLookupTables::convertTo9Bits(float et){
  unsigned long etBits = (unsigned short)(et/0.5);
  unsigned long nineBits = pow(2,9)-1;
  if(etBits > nineBits)
    return nineBits;
  else
    return etBits;
}

unsigned long L1RCTLookupTables::convertTo10Bits(float et){
  unsigned long etBits = (unsigned short)(et/0.5);
  unsigned long tenBits = pow(2,10)-1;
  if(etBits > tenBits)
    return tenBits;
  else
    return etBits;
}
