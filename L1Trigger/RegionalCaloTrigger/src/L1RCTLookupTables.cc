#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
#include <string>

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/RegionalCaloTrigger/interface/L1RCTLookupTables.h"

#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/L1TObjects/interface/L1RCTChannelMask.h"
#include "CondFormats/L1TObjects/interface/L1CaloEcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"

unsigned int L1RCTLookupTables::lookup(unsigned short ecalInput,
				       unsigned short hcalInput,
				       unsigned short fgbit,
				       unsigned short crtNo,
				       unsigned short crdNo,
				       unsigned short twrNo) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  if(channelMask_ == 0)
    throw cms::Exception("L1RCTChannelMask Invalid")
      << "L1RCTChannelMask should be set every event" << channelMask_;
  if(ecalInput > 0xFF) 
    throw cms::Exception("Invalid Data") 
      << "ECAL compressedET should be less than 0xFF, is " << ecalInput;
  if(hcalInput > 0xFF) 
    throw cms::Exception("Invalid Data") 
      << "HCAL compressedET should be less than 0xFF, is " << hcalInput;
  if(fgbit > 1) 
    throw cms::Exception("Invalid Data") 
      << "ECAL finegrain should be a single bit, is " << fgbit;
  short iEta = (short) rctParameters_->calcIEta(crtNo, crdNo, twrNo);
  unsigned short iAbsEta = (unsigned short) abs(iEta);
  short sign = iEta/iAbsEta;
  unsigned short iPhi = rctParameters_->calcIPhi(crtNo, crdNo, twrNo);
  unsigned short phiSide = (iPhi/4)%2;
  if(iAbsEta < 1 || iAbsEta > 28) 
    throw cms::Exception("Invalid Data") 
      << "1 <= |IEta| <= 28, is " << iAbsEta;
  float ecal;
  float hcal;
  // using channel mask to mask off ecal channels
  if (channelMask_->ecalMask[crtNo][phiSide][iAbsEta])
    {
      ecal = 0;
    }
  else
    {
      ecal = convertEcal(ecalInput, iAbsEta, sign);
    }
  // masking off hcal for channels in channel mask
  if (channelMask_->hcalMask[crtNo][phiSide][iAbsEta])
    {
      hcal = 0;
    }
  else
    {
      hcal = convertHcal(hcalInput, iAbsEta, sign);
    }
  // couts!
  //std::cout << "LUTs: ecalInput=" << ecalInput << " ecalConverted="
  //	    << ecal << std::endl;
  unsigned long etIn7Bits;
  unsigned long etIn9Bits;
  // Saturated input towers cause tower ET pegging at the highest value
  /*if(ecalInput == 0xFF || hcalInput == 0xFF)
    {
      etIn7Bits = 0x7F;
      etIn9Bits = 0x1FF;
      }*/
  /*else*/ if((ecalInput == 0 && hcalInput > 0) &&
  	  ((rctParameters_->noiseVetoHB() && iAbsEta > 0 && iAbsEta < 18)
  	   || (rctParameters_->noiseVetoHEplus() && iAbsEta>17 && crtNo>8)
  	   || (rctParameters_->noiseVetoHEminus() && iAbsEta>17 && crtNo<9)))
   {
      etIn7Bits = 0;
      etIn9Bits = 0;
    }
  else
    {
      etIn7Bits = eGammaETCode(ecal, hcal, iAbsEta);
      etIn9Bits = jetMETETCode(ecal, hcal, iAbsEta);
    }
  // Saturated input towers cause tower ET pegging at the highest value
  if((ecalInput == 0xFF && 
      rctParameters_->eGammaECalScaleFactors()[iAbsEta-1] != 0. ) 
     || (hcalInput == 0xFF &&
	 rctParameters_->eGammaHCalScaleFactors()[iAbsEta-1] != 0. )
     )
    {
      etIn7Bits = 0x7F; // egamma path
    }
  if((ecalInput == 0xFF &&
      rctParameters_->jetMETECalScaleFactors()[iAbsEta-1] != 0. )
     || (hcalInput == 0xFF &&
	 rctParameters_->jetMETHCalScaleFactors()[iAbsEta-1] != 0. ))
    {
      etIn9Bits = 0x1FF; // sums path
    }

  unsigned long shiftEtIn9Bits = etIn9Bits<<8;
  unsigned long shiftHE_FGBit = hOeFGVetoBit(ecal, hcal, fgbit)<<7;
  unsigned long shiftActivityBit = 0;
  if ( rctParameters_->jetMETECalScaleFactors()[iAbsEta-1] == 0.
       && rctParameters_->jetMETHCalScaleFactors()[iAbsEta-1] == 0. )
    {
      // do nothing, it's already zero
    }
  else if (rctParameters_->jetMETECalScaleFactors()[iAbsEta-1] == 0. )
    {
      shiftActivityBit = activityBit(0., hcal)<<17;
    }
  else if (rctParameters_->jetMETHCalScaleFactors()[iAbsEta-1] == 0. )
    {
      shiftActivityBit = activityBit(ecal, 0.)<<17;
    }
  else
    {
      shiftActivityBit = activityBit(ecal, hcal)<<17;
    }
  unsigned long output=etIn7Bits+shiftHE_FGBit+shiftEtIn9Bits+shiftActivityBit;
  return output;
}

unsigned int L1RCTLookupTables::lookup(unsigned short hfInput,
				       unsigned short crtNo,
				       unsigned short crdNo,
				       unsigned short twrNo
				       ) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  if(channelMask_ == 0)
    throw cms::Exception("L1RCTChannelMask Invalid")
      << "L1RCTChannelMask should be set every event" << channelMask_;
  if(hfInput > 0xFF) 
    throw cms::Exception("Invalid Data") 
      << "HF compressedET should be less than 0xFF, is " << hfInput;
  short iEta = rctParameters_->calcIEta(crtNo, crdNo, twrNo);
  unsigned short iAbsEta = abs(iEta);
  short sign = (iEta/iAbsEta);
  unsigned short phiSide = twrNo/4;
  if(iAbsEta < 29 || iAbsEta > 32) 
    throw cms::Exception("Invalid Data") 
      << "29 <= |iEta| <= 32, is " << iAbsEta;
  float et;
  if (channelMask_->hfMask[crtNo][phiSide][iAbsEta])
    {
      et = 0;
    }
  else
    {
      et = convertHcal(hfInput, iAbsEta, sign);
    }
  return convertToInteger(et, rctParameters_->jetMETLSB(), 8);
}

bool L1RCTLookupTables::hOeFGVetoBit(float ecal, float hcal, bool fgbit) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  bool veto = false;
  if(ecal > rctParameters_->eMinForFGCut() && 
     ecal < rctParameters_->eMaxForFGCut())
    {
      if(fgbit) veto = true;
    }
  if(ecal > rctParameters_->eMinForHoECut() && 
     ecal < rctParameters_->eMaxForHoECut())
    {
      if((hcal / ecal) > rctParameters_->hOeCut()) veto = true;
    }
  else 
    {
      if(hcal > rctParameters_->hMinForHoECut()) veto = true;  // Changed from eMinForHoECut() - JLL 2008-Feb-13
    }
  return veto;
}

bool L1RCTLookupTables::activityBit(float ecal, float hcal) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  return ((ecal > rctParameters_->eActivityCut()) || 
	  (hcal > rctParameters_->hActivityCut()));
}

// uses etScale
unsigned int L1RCTLookupTables::emRank(unsigned short energy) const 
{
  if(etScale_)
    {
      return etScale_->rank(energy);
    }
  else
    //    edm::LogInfo("L1RegionalCaloTrigger") 
    //      << "CaloEtScale was not used - energy instead of rank" << endl;
  return energy;
}

// converts compressed ecal energy to linear (real) scale
float L1RCTLookupTables::convertEcal(unsigned short ecal, unsigned short iAbsEta, short sign) const
{
  if(ecalScale_)
    {
      //std::cout << "[luts] energy " << ecal << " sign " << sign 
      //<< " iAbsEta " << iAbsEta << " iPhi "	<< iPhi << std::endl;
      float dummy = 0;
      dummy = float (ecalScale_->et( ecal, iAbsEta, sign ));
      /*
      if (ecal > 0)
	{
	  std::cout << "[luts] ecal converted from " << ecal << " to " 
		    << dummy << " with iAbsEta " << iAbsEta << std::endl;
	}
      */
      return dummy;
    }
  //else if(rctParameters_ == 0)
  //  {
  //    throw cms::Exception("L1RCTParameters Invalid")
  //	<< "L1RCTParameters should be set every event" << rctParameters_;
  //  }
  else
    {
      return ((float) ecal) * rctParameters_->eGammaLSB();
    }
}

// converts compressed hcal energy to linear (real) scale
float L1RCTLookupTables::convertHcal(unsigned short hcal, unsigned short iAbsEta, short sign) const
{
  if (hcalScale_ != 0)
    {
      return (hcalScale_->et( hcal, iAbsEta, sign ));
    }
  else
    {
      //      edm::LogInfo("L1RegionalCaloTrigger") 
      //	<< "CaloTPGTranscoder was not used" << std::endl;
      return ((float) hcal) * rctParameters_->jetMETLSB();
    }
}

// integerize given an LSB and set maximum value of 2^precision-1
unsigned long L1RCTLookupTables::convertToInteger(float et, 
						  float lsb, 
						  int precision) const
{
  unsigned long etBits = (unsigned long)(et/lsb);
  unsigned long maxValue = (1 << precision) - 1;
  if(etBits > maxValue)
    return maxValue;
  else
    return etBits;
}

unsigned int L1RCTLookupTables::eGammaETCode(float ecal, float hcal, int iAbsEta) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  float etLinear = 
    rctParameters_->eGammaECalScaleFactors()[iAbsEta-1] * ecal +
    rctParameters_->eGammaHCalScaleFactors()[iAbsEta-1] * hcal;
  return convertToInteger(etLinear, rctParameters_->eGammaLSB(), 7);
}

unsigned int L1RCTLookupTables::jetMETETCode(float ecal, float hcal, int iAbsEta) const
{
  if(rctParameters_ == 0)
    throw cms::Exception("L1RCTParameters Invalid")
      << "L1RCTParameters should be set every event" << rctParameters_;
  float etLinear = 
    rctParameters_->jetMETECalScaleFactors()[iAbsEta-1] * ecal +
    rctParameters_->jetMETHCalScaleFactors()[iAbsEta-1] * hcal;
  return convertToInteger(etLinear, rctParameters_->jetMETLSB(), 9);
}
