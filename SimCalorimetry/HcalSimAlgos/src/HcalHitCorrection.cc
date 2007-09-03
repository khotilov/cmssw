#include "SimCalorimetry/HcalSimAlgos/interface/HcalHitCorrection.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloSimParameters.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

HcalHitCorrection::HcalHitCorrection(const CaloVSimParameterMap * parameterMap)
: theParameterMap(parameterMap)
{
}


void HcalHitCorrection::fillChargeSums(MixCollection<PCaloHit> & hits)
{
  clear();
  for(MixCollection<PCaloHit>::MixItr hitItr = hits.begin();
      hitItr != hits.end(); ++hitItr)
  {
    LogDebug("HcalHitCorrection") << "HcalHitCorrection::Hit 0x" << std::hex
				  << hitItr->id() << std::dec;
    int tbin = timeBin(*hitItr);
    LogDebug("HcalHitCorrection") << "HcalHitCorrection::Hit tbin" << tbin;
    if(tbin >= 0 && tbin < 10) 
    {  
      theChargeSumsForTimeBin[tbin][DetId(hitItr->id())] += charge(*hitItr);
    }
  }
}


void HcalHitCorrection::clear()
{
  for(int i = 0; i < 10; ++i)
  {
    theChargeSumsForTimeBin[i].clear();
  }
}

double HcalHitCorrection::charge(const PCaloHit & hit) const
{
  DetId detId(hit.id());
  const CaloSimParameters & parameters = theParameterMap->simParameters(detId);
  double simHitToCharge = parameters.simHitToPhotoelectrons()
                        * parameters.photoelectronsToAnalog();
  return hit.energy() * simHitToCharge;
}


double HcalHitCorrection::delay(const PCaloHit & hit) const 
{
  // HO goes slow, HF shouldn't be used at all
  //ZDC not used for the moment

  DetId detId(hit.id());
  if(detId.det()==DetId::Calo && detId.subdetId()==HcalZDCDetId::SubdetectorId) return 0;
  HcalDetId hcalDetId(hit.id());
  if(hcalDetId.subdet() == HcalForward) return 0;  
  HcalTimeSlew::BiasSetting biasSetting = (hcalDetId.subdet() == HcalOuter) ?
                                          HcalTimeSlew::Slow :
                                          HcalTimeSlew::Medium;
  double delay = 0.;
  int tbin = timeBin(hit);
  if(tbin >= 0 && tbin < 10)
  {
    ChargeSumsByChannel::const_iterator totalChargeItr = theChargeSumsForTimeBin[tbin].find(detId);
    if(totalChargeItr == theChargeSumsForTimeBin[tbin].end())
    {
      throw cms::Exception("HcalHitCorrection") << "Cannot find HCAL charge sum for hit " << hit;
    }
    double totalCharge = totalChargeItr->second;
    delay = HcalTimeSlew::delay(totalCharge, biasSetting);
    LogDebug("HcalHitCorrection") << "TIMESLEWcharge " << charge(hit) 
				  << "  totalcharge " << totalCharge 
				  << " olddelay "  << HcalTimeSlew::delay(charge(hit), biasSetting) 
				  << " newdelay " << delay;
  }

  return delay;
}


void HcalHitCorrection::correct(PCaloHit & hit) const {
  // replace the hit with a new one, with a time delay
  hit = PCaloHit(hit.id(), hit.energyEM(), hit.energyHad(), hit.time()+delay(hit), hit.geantTrackId());
}


int HcalHitCorrection::timeBin(const PCaloHit & hit) const
{
  const CaloSimParameters & parameters = theParameterMap->simParameters(DetId(hit.id()));
  double t = hit.time() - timeOfFlight(DetId(hit.id())) + parameters.timePhase();
  return static_cast<int> (t / 25) + parameters.binOfMaximum() - 1;
}


double HcalHitCorrection::timeOfFlight(const DetId & detId) const
{
  if(detId.det()==DetId::Calo && detId.subdetId()==HcalZDCDetId::SubdetectorId)
    return 37.666; 
  switch(detId.subdetId())
    {
    case HcalBarrel:
      return 8.4;
      break;
    case HcalEndcap:
      return 13.;
      break;
    case HcalOuter:
      return 18.7;
      break;
    case HcalForward:
      return 37.;
      break;
    default:
      throw cms::Exception("HcalHitCorrection") << "Bad Hcal subdetector " << detId.subdetId();
      break;
    }
}

