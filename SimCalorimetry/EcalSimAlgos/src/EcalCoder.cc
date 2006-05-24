#include "SimCalorimetry/EcalSimAlgos/interface/EcalCoder.h"
#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CLHEP/Random/RandGaussQ.h"
#include <iostream>


EcalCoder::EcalCoder(bool addNoise)
:  thePedestals(0),
   addNoise_(addNoise)

{

  // 4095(MAXADC)*12(gain 2)*0.035(GeVtoADC)
  
  m_maxEneEB = 1719.9 ; 
  
  // 4095(MAXADC)*12(gain 2)*0.060(GeVtoADC)
  
  m_maxEneEE = 2948.4 ; 
  
}  


double EcalCoder::fullScaleEnergy(const DetId & detId) const 
{

  if (detId.subdetId() == EcalBarrel) 
    return m_maxEneEB ;
  else 
    return m_maxEneEE ;
  
}


void EcalCoder::digitalToAnalog(const EBDataFrame& df, CaloSamples& lf) const {
  for(int i = 0; i < df.size(); ++i) {
    lf[i] = decode(df[i], df.id());
  }
}  


void EcalCoder::digitalToAnalog(const EEDataFrame& df, CaloSamples& lf) const {
  for(int i = 0; i < df.size(); ++i) {
    lf[i] = decode(df[i], df.id());
  }
}


void EcalCoder::analogToDigital(const CaloSamples& clf, EBDataFrame& df) const {
  std::vector<EcalMGPASample> mgpaSamples = encode(clf);

  df.setSize(clf.size());
  for(int i = 0; i < df.size(); ++i) {
    df.setSample(i, mgpaSamples[i]);
  }

}


void EcalCoder::analogToDigital(const CaloSamples& clf, EEDataFrame& df) const {
  std::vector<EcalMGPASample> mgpaSamples = encode(clf);

  df.setSize(clf.size());
  for(int i = 0; i < df.size(); ++i) {
    df.setSample(i, mgpaSamples[i]);
  }
}


std::vector<EcalMGPASample>  
EcalCoder::encode(const CaloSamples& caloSamples) const
{
  assert(thePedestals != 0);
  std::vector<EcalMGPASample> results;
  results.reserve(caloSamples.size());

  DetId detId = caloSamples.id();
  double Emax = fullScaleEnergy(detId);
  //....initialisation

  if ( caloSamples[5] > 0. ) 
    LogDebug("DigiCoder") << "Input caloSample" << "\n" << caloSamples;

  double LSB[NGAINS];
  double pedestals[NGAINS];
  double widths[NGAINS];
  double gains[NGAINS];
  double threeSigmaADCNoise[NGAINS];
  for(int igain = 0; igain < NGAINS; ++igain) {
    // fill in the pedestal and width
    findPedestal(detId, igain, pedestals[igain], widths[igain]);
    // set nominal value first
    findGains(detId, gains);
    LSB[igain]= Emax/(MAXADC*gains[igain]);
    threeSigmaADCNoise[igain] = widths[igain]/LSB[igain] * 3.;
  }

  int wait = 0 ;
  int gainId = NGAINS - 1 ;
  for (int i = 0 ; i < caloSamples.size() ; ++i)
  {    
     int adc = MAXADC;
     if (wait == 0) gainId = NGAINS - 1;

     // see which gain bin it fits in
     int igain = gainId + 1 ;
     while (igain != 0) {
       --igain;

       double ped = pedestals[igain];
       int tmpadc = int (ped + caloSamples[i] / LSB[igain]) ;

       // see if it's close enough to the boundary that we have to throw noise
       if(addNoise_ && (tmpadc <= MAXADC+threeSigmaADCNoise[igain]) ) {
          ped = RandGauss::shoot(ped, widths[igain]);
          tmpadc = int (ped + caloSamples[i] / LSB[igain]) ;
       }
       //std::cout << "DetId " << detId.rawId() << " gain " << igain << " " << caloSamples[i] << " " << pedestals[igain] << " " << widths[igain] << " " << LSB[igain] << " result = " << ped << " " << tmpadc <<std::endl;
         
       if(tmpadc <= MAXADC ) {
         adc = tmpadc;
         break ;
       }
     }
     
     if (igain == NGAINS - 1) 
       {
         wait = 0 ;
         gainId = igain ;
       }
     else 
       {
         if (igain == gainId) --wait ;
         else 
           {
             wait = 5 ;
             gainId = igain ;
           }
       }

     results.push_back(EcalMGPASample(adc, gainId));
  }
  return results;
}

double EcalCoder::decode(const EcalMGPASample & sample, const DetId & id) const
{
  double Emax = fullScaleEnergy(id); 
  int gainNumber  = sample.gainId();
  assert( gainNumber >=0 && gainNumber <=2);
  double gains[NGAINS];
  findGains(id, gains);
  double LSB = Emax/(MAXADC*gains[gainNumber]) ;
  double pedestal = 0.;
  double width = 0.;
  findPedestal(id, gainNumber, pedestal, width);
  // we shift by LSB/2 to be centered
  return LSB * (sample.adc() + 0.5 - pedestal) ;
}

void EcalCoder::findPedestal(const DetId & detId, int gainId, 
                             double & ped, double & width) const
{
  EcalPedestalsMapIterator mapItr 
    = thePedestals->m_pedestals.find(detId.rawId());
  // should I care if it doesn't get found?
  if(mapItr == thePedestals->m_pedestals.end()) {
    edm::LogError("SetupInfo") << "Could not find pedestal for " << detId.rawId() << " among the " << thePedestals->m_pedestals.size();
  } else {
    EcalPedestals::Item item = mapItr->second;

      switch(gainId) {
    case 0:
      ped = item.mean_x1;
      width = item.rms_x1;
      break;
    case 1:
      ped = item.mean_x6;
      width = item.rms_x6;
      break;
    case 2:
      ped = item.mean_x12;
      width = item.rms_x12;
      break;
    default:
      edm::LogError("SetupInfo") << "Bad Pedestal " << gainId;
      break;
    }
    //LogDebug("SetupInfo") << "Pedestals for " << detId.rawId() << " gain range " << gainId << " : \n" << "Mean = " << ped << " rms = " << width;
  }
}

void EcalCoder::findGains(const DetId & detId, double Gains[]) const
{
  EcalGainRatios::EcalGainRatioMap::const_iterator grit=theGainRatios->getMap().find(detId.rawId());
  EcalMGPAGainRatio mgpa;
  if( grit!=theGainRatios->getMap().end() ){
    mgpa = grit->second;
    Gains[0] = 1.;
    Gains[1] = mgpa.gain6Over1() ;
    Gains[2] = Gains[1]*(mgpa.gain12Over6()) ;
    //LogDebug("SetupInfo") << "Gains for " << detId.rawId() << "\n" << " 0 = " << Gains[0] << "\n" << " 1 = " << Gains[1] << "\n" << " 2 = " << Gains[2];
  } else {
    edm::LogError("SetupInfo") << "Could not find gain ratios for " << detId.rawId() << " among the " << theGainRatios->getMap().size();
  }
}
