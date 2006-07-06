
#ifndef EcalSimAlgos_EcalCoder_h
#define EcalSimAlgos_EcalCoder_h 1

#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloCorrelatedNoisifier.h"

#include "CalibCalorimetry/EcalTrivialCondModules/interface/EcalTrivialConditionRetriever.h"


class EcalMGPASample;
class EcalPedestals;
class EBDataFrame;
class EEDataFrame;
class DetId;
#include<vector>

/* \class EEDigitizerTraits
 * \brief Converts CaloDataFrame in CaloTimeSample and vice versa.
 *
 */
class EcalCoder
{
 public:
  /// number of available bits
  enum {NBITS = 12};
  // 2^12 -1
  /// adc max range
  enum {MAXADC = 4095}; 
  /// adc gain switch
  enum {ADCGAINSWITCH = 4080};
  /// number of electronic gains
  enum {NGAINS = 3};

  /// ctor
  EcalCoder(bool addNoise, CaloCorrelatedNoisifier * theCorrNoise) ;
  /// dtor
  virtual ~EcalCoder() {}

  /// can be fetched every event from the EventSetup
  void setPedestals(const EcalPedestals * pedestals) {thePedestals = pedestals;}

  void setGainRatios(const EcalGainRatios * gainRatios) {theGainRatios = gainRatios; }

  void setFullScaleEnergy(const double EBscale , const double EEscale) {m_maxEneEB = EBscale; m_maxEneEE = EEscale; }

 
  /// from EBDataFrame to CaloSamples
  virtual void digitalToAnalog(const EBDataFrame& df, CaloSamples& lf) const;
  /// from EEDataFrame to CaloSamples
  virtual void digitalToAnalog(const EEDataFrame& df, CaloSamples& lf) const;
  /// from CaloSamples to EBDataFrame
  virtual void analogToDigital(const CaloSamples& clf, EBDataFrame& df) const;
  /// from CaloSamples to EEDataFrame
  virtual void analogToDigital(const CaloSamples& clf, EEDataFrame& df) const;
 
  ///  anything that needs to be done once per event
  void newEvent() {}

 private:

  /// limit on the energy scale due to the electronics range
  double fullScaleEnergy (const DetId & ) const ;

  /// produce the pulse-shape
  std::vector<EcalMGPASample> encode(const CaloSamples& timeframe) const;

  double decode(const EcalMGPASample & sample, const DetId & detId) const;

  /// not yet implemented
  void noisify(float * values, int size) const;

  /// look for the right pedestal according to the electronics gain
  void findPedestal(const DetId & detId, int gainId, 
                    double & pedestal, double & width) const;

  double theGains[NGAINS+1];
   
  void findGains(const DetId & detId, double theGains[] ) const;
   
  /// the pedestals
  const EcalPedestals * thePedestals;
  /// the electronics gains
  const EcalGainRatios * theGainRatios;
  /// max attainable energy in the ecal barrel
  double m_maxEneEB ;
  /// max attainable energy in the ecal endcap
  double m_maxEneEE ;
  /// whether add noise to the pedestals and the gains
  bool addNoise_;
  /// Correlated noisifier
  CaloCorrelatedNoisifier * theCorrNoise_;

};


#endif
