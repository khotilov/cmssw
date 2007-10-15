#ifndef HcalSimAlgos_HcalSimParameters_h
#define HcalSimAlgos_HcalSimParameters_h

#include "SimCalorimetry/CaloSimAlgos/interface/CaloSimParameters.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

class HcalSimParameters : public CaloSimParameters
{
public:
  HcalSimParameters(double simHitToPhotoelectrons, double photoelectronsToAnalog,
                 double samplingFactor, double timePhase,
                 int readoutFrameSize, int binOfMaximum,
                 bool doPhotostatistics, bool syncPhase,
                 int firstRing, const std::vector<double> & samplingFactors);
  HcalSimParameters(const edm::ParameterSet & p);

  virtual ~HcalSimParameters() {}

  void setDbService(const HcalDbService * service) {theDbService = service;}

  virtual double simHitToPhotoelectrons(const DetId & detId) const;
  virtual double photoelectronsToAnalog(const DetId & detId) const {
    return CaloSimParameters::photoelectronsToAnalog();
  }

  double fCtoGeV(const DetId & detId) const;

  /// the ratio of actual incident energy to deposited energy
  /// in the SimHit
  virtual double samplingFactor(const DetId & detId) const;


private:
  const HcalDbService * theDbService;
  int theFirstRing;
  std::vector<double> theSamplingFactors;
};

#endif
  
