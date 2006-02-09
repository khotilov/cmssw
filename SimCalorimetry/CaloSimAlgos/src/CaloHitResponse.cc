using namespace std;
#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h" 
#include "SimDataFormats/CaloHit/interface/PCaloHit.h" 
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVSimParameterMap.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloSimParameters.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVShape.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVHitCorrection.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloVHitFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include<iostream>


CaloHitResponse::CaloHitResponse(const CaloVSimParameterMap * parametersMap, 
                                 const CaloVShape * shape)
: theParameterMap(parametersMap), 
  theShape(shape),  
  theHitCorrection(0),
  theHitFilter(0),
  theGeometry(0),
  theMinBunch(-10), 
  theMaxBunch(10)
{
}


void CaloHitResponse::setBunchRange(int minBunch, int maxBunch) {
  theMinBunch = minBunch;
  theMaxBunch = maxBunch;
}


void CaloHitResponse::run(MixCollection<PCaloHit> & hits) {
  for(MixCollection<PCaloHit>::MixItr hitItr = hits.begin();
      hitItr != hits.end(); ++hitItr)
  {
    // maybe it's not from this subdetector
    if(theHitFilter == 0 || theHitFilter->accepts(*hitItr)) {
      LogDebug("CaloHitResponse") << *hitItr;
      CaloSamples signal = makeAnalogSignal(*hitItr);
      LogDebug("CaloHitResponse") << signal;
      // if there's already a frame for this in the map, superimpose it
      DetId id(hitItr->id());
      CaloSamples * oldSignal = findSignal(id);
      if (oldSignal == 0) {
        theAnalogSignalMap[id] = signal;
      } else  {
        // need a "+=" to CaloSamples
        int sampleSize =  oldSignal->size();
        assert(sampleSize == signal.size());
        assert(signal.presamples() == oldSignal->presamples());
 
        for(int i = 0; i < sampleSize; ++i) {
          (*oldSignal)[i] += signal[i];
        }
      }
    } //  filter accepts 
  } // loop over hits
}


CaloSamples CaloHitResponse::makeAnalogSignal(const PCaloHit & inputHit) const {

  // see if we need to correct the hit 
  PCaloHit hit = inputHit;
  if(theHitCorrection != 0) {
    theHitCorrection->correct(hit);
  }

  DetId detId(hit.id());
  const CaloSimParameters & parameters = theParameterMap->simParameters(detId);


  double signal = analogSignalAmplitude(hit, parameters);

  double jitter = hit.time() - timeOfFlight(detId);

  // assume bins count from zero, go for center of bin
  const double tzero = parameters.timePhase() -jitter -
     BUNCHSPACE*(parameters.binOfMaximum()-1);
  double binTime = tzero;

  CaloSamples result(detId, parameters.readoutFrameSize());
  for(int bin = 0; bin < result.size(); bin++) {
    result[bin] += (*theShape)(binTime)* signal;
    binTime += BUNCHSPACE;
  }

  result.setPresamples(parameters.binOfMaximum()-1);
  return result;
} 


double CaloHitResponse::analogSignalAmplitude(const PCaloHit & hit, const CaloSimParameters & parameters) const {
  // OK, the "energy" in the hit could be a real energy, deposited energy,
  // or pe count.  This factor converts to photoelectrons
  double npe = hit.energy() * parameters.simHitToPhotoelectrons();
  // do we need to doPoisson statistics for the photoelectrons?
  if(parameters.doPhotostatistics()) {
    npe = RandPoisson::shoot(static_cast<int>(npe));
  }
  return npe;
}


CaloSamples * CaloHitResponse::findSignal(const DetId & detId) {
  CaloSamples * result = 0;
  AnalogSignalMap::iterator signalItr = theAnalogSignalMap.find(detId);
  if(signalItr == theAnalogSignalMap.end()) {
     result = 0;
  } 
  else {
     result = &(signalItr->second);
  }
  return result;
}


CaloSamples * CaloHitResponse::makeNewSignal(const DetId & detId) {
  const CaloSimParameters & parameters = theParameterMap->simParameters(detId);
  CaloSamples result(detId, parameters.readoutFrameSize());
  result.setPresamples(parameters.binOfMaximum()-1);
  theAnalogSignalMap[detId] = result;
  return &(theAnalogSignalMap[detId]);
}


double CaloHitResponse::timeOfFlight(const DetId & detId) const {
  // not going to assume there's one of these per subdetector.
  // Take the whole CaloGeometry and find the right subdet
  double result = 0.;
  if(theGeometry == 0) {
    edm::LogWarning("CaloHitResponse") << "No Calo Geometry set, so no time of flight correction";
  } 
  else {
    const CaloCellGeometry* cellGeometry = theGeometry->getSubdetectorGeometry(detId)->getGeometry(detId);
    double distance = cellGeometry->getPosition().mag();
    result =  distance * cm / c_light; // Units of c_light: mm/ns
  }
  return result;
}


