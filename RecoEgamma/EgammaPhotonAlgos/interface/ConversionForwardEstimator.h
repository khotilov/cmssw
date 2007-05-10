#ifndef  RecoEGAMMA_ConversionForwardEstimator_H
#define  RecoEGAMMA_ConversionForwardEstimator_H

/**
 * \class ConversionForwardEstimator
 *  Defines the search area in the  forward 
 *
 *   $Date: 2007/03/07 16:48:53 $ 
 *   $Revision: 1.3 $
 *  \author Nancy Marinelli, U. of Notre Dame, US
 */

#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h" 
#include "DataFormats/GeometryVector/interface/Vector2DBase.h"
#include "DataFormats/GeometryVector/interface/LocalTag.h"


#include <iostream> 
class RecHit;
class TrajectoryStateOnSurface;
class BoundPlane;

class ConversionForwardEstimator : public MeasurementEstimator {
public:
  ConversionForwardEstimator() {};
  ConversionForwardEstimator( float phiRangeMin, float phiRangeMax, float dr, double nSigma = 3.  ) :
                           thePhiRangeMin( phiRangeMin), thePhiRangeMax( phiRangeMax), dr_(dr),  theNSigma(nSigma) {
    //std::cout << " ConversionForwardEstimator CTOR " << std::endl;
}

  // zero value indicates incompatible ts - hit pair
  virtual std::pair<bool,double> estimate( const TrajectoryStateOnSurface& ts, 
			   const TransientTrackingRecHit& hit) const;
  virtual bool estimate( const TrajectoryStateOnSurface& ts, 
			   const BoundPlane& plane) const;
  virtual ConversionForwardEstimator* clone() const {
    return new ConversionForwardEstimator(*this);
  } 


virtual Local2DVector maximalLocalDisplacement( const TrajectoryStateOnSurface& ts, const BoundPlane& plane) const;

 double nSigmaCut() const {return theNSigma;}

private:

  float thePhiRangeMin;
  float thePhiRangeMax;
  float dr_;
  double theNSigma;
};

#endif // ConversionForwardEstimator_H












