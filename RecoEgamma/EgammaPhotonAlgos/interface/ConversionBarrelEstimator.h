#ifndef RecoEGAMMA_ConversionBarrelEstimator_H
#define RecoEGAMMA_ConversionBarrelEstimator_H
/**
 * \class ConversionBarrelEstimator
 *  Defines the search area in the barrel 
 *
 *   $Date: 2006/11/14 11:54:10 $
 *   $Revision: 1.2 $
 *   \author Nancy Marinelli, U. of Notre Dame, US
 */

#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h" 
#include "Geometry/Vector/interface/Vector2DBase.h"
#include "Geometry/Vector/interface/LocalTag.h"


class TrajectoryStateOnSurface;
class RecHit;
class BoundPlane;

class ConversionBarrelEstimator : public MeasurementEstimator {
public:

  ConversionBarrelEstimator() {};
  ConversionBarrelEstimator( float phiRangeMin, float phiRangeMax, 
                                 float zRangeMin, float zRangeMax, double nSigma = 3. ) : 
                           thePhiRangeMin( phiRangeMin), thePhiRangeMax( phiRangeMax),
                           theZRangeMin( zRangeMin), theZRangeMax( zRangeMax), theNSigma(nSigma) {
    //    std::cout << " ConversionBarrelEstimator CTOR " << std::endl;
}

  // zero value indicates incompatible ts - hit pair
  virtual std::pair<bool,double> estimate( const TrajectoryStateOnSurface& ts, 
                               const TransientTrackingRecHit& hit	) const;
  virtual bool  estimate( const TrajectoryStateOnSurface& ts, 
				       const BoundPlane& plane) const;
  virtual ConversionBarrelEstimator* clone() const {
    return new ConversionBarrelEstimator(*this);
  } 





  virtual Local2DVector maximalLocalDisplacement( const TrajectoryStateOnSurface& ts, const BoundPlane& plane) const;


  double nSigmaCut() const {return theNSigma;}

private:

  float thePhiRangeMin;
  float thePhiRangeMax;
  float theZRangeMin;
  float theZRangeMax;
  double theNSigma;



};

#endif // ConversionBarrelEstimator_H
