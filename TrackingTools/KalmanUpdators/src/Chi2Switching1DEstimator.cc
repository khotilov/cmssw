#include "TrackingTools/KalmanUpdators/interface/Chi2Switching1DEstimator.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"

using namespace std;

pair<bool,double> 
Chi2Switching1DEstimator::estimate (const TrajectoryStateOnSurface& aTsos,
				    const TransientTrackingRecHit& aHit) const {
  if(//aHit.isMatched() || 
     aHit.detUnit()->type().isTrackerPixel()) {
    return localEstimator().estimate(aTsos, aHit);
  } else {
    return stripEstimator().estimate(aTsos, aHit);
  }    
}
