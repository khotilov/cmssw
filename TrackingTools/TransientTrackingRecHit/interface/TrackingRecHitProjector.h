#ifndef TrackingRecHitProjector_H
#define TrackingRecHitProjector_H

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
//#include <iostream>

template <class ResultingHit>
class TrackingRecHitProjector {
 public:

  typedef  TransientTrackingRecHit::RecHitPointer      RecHitPointer;

  RecHitPointer project( const TransientTrackingRecHit& hit,
			 const GeomDet& det, 
			 const TrajectoryStateOnSurface& ts) const {
    using namespace std;

    GlobalVector gdir = ts.globalParameters().momentum();
    const BoundPlane& gluedPlane = det.surface();
    const BoundPlane& hitPlane = hit.det()->surface();

    // check if the planes are parallel
    const float epsilon = 1.e-7; // corresponds to about 0.3 miliradian but cannot be reduced
                                 // because of float precision
    if (fabs(gluedPlane.normalVector().dot( hitPlane.normalVector())) < 1-epsilon) {
//       cout << "TkGluedMeasurementDet plane not parallel to DetUnit plane: dot product is " 
// 	   << gluedPlane.normalVector().dot( hitPlane.normalVector()) << endl;
// FIXME: throw the appropriate exception here...
      //throw MeasurementDetException("TkGluedMeasurementDet plane not parallel to DetUnit plane");
    }

    double delta = gluedPlane.localZ( hitPlane.position());
    LocalVector ldir = gluedPlane.toLocal(gdir);
    LocalPoint lhitPos = gluedPlane.toLocal( hit.globalPosition());
    LocalPoint projectedHitPos = lhitPos - ldir * delta/ldir.z();

    LocalVector hitXAxis = gluedPlane.toLocal( hitPlane.toGlobal( LocalVector(1,0,0)));
    LocalError hitErr = hit.localPositionError();
    if (gluedPlane.normalVector().dot( hitPlane.normalVector()) < 0) {
      // the two planes are inverted, and the correlation element must change sign
      hitErr = LocalError( hitErr.xx(), -hitErr.xy(), hitErr.yy());
    }
    LocalError rotatedError = hitErr.rotate( hitXAxis.x(), hitXAxis.y());

    return ResultingHit::build( projectedHitPos, rotatedError, &det, hit);
  }

};

#endif
