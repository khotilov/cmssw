#ifndef MeasurementDet_H
#define MeasurementDet_H

#include <vector>

#include "Geometry/CommonDetUnit/interface/GeomDet.h"

class TrajectoryStateOnSurface;
class TrajectoryMeasurement;
class GeomDet;
class Propagator;
class MeasurementEstimator;
class TransientTrackingRecHit;
class BoundSurface;

class MeasurementDet {
public:

  typedef std::vector<TransientTrackingRecHit*>         RecHitContainer;

  MeasurementDet( const GeomDet* gdet) : theGeomDet(gdet) {}

  virtual RecHitContainer recHits( const TrajectoryStateOnSurface&) const = 0;

  /** faster version in case the TrajectoryState on the surface of the
   *  Det is already available. The first TrajectoryStateOnSurface is on the surface of this 
   *  Det, and the second TrajectoryStateOnSurface is the statrting state, usually
   *  not on the surface of this Det. The stateOnThisDet should the result of <BR>
   *  prop.propagate( startingState, this->surface())
   */
  virtual std::vector<TrajectoryMeasurement> 
  fastMeasurements( const TrajectoryStateOnSurface& stateOnThisDet, 
		    const TrajectoryStateOnSurface& startingState, 
		    const Propagator&, 
		    const MeasurementEstimator&) const = 0;
  

  virtual const GeomDet& geomDet() const { return *theGeomDet;}

  const BoundSurface& surface() const {return  geomDet().surface();}

  const Surface::PositionType& position() const {return geomDet().position();}

 private:

  const GeomDet* theGeomDet;

};


#endif
