#ifndef DetLayers_DetRod_H
#define DetLayers_DetRod_H

/** \class DetRod
 *  Abstract interface for a rod of detectors sitting on a BoundPlane.
 */

#include "TrackingTools/DetLayers/interface/CompositeGSD.h"
#include "Geometry/Surface/interface/BoundPlane.h"

class MeasurementEstimator;

class DetRod : public CompositeGSD   {
public:

  virtual ~DetRod();

  /*  It needs that propagator is committed
  virtual vector<DetWithState> 
  compatibleDets( const TrajectoryStateOnSurface& startingState,
		  const Propagator& prop, 
		  const MeasurementEstimator& est) const;
  */

  virtual const BoundSurface& surface() const {return *thePlane;}


  //--- Extension of the interface
  
  /// Return the rod surface as a BoundPlane
  virtual const BoundPlane& specificSurface() const {return *thePlane;}

protected:
  /// Set the rod's plane
  void setPlane( BoundPlane* plane) { thePlane = plane;}

  /// Return the range in Z to be checked for compatibility
  float zError( const TrajectoryStateOnSurface& tsos,
		const MeasurementEstimator& est) const;

 private:
  ReferenceCountingPointer<BoundPlane>  thePlane;

};

#endif
