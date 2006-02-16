#ifndef ForwardDetRing_H
#define ForwardDetRing_H

/** \class ForwardDetRing
 *  Abstract interface for a ring of detectors sitting on a BoundDisk.
 */

#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "Geometry/Surface/interface/BoundDisk.h"

class ForwardDetRing : public GeometricSearchDet {
 public:

  virtual ~ForwardDetRing();

  
  virtual vector<DetWithState> 
  compatibleDets( const TrajectoryStateOnSurface& tsos,
		  const Propagator& prop, 
		  const MeasurementEstimator& est) const;
  
  virtual const BoundSurface& surface() const {return *theDisk;}

  
  //--- Extension of the interface

  /// Return the ring surface as a BoundDisk
  const BoundDisk& specificSurface() const {return *theDisk;}


protected:

  /// Set the rod's disk
  void setDisk( BoundDisk* disk) { theDisk = disk;}

  
 private:
  ReferenceCountingPointer<BoundDisk>  theDisk;

};
#endif

