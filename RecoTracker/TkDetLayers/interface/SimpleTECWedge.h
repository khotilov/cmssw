#ifndef TkDetLayers_SimpleTECWedge_h
#define TkDetLayers_SimpleTECWedge_h


#include "RecoTracker/TkDetLayers/interface/TECWedge.h"


/** A concrete implementation for TEC wedge
 *  built out of only one det
 */

class SimpleTECWedge : public TECWedge{
 public:
  SimpleTECWedge(const GeomDet* theDet);

  ~SimpleTECWedge();
  
  // GeometricSearchDet interface
  virtual const vector<const GeomDet*>& basicComponents() const {return theDets;}

  virtual const vector<const GeometricSearchDet*>& components() const;
  
  virtual pair<bool, TrajectoryStateOnSurface>
  compatible( const TrajectoryStateOnSurface& ts, const Propagator&, 
	      const MeasurementEstimator&) const;

  virtual vector<DetGroup> 
  groupedCompatibleDets( const TrajectoryStateOnSurface& startingState,
			 const Propagator& prop,
			 const MeasurementEstimator& est) const;

 private:
  const GeomDet* theDet;
  vector<const GeomDet*> theDets;

};


#endif 
