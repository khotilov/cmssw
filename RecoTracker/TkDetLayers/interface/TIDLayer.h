#ifndef TkDetLayers_TIDLayer_h
#define TkDetLayers_TIDLayer_h


#include "TrackingTools/DetLayers/interface/RingedForwardLayer.h"
#include "RecoTracker/TkDetLayers/interface/TIDRing.h"


/** A concrete implementation for TID layer 
 *  built out of TIDRings
 */

class TIDLayer : public RingedForwardLayer{
 public:
  TIDLayer(std::vector<const TIDRing*>& rings);
  ~TIDLayer();
  
  // GeometricSearchDet interface
  
  virtual const std::vector<const GeomDet*>& basicComponents() const {return theBasicComps;}
  
  virtual const std::vector<const GeometricSearchDet*>& components() const {return theComps;}

  virtual std::vector<DetWithState> 
    compatibleDets( const TrajectoryStateOnSurface& startingState,
		    const Propagator& prop, 
		    const MeasurementEstimator& est) const;
  
  virtual std::vector<DetGroup> 
  groupedCompatibleDets( const TrajectoryStateOnSurface& startingState,
			 const Propagator& prop,
			 const MeasurementEstimator& est) const;


  virtual bool hasGroups() const {return true;}

  // DetLayer interface
  virtual Module   module()   const { return silicon;}


 private:
  // private methods for the implementation of groupedCompatibleDets()
  virtual BoundDisk* computeDisk( const std::vector<const TIDRing*>& rings) const;

  virtual std::vector<int> ringIndicesByCrossingProximity(const TrajectoryStateOnSurface& startingState,
						     const Propagator& prop ) const;

 protected:  
  //  bool isCompatible( const TrajectoryStateOnSurface& ms,
  //	     const MeasurementEstimator& est) const;

  int findClosest( const std::vector<GlobalPoint>& ) const;
  
  int findNextIndex( const std::vector<GlobalPoint>& , int ) const;
  
  bool overlapInR( const TrajectoryStateOnSurface& tsos, int i, double ymax) const;
  
  
  float computeWindowSize( const GeomDet* det, 
  			   const TrajectoryStateOnSurface& tsos, 
			   const MeasurementEstimator& est) const;
  
  std::vector<DetGroup> orderAndMergeLevels(const TrajectoryStateOnSurface& tsos,
					    const Propagator& prop,
					    const std::vector<std::vector<DetGroup> > groups,
					    const std::vector<int> indices ) const;




 protected:
  std::vector<const GeometricSearchDet*> theComps;
  std::vector<const GeomDet*> theBasicComps;
  
};


#endif 
