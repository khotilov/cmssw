#ifndef TkDetLayers_TIDRing_h
#define TkDetLayers_TIDRing_h


#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "Utilities/BinningTools/interface/PeriodicBinFinderInPhi.h"
#include "RecoTracker/TkDetLayers/interface/SubLayerCrossings.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"

/** A concrete implementation for TID rings 
 */

class TIDRing : public GeometricSearchDetWithGroups{
 public:
  TIDRing(std::vector<const GeomDet*>& innerDets,
	  std::vector<const GeomDet*>& outerDets);
  ~TIDRing();
  
  // GeometricSearchDet interface
  virtual const BoundSurface& surface() const {return *theDisk;}
  
  virtual const std::vector<const GeomDet*>& basicComponents() const {return theDets;}
  
  virtual const std::vector<const GeometricSearchDet*>& components() const;

  void groupedCompatibleDetsV( const TrajectoryStateOnSurface& tsos,
			       const Propagator& prop,
			       const MeasurementEstimator& est,
			       std::vector<DetGroup> & result) const;
  
 
  //Extension of interface
  virtual const BoundDisk& specificSurface() const {return *theDisk;}
  

 private:
  // private methods for the implementation of groupedCompatibleDets()

  SubLayerCrossings computeCrossings( const TrajectoryStateOnSurface& tsos,
				      PropagationDirection propDir) const;
  
  bool addClosest( const TrajectoryStateOnSurface& tsos,
		   const Propagator& prop,
		   const MeasurementEstimator& est,
		   const SubLayerCrossing& crossing,
		   std::vector<DetGroup>& result) const;

  float computeWindowSize( const GeomDet* det, 
			   const TrajectoryStateOnSurface& tsos, 
			   const MeasurementEstimator& est) const;

  float calculatePhiWindow( const MeasurementEstimator::Local2DVector&  maxDistance, 
  			    const TrajectoryStateOnSurface& ts, 
			    const BoundPlane& plane) const;

  std::pair<float, float> computeDetPhiRange( const BoundPlane& plane) const;
  

  void searchNeighbors( const TrajectoryStateOnSurface& tsos,
			const Propagator& prop,
			const MeasurementEstimator& est,
			const SubLayerCrossing& crossing,
			float window, 
			std::vector<DetGroup>& result,
			bool checkClosest) const;

  bool overlapInPhi( const GlobalPoint& startPoint,const GeomDet* det, float phiWin ) const;  

  const std::vector<const GeomDet*>& subLayer( int ind) const {
    return (ind==0 ? theFrontDets : theBackDets);
  }


 private:
  std::vector<const GeomDet*> theDets;
  std::vector<const GeomDet*> theFrontDets;
  std::vector<const GeomDet*> theBackDets;

  ReferenceCountingPointer<BoundDisk> theDisk;
  ReferenceCountingPointer<BoundDisk> theFrontDisk;
  ReferenceCountingPointer<BoundDisk> theBackDisk;

  typedef PeriodicBinFinderInPhi<double>   BinFinderType;

  BinFinderType theFrontBinFinder;
  BinFinderType theBackBinFinder;


  
  };


#endif 
