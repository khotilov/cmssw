#ifndef TkDetLayers_PixelForwardLayer_h
#define TkDetLayers_PixelForwardLayer_h


#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "RecoTracker/TkDetLayers/interface/PixelBlade.h"
#include "Utilities/BinningTools/interface/PeriodicBinFinderInPhi.h"


/** A concrete implementation for PixelForward layer 
 *  built out of ForwardPixelBlade
 */

class PixelForwardLayer : public ForwardDetLayer, public GeometricSearchDetWithGroups {
 public:
  PixelForwardLayer(std::vector<const PixelBlade*>& blades);
  ~PixelForwardLayer();
  
  // GeometricSearchDet interface
  
  virtual const std::vector<const GeomDet*>& basicComponents() const {return theBasicComps;}

  virtual const std::vector<const GeometricSearchDet*>& components() const {return theComps;}
  
  void groupedCompatibleDetsV( const TrajectoryStateOnSurface& tsos,
			       const Propagator& prop,
			       const MeasurementEstimator& est,
			       std::vector<DetGroup> & result) const;

  // DetLayer interface
  virtual SubDetector subDetector() const {return GeomDetEnumerators::PixelEndcap;}
  

 private:  
  // methods for groupedCompatibleDets implementation
  int computeHelicity(const GeometricSearchDet* firstBlade,const GeometricSearchDet* secondBlade) const;

  struct SubTurbineCrossings {
    SubTurbineCrossings(){};
    SubTurbineCrossings( int ci, int ni, float nd) : 
      closestIndex(ci), nextIndex(ni), nextDistance(nd) {}
    
    int   closestIndex;
    int   nextIndex;
    float nextDistance;
  };
  
  void searchNeighbors( const TrajectoryStateOnSurface& tsos,
			const Propagator& prop,
			const MeasurementEstimator& est,
			const SubTurbineCrossings& crossings,
			float window, 
			std::vector<DetGroup>& result) const;
  
  SubTurbineCrossings 
    computeCrossings( const TrajectoryStateOnSurface& startingState,
		      PropagationDirection propDir) const;

  float computeWindowSize( const GeomDet* det, 
			   const TrajectoryStateOnSurface& tsos, 
			   const MeasurementEstimator& est) const;
  
 private:
  typedef PeriodicBinFinderInPhi<double>   BinFinderType;
  BinFinderType    theBinFinder;

  std::vector<const GeometricSearchDet*> theComps;
  std::vector<const GeomDet*> theBasicComps;
};


#endif 
