#ifndef TkDetLayers_PixelRod_h
#define TkDetLayers_PixelRod_h


#include "TrackingTools/DetLayers/interface/DetRodOneR.h"
#include "TrackingTools/DetLayers/interface/PeriodicBinFinderInZ.h"

/** A concrete implementation for PixelRod
 */

class PixelRod : public DetRodOneR{
 public:
    typedef PeriodicBinFinderInZ<float>   BinFinderType;

  PixelRod(std::vector<const GeomDet*>& theDets);
  ~PixelRod();
  
  // GeometricSearchDet interface

  virtual const std::vector<const GeometricSearchDet*>& components() const;
  
  virtual std::pair<bool, TrajectoryStateOnSurface>
  compatible( const TrajectoryStateOnSurface& ts, const Propagator&, 
	      const MeasurementEstimator&) const;

  virtual void
  compatibleDetsV( const TrajectoryStateOnSurface& startingState,
		  const Propagator& prop, 
		  const MeasurementEstimator& est,
		  std::vector<DetWithState> & result) const;

  virtual void  
  groupedCompatibleDetsV( const TrajectoryStateOnSurface&,
			 const Propagator&,
			 const MeasurementEstimator&,
			 std::vector<DetGroup> &) const;


  virtual bool hasGroups() const {return false;}

 private:
  BinFinderType theBinFinder;
      
  
};


#endif 
