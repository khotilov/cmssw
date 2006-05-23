#ifndef TkNavigation_StartingLayerFinder_H_
#define TkNavigation_StartingLayerFinder_H_


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"


#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

#include <vector>

/** Finds the nearest navigable layer.
 *  Needed to start trajectory building in case the seed does not 
 *  have a DetLayer
 */

class PTrajectoryStateOnDet;

class StartingLayerFinder {

private:

  typedef FreeTrajectoryState FTS;
  typedef TrajectoryStateOnSurface TSOS;
  typedef pair<float, float> Range;

public: 

  StartingLayerFinder(const Propagator* aPropagator, const MeasurementTracker*  tracker ) :

    thePropagator(aPropagator),
    theMeasurementTracker(tracker),
    thePixelLayersValid(false),
    theFirstPixelBarrelLayer(0),
    theFirstNegPixelFwdLayer(0),
    theFirstPosPixelFwdLayer(0) { }

  ~StartingLayerFinder() {}

  vector<const DetLayer*> startingLayers(const FTS& aFts, float dr, float dz) const;
  

  vector<const DetLayer*> startingLayers(const TrajectorySeed& aSeed) const;

  const BarrelDetLayer* firstPixelBarrelLayer() const;
  const vector<ForwardDetLayer*> firstNegPixelFwdLayer() const;
  const vector<ForwardDetLayer*> firstPosPixelFwdLayer() const;

  const Propagator* propagator() const {return thePropagator;}



private:

  const Propagator* thePropagator;
  const MeasurementTracker*     theMeasurementTracker;
  mutable bool thePixelLayersValid;
  mutable BarrelDetLayer* theFirstPixelBarrelLayer;
  mutable vector<ForwardDetLayer*> theFirstNegPixelFwdLayer;
  mutable vector<ForwardDetLayer*> theFirstPosPixelFwdLayer;


  void checkPixelLayers() const;
  
  
  
  inline bool rangesIntersect( const Range& a, const Range& b) const {
    if ( a.first > b.second || b.first > a.second) return false;
    else return true;
  }




};
#endif //TR_StartingLayerFinder_H_
