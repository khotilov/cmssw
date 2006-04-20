#ifndef GeometricSearchDetMeasurements_H
#define GeometricSearchDetMeasurements_H

#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDetException.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDetSystem.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDet.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"
#include "TrackingTools/TransientTrackingRecHit/interface/InvalidTransientRecHit.h"


#include <algorithm>

/** Concrete implementation of the Det::measurements method
 *  for GeometricSearchDets. It is done in a separate class to facilitate
 *  it's reuse from the various GeometricSearchDets.
 */

class GeometricSearchDetMeasurements {
public:

  typedef GeometricSearchDet::DetWithState DetWithState;
  typedef TrajectoryStateOnSurface   TSOS;
  typedef TrajectoryMeasurement      TM;
  typedef TrajectoryMeasurementGroup TMG;

  GeometricSearchDetMeasurements( const MeasurementDetSystem* detSysytem) :
    theDetSystem(detSysytem) {}



  /*
  template <class TrajectoryState> 
  std::vector<DetWithState> getDets( const GeometricSearchDet& det,
				const TrajectoryState& ts, 
				const Propagator& prop, 
				const MeasurementEstimator& est) const {
    pair<bool, TSOS> compat = det.compatible( ts, prop, est);
    if ( compat.first) return det.fastCompatibleDets( compat.second, ts, prop, est);
    else return std::vector<DetWithState>();
  }

  template <class TrajectoryState> 
  std::vector<DetWithState> getDets( const GeometricSearchDet& det,
				const TrajectoryStateOnSurface& stateOnDet,
				const TrajectoryState& ts, 
				const Propagator& prop, 
				const MeasurementEstimator& est) const {
    return det.fastCompatibleDets( stateOnDet, ts, prop, est);
  }

  template <class TrajectoryState>
  std::vector<TrajectoryMeasurement> get( const GeometricSearchDet& det,
				     const TrajectoryState& ts, 
				     const Propagator& prop, 
				     const MeasurementEstimator& est) const {
    pair<bool, TSOS> compat = det.compatible( ts, prop, est);
    if ( compat.first) {
      return det.fastMeasurements( compat.second, ts, prop, est);
    }
    else return std::vector<TrajectoryMeasurement>();
  }
    
  template <class TrajectoryState>
  std::vector<TrajectoryMeasurement> get( const GeometricSearchDet& det,
				     const TrajectoryStateOnSurface& stateOnDet,
				     const TrajectoryState& ts, 
				     const Propagator& prop, 
				     const MeasurementEstimator& est) const {
    std::vector<DetWithState> compatDets = det.fastCompatibleDets( stateOnDet, ts, prop, est);
    if (!compatDets.empty()) {
      return get( det, compatDets, ts, prop, est);
    }
    else {
      std::vector<TrajectoryMeasurement> result;
      addInvalidMeas( result, stateOnDet, &det);
      return result;
    }
  }
  */

  /** The std::vector<DetWithState> passed to this method should not be empty.
   *  In case of no compatible dets the result should be either an empty container if 
   *  the det is itself incompatible, or a container with one invalid measurement
   *  on the det surface. The method does not have enough information to do
   *  this efficiently, so it should be done by the caller, or an exception will
   *  be thrown (DetLogicError).
   */
  template <class TrajectoryState>
  std::vector<TrajectoryMeasurement> get( const GeometricSearchDet& det,
					  const std::vector<DetWithState>& compatDets,
					  const TrajectoryState& ts, 
					  const Propagator& prop, 
					  const MeasurementEstimator& est) const;
  /*
  template <class TrajectoryState>
  std::vector<TMG>
  getGrouped( const GeometricSearchDet& det, const std::vector<DetWithState>& dg,
	      const TrajectoryState& ts, const Propagator& prop,
	      const MeasurementEstimator& est) const {
    if (!dg.empty()) {
      std::vector<TMG> result(1);
      result[0] = TMG( get( det, dg, ts, prop, est), dg);
      return result;
    }
    else return std::vector<TMG>();
  }
  */

  
  void addInvalidMeas( std::vector<TrajectoryMeasurement>& result, 
		       const TrajectoryStateOnSurface& ts, const GeomDet* det) const {
    result.push_back( TM( ts, new InvalidTransientRecHit(det), 0.F,0));
  }
  

private:

  const MeasurementDetSystem* theDetSystem;

};


template <class TrajectoryState>
std::vector<TrajectoryMeasurement> 
GeometricSearchDetMeasurements::get( const GeometricSearchDet& det,
				     const std::vector<DetWithState>& compatDets,
				     const TrajectoryState& ts, 
				     const Propagator& prop, 
				     const MeasurementEstimator& est) const
{
  std::vector<TrajectoryMeasurement> result;
  if (!compatDets.empty()) {
    for ( std::vector<DetWithState>::const_iterator i=compatDets.begin();
	  i != compatDets.end(); i++) {
      const MeasurementDet* mdet = theDetSystem->idToDet(i->first->geographicalId());
      if (mdet == 0) {
	throw MeasurementDetException( "MeasurementDet not found");
      }

      std::vector<TM> tmp = mdet->fastMeasurements( i->second, ts, prop, est);
      if ( !tmp.empty()) {
	// only collect valid RecHits
	std::vector<TM>::iterator end = (tmp.back().recHit()->isValid() ? tmp.end() : tmp.end()-1);
	result.insert( result.end(), tmp.begin(), end);
      }
    }
    // sort the final result
    if ( result.size() > 1) {
      sort( result.begin(), result.end(), TrajMeasLessEstim());
    }

    
    if ( !result.empty()) {
      // invalidMeas on Det of most compatible hit
      addInvalidMeas( result, result.front().predictedState(), result.front().recHit()->det());
    }
    else {
      cout << "no compatible meas found in mdet->fastMeasurements" << endl;
      // invalid state on first compatible Det
      addInvalidMeas( result, compatDets.front().second, compatDets.front().first);
    }
    

  }
  else {
    // this case should have been handled by the caller!
    throw MeasurementDetException("GeometricSearchDetMeasurements::get called with empty std::vector<DetWithState>");
  }
  return result;
}


#endif
