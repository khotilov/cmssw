#ifndef RecoMuon_L3TrackFinder_MuonCkfTrajectoryBuilder_H
#define RecoMuon_L3TrackFinder_MuonCkfTrajectoryBuilder_H

#include "RecoTracker/CkfPattern/interface/CkfTrajectoryBuilder.h"

class MuonCkfTrajectoryBuilder : public CkfTrajectoryBuilder {
 public:
  MuonCkfTrajectoryBuilder(const edm::ParameterSet&              conf,
			   const TrajectoryStateUpdator*         updator,
			   const Propagator*                     propagatorAlong,
			   const Propagator*                     propagatorOpposite,
			   const Propagator*                     propagatorProximity,
			   const Chi2MeasurementEstimatorBase*   estimator,
			   const TransientTrackingRecHitBuilder* RecHitBuilder,
			   const MeasurementTracker*             measurementTracker);
  virtual ~MuonCkfTrajectoryBuilder();
  
 protected:
  void collectMeasurement(const std::vector<const DetLayer*>& nl,const TrajectoryStateOnSurface & currentState, std::vector<TM>& result,int& invalidHits,const Propagator *) const;

  virtual void findCompatibleMeasurements( const TempTrajectory& traj, std::vector<TrajectoryMeasurement> & result) const;
  
  //and other fields
  bool theUseSeedLayer;
  double theRescaleErrorIfFail;
  mutable const Propagator * theProximityPropagator;
};


#endif
