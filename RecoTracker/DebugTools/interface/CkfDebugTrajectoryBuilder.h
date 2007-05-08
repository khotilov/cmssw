#ifndef CkfDebugTrajectoryBuilder_H
#define CkfDebugTrajectoryBuilder_H

#include "RecoTracker/CkfPattern/interface/CkfTrajectoryBuilder.h"
#include "RecoTracker/DebugTools/interface/CkfDebugger.h"


class CkfDebugTrajectoryBuilder: public CkfTrajectoryBuilder{
 public:

  CkfDebugTrajectoryBuilder(const edm::ParameterSet&              conf,
                       const TrajectoryStateUpdator*         updator,
                       const Propagator*                     propagatorAlong,
                       const Propagator*                     propagatorOpposite,
                       const Chi2MeasurementEstimatorBase*   estimator,
                       const TransientTrackingRecHitBuilder* RecHitBuilder,
                       const MeasurementTracker*             measurementTracker) : 
    CkfTrajectoryBuilder( conf,updator,propagatorAlong,propagatorOpposite,estimator,RecHitBuilder,measurementTracker) {}

  virtual void setDebugger( CkfDebugger * dbg) { theDbg = dbg;}
  virtual CkfDebugger * debugger() const{ return theDbg;}

 private:
  mutable CkfDebugger * theDbg;
  bool analyzeMeasurementsDebugger(Trajectory& traj, std::vector<TM> meas,
				   const MeasurementTracker* theMeasurementTracker, const Propagator* theForwardPropagator, 
				   const Chi2MeasurementEstimatorBase* theEstimator, 
				   const TransientTrackingRecHitBuilder * theTTRHBuilder) const { 
    return theDbg->analyseCompatibleMeasurements(traj,meas,theMeasurementTracker,theForwardPropagator,theEstimator,theTTRHBuilder);
  };
  void fillSeedHistoDebugger(std::vector<TrajectoryMeasurement>::iterator begin, 
			     std::vector<TrajectoryMeasurement>::iterator end) const {
    if (end-begin>=2)
      theDbg->fillSeedHist(begin->recHit(),(begin+1)->recHit(),(begin+1)->updatedState());
  }; 

};
#endif
