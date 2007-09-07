#include "TrackingTools/GsfTracking/interface/GsfTrajectorySmoother.h"

#include "TrackingTools/GsfTracking/interface/MultiTrajectoryStateMerger.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

GsfTrajectorySmoother::GsfTrajectorySmoother(const GsfPropagatorWithMaterial& aPropagator,
					     const TrajectoryStateUpdator& aUpdator,
					     const MeasurementEstimator& aEstimator,
					     const MultiTrajectoryStateMerger& aMerger,
					     float errorRescaling,
					     const bool materialBeforeUpdate) :
  thePropagator(aPropagator.clone()),
  theGeomPropagator(0),
  theConvolutor(0),
  theUpdator(aUpdator.clone()),
  theEstimator(aEstimator.clone()),
  theMerger(aMerger.clone()),
  theMatBeforeUpdate(materialBeforeUpdate),
  theErrorRescaling(errorRescaling)
{
  if ( !theMatBeforeUpdate ) {
    theGeomPropagator = new GsfPropagatorAdapter(thePropagator->geometricalPropagator());
    theConvolutor = thePropagator->convolutionWithMaterial().clone();
  }

  //   static SimpleConfigurable<bool> timeConf(false,"GsfTrajectorySmoother:activateTiming");
  //   theTiming = timeConf.value();
}

GsfTrajectorySmoother::~GsfTrajectorySmoother() {
  delete thePropagator;
  delete theGeomPropagator;
  delete theConvolutor;
  delete theUpdator;
  delete theEstimator;
  delete theMerger;
}

std::vector<Trajectory> 
GsfTrajectorySmoother::trajectories(const Trajectory& aTraj) const {

  //   static TimingReport::Item* propTimer =
  //     &(*TimingReport::current())[string("GsfTrajectorySmoother:propagation")];
  //   propTimer->switchCPU(false);
  //   if ( !theTiming )  propTimer->switchOn(false);
  //   static TimingReport::Item* updateTimer =
  //     &(*TimingReport::current())[string("GsfTrajectorySmoother:update")];
  //   updateTimer->switchCPU(false);
  //   if ( !theTiming )  updateTimer->switchOn(false);
  
  if(aTraj.empty()) return std::vector<Trajectory>();
  
  if (  aTraj.direction() == alongMomentum) {
    thePropagator->setPropagationDirection(oppositeToMomentum);
  }
  else {
    thePropagator->setPropagationDirection(alongMomentum);
  }

  Trajectory myTraj(aTraj.seed(), propagator()->propagationDirection());
  
  std::vector<TM> avtm = aTraj.measurements();
  
  TSOS predTsos = avtm.back().forwardPredictedState();
  predTsos.rescaleError(theErrorRescaling);

  if(!predTsos.isValid()) {
    edm::LogInfo("GsfTrajectorySmoother") 
      << "GsfTrajectorySmoother: predicted tsos of last measurement not valid!";
    return std::vector<Trajectory>();
  }
  TSOS currTsos;

  //first smoothed tm is last fitted
  if(avtm.back().recHit()->isValid()) {
    {
      //       TimeMe t(*updateTimer,false);
      currTsos = updator()->update(predTsos, *avtm.back().recHit());
    }
    if(!currTsos.isValid()) {
      edm::LogInfo("GsfTrajectorySmoother") << "GsfTrajectorySmoother: tsos not valid after update!";
      return std::vector<Trajectory>();
    }
    myTraj.push(TM(avtm.back().forwardPredictedState(), 
		   predTsos,
		   avtm.back().updatedState(),
		   avtm.back().recHit(),
		   avtm.back().estimate()), 
		avtm.back().estimate());
  } else {
    currTsos = predTsos;
    myTraj.push(TM(avtm.back().forwardPredictedState(),
		   avtm.back().recHit()));
  }
  
  TrajectoryStateCombiner combiner;

  for(std::vector<TM>::reverse_iterator itm = avtm.rbegin() + 1; 
      itm < avtm.rend() - 1; ++itm) {
    {
      //       TimeMe t(*propTimer,false);
      //       //
      //       // check type of surface in case of invalid hit
      //       // (in this version only propagations to planes are
      //       // supported for multi trajectory states)
      //       //
      //       if ( !(*itm).recHit().isValid() ) {
      // 	const BoundPlane* plane = 
      // 	  dynamic_cast<const BoundPlane*>(&(*itm).recHit().det().surface());
      // 	//
      // 	// no plane: insert invalid 
      // 	if ( plane==0 ) {
      // 	  myTraj.push(TM(TrajectoryStateOnSurface(),
      // 			 (*itm).recHit());
      // 	  continue;
      // 	}
      //     }
      predTsos = propagator()->propagate(currTsos,
					 (*itm).recHit()->det()->surface());
    }
    if ( predTsos.isValid() && theConvolutor && theMatBeforeUpdate )
      predTsos = (*theConvolutor)(predTsos,
				  propagator()->propagationDirection());
    if(!predTsos.isValid()) {
      edm::LogInfo("GsfTrajectorySmoother") << "GsfTrajectorySmoother: predicted tsos not valid!";
      return std::vector<Trajectory>();
    }
    if ( theMerger )  predTsos = theMerger->merge(predTsos);

    if((*itm).recHit()->isValid()) {
      //update
      {
	// 	TimeMe t(*updateTimer,false);
	currTsos = updator()->update(predTsos, *(*itm).recHit());
      }
      if ( currTsos.isValid() && theConvolutor && !theMatBeforeUpdate )
	currTsos = (*theConvolutor)(currTsos,
				    propagator()->propagationDirection());
      if(!currTsos.isValid()) {
	edm::LogInfo("GsfTrajectorySmoother") 
	  << "GsfTrajectorySmoother: tsos not valid after update / material effects!";
	return std::vector<Trajectory>();
      }
      //3 different possibilities to calculate smoothed state:
      //1: update combined predictions with hit
      //2: combine fwd-prediction with bwd-filter
      //3: combine bwd-prediction with fwd-filter
      TSOS combTsos = combiner(predTsos, (*itm).forwardPredictedState());
      if(!combTsos.isValid()) {
	LogDebug("GsfTrajectorySmoother") << 
	  "KFTrajectorySmoother: combined tsos not valid!\n"<<
	  "pred Tsos pos: "<<predTsos.globalPosition()<< "\n" <<
	  "pred Tsos mom: "<<predTsos.globalMomentum()<< "\n" <<
	  "TrackingRecHit: "<<(*itm).recHit()->surface()->toGlobal((*itm).recHit()->localPosition())<< "\n" ;
	return std::vector<Trajectory>();
      }

      TSOS smooTsos = combiner((*itm).updatedState(), predTsos);

      if(!smooTsos.isValid()) {
	LogDebug("GsfTrajectorySmoother") <<
	  "KFTrajectorySmoother: smoothed tsos not valid!";
	return std::vector<Trajectory>();
      }

      myTraj.push(TM((*itm).forwardPredictedState(),
		     predTsos,
		     smooTsos,
		     (*itm).recHit(),
		     estimator()->estimate(combTsos, *(*itm).recHit()).second),
		  (*itm).estimate());
    } 
    else {
      currTsos = predTsos;
      TSOS combTsos = combiner(predTsos, (*itm).forwardPredictedState());
      
      if(!combTsos.isValid()) {
    	LogDebug("GsfTrajectorySmoother") << 
    	  "KFTrajectorySmoother: combined tsos not valid!";
    	return std::vector<Trajectory>();
      }

      myTraj.push(TM((*itm).forwardPredictedState(),
    		     predTsos,
    		     combTsos,
    		     (*itm).recHit()));
    }
    if ( theMerger )  currTsos = theMerger->merge(currTsos);
  }

  //last smoothed tm is last filtered
  {
    //     TimeMe t(*propTimer,false);
    predTsos = propagator()->propagate(currTsos,
				       avtm.front().recHit()->det()->surface());
  }
  if ( predTsos.isValid() && theConvolutor && theMatBeforeUpdate )
    predTsos = (*theConvolutor)(predTsos,
				propagator()->propagationDirection());
  if(!predTsos.isValid()) {
    edm::LogInfo("GsfTrajectorySmoother") << "GsfTrajectorySmoother: predicted tsos not valid!";
    return std::vector<Trajectory>();
  }
  if ( theMerger )  predTsos = theMerger->merge(predTsos);

  if(avtm.front().recHit()->isValid()) {
    //update
    {
      //       TimeMe t(*updateTimer,false);
      currTsos = updator()->update(predTsos, *avtm.front().recHit());
    }
    if ( currTsos.isValid() && theConvolutor && !theMatBeforeUpdate )
      currTsos = (*theConvolutor)(currTsos,
				  propagator()->propagationDirection());
    if(!currTsos.isValid()) {
      edm::LogInfo("GsfTrajectorySmoother") 
	<< "GsfTrajectorySmoother: tsos not valid after update / material effects!";
      return std::vector<Trajectory>();
    }
  
    myTraj.push(TM(avtm.front().forwardPredictedState(),
		   predTsos,
		   currTsos,
		   avtm.front().recHit(),
		   estimator()->estimate(predTsos, *avtm.front().recHit()).second),
		avtm.front().estimate());
    //estimator()->estimate(predTsos, avtm.front().recHit()));
  } 
  else {
    myTraj.push(TM(avtm.front().forwardPredictedState(),
		   avtm.front().recHit()));
  }

  return std::vector<Trajectory>(1, myTraj); 
}
