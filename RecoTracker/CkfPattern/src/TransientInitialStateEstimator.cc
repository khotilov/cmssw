#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"

#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "TrackingTools/TrackFitters/interface/KFFittingSmoother.h"



TransientInitialStateEstimator::TransientInitialStateEstimator( const edm::EventSetup& es)
{
  edm::ESHandle<MagneticField>                theMagField;
  es.get<IdealMagneticFieldRecord>().get(theMagField);

  theReversePropagator = new  PropagatorWithMaterial( oppositeToMomentum, 0.105,  
						      &(*theMagField));

  theForwardPropagator = new  PropagatorWithMaterial( alongMomentum, 0.105,  
						      &(*theMagField));
}


std::pair<TrajectoryStateOnSurface, const GeomDet*> 
TransientInitialStateEstimator::innerState( const Trajectory& traj) const
{
  int lastFitted = 4;
  int nhits = traj.foundHits();
  if (nhits < lastFitted+1) lastFitted = nhits-1;

  std::vector<TrajectoryMeasurement> measvec = traj.measurements();
  edm::OwnVector<TransientTrackingRecHit> firstHits;
  for (int i=lastFitted; i >= 0; i--) {
    firstHits.push_back( measvec[i].recHit()->clone());
  }
  TSOS unscaledState = measvec[lastFitted].updatedState();
  AlgebraicSymMatrix C(5,1);
  // C *= 100.;

  TSOS startingState( unscaledState.localParameters(), LocalTrajectoryError(C),
		      unscaledState.surface(),
		      theReversePropagator->magneticField());

  // cout << endl << "FitTester starts with state " << startingState << endl;

  KFTrajectoryFitter backFitter( *theReversePropagator, KFUpdator(), 
				 Chi2MeasurementEstimator( 100., 3));

  vector<Trajectory> fitres = backFitter.fit( traj.seed(), firstHits, startingState);

  if (fitres.size() != 1) {
    // cout << "FitTester: first hits fit failed!" << endl;
    return std::pair<TrajectoryStateOnSurface, const GeomDet*>();
  }

  TrajectoryMeasurement firstMeas = fitres[0].lastMeasurement();
  TSOS firstState = firstMeas.updatedState();

  //  cout << "FitTester: Fitted first state " << firstState << endl;
  //cout << "FitTester: chi2 = " << fitres[0].chiSquared() << endl;

  TSOS initialState( firstState.localParameters(), LocalTrajectoryError(C),
		     firstState.surface(),
		     theReversePropagator->magneticField());

  return std::pair<TrajectoryStateOnSurface, const GeomDet*>( initialState, 
							      firstMeas.recHit()->det());
}

