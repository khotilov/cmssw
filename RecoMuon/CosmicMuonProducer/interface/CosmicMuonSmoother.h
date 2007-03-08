#ifndef CosmicMuonProducer_CosmicMuonSmoother_H
#define CosmicMuonProducer_CosmicMuonSmoother_H

/** \file CosmicMuonSmoother
 *
 *  $Date: 2007/03/02 13:19:37 $
 *  $Revision: 1.1 $
 *  \author Chang Liu  -  Purdue University
 */

#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonUtilities.h"

class Propagator;
class KFUpdator;
class MuonServiceProxy;
class Chi2MeasurementEstimator;

namespace edm {class ParameterSet; class Event; class EventSetup;}

class Trajectory;
class TrajectoryMeasurement;

typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
typedef TransientTrackingRecHit::ConstRecHitPointer ConstRecHitPointer;
typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;


class CosmicMuonSmoother : public TrajectorySmoother {
public:


  CosmicMuonSmoother(const edm::ParameterSet&,const MuonServiceProxy* service);
  virtual ~CosmicMuonSmoother();

  virtual std::vector<Trajectory> trajectories(const Trajectory&) const;

  virtual CosmicMuonSmoother* clone() const {
    return new CosmicMuonSmoother(*this);
  }

 /// refit trajectory
    virtual TrajectoryContainer trajectories(const TrajectorySeed& seed,
				             const ConstRecHitContainer& hits, 
				             const TrajectoryStateOnSurface& firstPredTsos) const;


  const Propagator* propagator() const {return &*theService->propagator(thePropagatorName);}

  KFUpdator* updator() const {return theUpdator;}

  CosmicMuonUtilities* utilities() const {return theUtilities; } 

  Chi2MeasurementEstimator* estimator() const {return theEstimator;}

private:

  std::vector<Trajectory> fit(const Trajectory&) const;
  std::vector<Trajectory> fit(const TrajectorySeed& seed,
                              const ConstRecHitContainer& hits,
                              const TrajectoryStateOnSurface& firstPredTsos) const;
  std::vector<Trajectory> smooth(const std::vector<Trajectory>& ) const;
  std::vector<Trajectory> smooth(const Trajectory&) const;

  KFUpdator* theUpdator;
  Chi2MeasurementEstimator* theEstimator;
  CosmicMuonUtilities* theUtilities; 

  const MuonServiceProxy* theService;

  std::string thePropagatorName;
  
};
#endif
