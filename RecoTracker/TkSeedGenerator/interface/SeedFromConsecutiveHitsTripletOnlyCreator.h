#ifndef RecoTracker_TkSeedGenerator_SeedFromConsecutiveHitsTripletOnlyCreator_H
#define RecoTracker_TkSeedGenerator_SeedFromConsecutiveHitsTripletOnlyCreator_H

#include "RecoTracker/TkSeedGenerator/interface/SeedFromConsecutiveHitsCreator.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
class FreeTrajectoryState;

class SeedFromConsecutiveHitsTripletOnlyCreator : public SeedFromConsecutiveHitsCreator {
public:

  SeedFromConsecutiveHitsTripletOnlyCreator( const edm::ParameterSet & cfg):
    SeedFromConsecutiveHitsCreator(cfg) { }

  virtual ~SeedFromConsecutiveHitsTripletOnlyCreator(){}

protected:

  virtual GlobalTrajectoryParameters initialKinematic(
      const SeedingHitSet::Hits & hits,
      const TrackingRegion & region,
      const edm::EventSetup& es) const;

};
#endif

