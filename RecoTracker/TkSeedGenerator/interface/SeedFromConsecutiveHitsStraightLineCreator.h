#ifndef RecoTracker_TkSeedGenerator_SeedFromConsecutiveHitsStraightLineCreator_H
#define RecoTracker_TkSeedGenerator_SeedFromConsecutiveHitsStraightLineCreator_H

#include "RecoTracker/TkSeedGenerator/interface/SeedFromConsecutiveHitsCreator.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"
class FreeTrajectoryState;

class SeedFromConsecutiveHitsStraightLineCreator : public SeedFromConsecutiveHitsCreator {
public:

  SeedFromConsecutiveHitsStraightLineCreator( const edm::ParameterSet & cfg):
    SeedFromConsecutiveHitsCreator(cfg) { }

  virtual ~SeedFromConsecutiveHitsStraightLineCreator(){}

protected:

  virtual GlobalTrajectoryParameters initialKinematic(
      const SeedingHitSet & hits,
      const TrackingRegion & region,
      const edm::EventSetup& es) const;

};
#endif

