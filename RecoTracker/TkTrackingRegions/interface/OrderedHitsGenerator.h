#ifndef TkTrackingRegions_OrderedHitsGenerator_H
#define TkTrackingRegions_OrderedHitsGenerator_H

#include "RecoTracker/TkSeedingLayers/interface/OrderedSeedingHits.h"
#include <vector>

class TrackingRegion;
namespace edm { class Event; class EventSetup; };

class OrderedHitsGenerator {
public:
  virtual ~OrderedHitsGenerator() {}

  virtual const OrderedSeedingHits & run( 
      const TrackingRegion& reg, const edm::Event & ev, const edm::EventSetup& es ) = 0;

//  virtual OrderedHitsGenerator * clone() const = 0;

};
#endif
