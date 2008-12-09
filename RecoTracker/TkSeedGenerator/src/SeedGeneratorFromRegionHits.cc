#include "RecoTracker/TkSeedGenerator/interface/SeedGeneratorFromRegionHits.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedComparitor.h"
#include "RecoTracker/TkSeedGenerator/interface/SeedCreator.h"
#include "RecoTracker/TkSeedGenerator/interface/SeedFromConsecutiveHitsCreator.h"
#include "RecoTracker/TkSeedGenerator/interface/SeedFromConsecutiveHits.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

template <class T> T sqr( T t) {return t*t;}

SeedGeneratorFromRegionHits::SeedGeneratorFromRegionHits(
    OrderedHitsGenerator *ohg, SeedComparitor* asc, SeedCreator* asp)
  : theHitsGenerator(ohg), theComparitor(asc), theSeedCreator(asp)
{ }

SeedGeneratorFromRegionHits::SeedGeneratorFromRegionHits(
    OrderedHitsGenerator * aGenerator,
    const edm::ParameterSet & cfg,
    SeedComparitor * aComparitor)
  : theHitsGenerator(aGenerator), theComparitor(aComparitor)
{
  theSeedCreator = new SeedFromConsecutiveHitsCreator();
}


SeedGeneratorFromRegionHits::~SeedGeneratorFromRegionHits()
{
  delete theHitsGenerator;
  delete theComparitor;
  delete theSeedCreator;
}

void SeedGeneratorFromRegionHits::run(TrajectorySeedCollection & seedCollection, 
    const TrackingRegion & region, const edm::Event& ev, const edm::EventSetup& es)
{
  const OrderedSeedingHits & hitss = theHitsGenerator->run(region, ev, es);

  unsigned int nHitss =  hitss.size();
  for (unsigned int iHits = 0; iHits < nHitss; ++iHits) { 
    const SeedingHitSet & hits =  hitss[iHits];
    if (!theComparitor || theComparitor->compatible( hits, es) ) {
      theSeedCreator->trajectorySeed(seedCollection, hits, region, es);
    }
  }
}
