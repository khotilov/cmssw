#ifndef PixelTrackCleaner_H
#define PixelTrackCleaner_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include <utility>
#include <vector>

/*
PixelTrackCleaner:

Discards tracks with more than one common recHit.

*/


class PixelTrackCleaner {

public:

  typedef std::pair<reco::Track*, std::vector<const TrackingRecHit *> > TrackWithRecHits;

	PixelTrackCleaner();

	std::vector<TrackWithRecHits> cleanTracks(std::vector<TrackWithRecHits> trackWithRecHits);
  void cleanTrack();
  bool recHitsAreEqual(const TrackingRecHit *recHit1, const TrackingRecHit *recHit2);

private:

  std::vector<bool> trackOk;
  reco::Track *track1, *track2;
  int iTrack1, iTrack2;

};

#endif
