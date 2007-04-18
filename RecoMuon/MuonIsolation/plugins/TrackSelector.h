#ifndef MuonIsolation_TrackSelector_H
#define MuonIsolation_TrackSelector_H

#include "RecoMuon/MuonIsolation/interface/Range.h"
#include "DataFormats/MuonReco/interface/Direction.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


namespace muonisolation {

class TrackSelector {
public:

  typedef muonisolation::Range<float> Range;

  TrackSelector(const Range & z, float r, const Direction &dir, float drMax);
  reco::TrackCollection operator()(const reco::TrackCollection & tracks) const;


private:
  Range theZ;
  Range theR;
  Direction theDir;
  float theDR_Max;

}; 

}

#endif 
