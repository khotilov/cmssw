//
// $Id: TrackerIsolationPt.cc,v 1.1 2008/01/07 11:48:27 lowette Exp $
//

#include "PhysicsTools/PatUtils/interface/TrackerIsolationPt.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include <vector>


using namespace pat;


/// constructor
TrackerIsolationPt::TrackerIsolationPt() {
}


/// destructor
TrackerIsolationPt::~TrackerIsolationPt() {
}


/// calculate the TrackIsoPt for the lepton object
float TrackerIsolationPt::calculate(const Electron & theElectron, const reco::TrackCollection & theTracks, float isoConeElectron) const {
  return this->calculate(*theElectron.gsfTrack(), theTracks, isoConeElectron);
}
float TrackerIsolationPt::calculate(const Muon & theMuon, const reco::TrackCollection & theTracks, float isoConeMuon) const {
  return this->calculate(*theMuon.track(), theTracks, isoConeMuon);
}


/// calculate the TrackIsoPt for the lepton's track
float TrackerIsolationPt::calculate(const reco::Track & theTrack, const reco::TrackCollection & theTracks, float isoCone) const {
  // initialize some variables
  float isoPtLepton = 0;
  const reco::Track * closestTrackDRPT = 0, * closestTrackDR = 0;
  float closestDRPT = 10000, closestDR = 10000;
  // use all these pointless vector conversions because the momenta from tracks
  // are completely unusable; bah, these math-vectors are worthless!
  HepLorentzVector lepton(theTrack.px(), theTrack.py(), theTrack.pz(), theTrack.p());
  for (std::vector<reco::Track>::const_iterator itTrack = theTracks.begin(); itTrack != theTracks.end(); itTrack++) {
    HepLorentzVector track(itTrack->px(), itTrack->py(), itTrack->pz(), itTrack->p());
    float dR = lepton.deltaR(track);
    if (dR < isoCone) {
      isoPtLepton += track.perp();
      // find the closest matching track
      // FIXME: we could association by hits or chi2 to match
      float pRatio = track.perp()/lepton.perp();
      if (dR < closestDRPT && pRatio > 0.5 && pRatio < 1.5) {
        closestDRPT = dR;
        closestTrackDRPT = &*itTrack;
      }
      if (dR < closestDR) {
        closestDR = dR;
        closestTrackDR = &*itTrack;
      }
    }
  }
  if (closestTrackDRPT) {
    GlobalVector closestTrackVector(closestTrackDRPT->px(), closestTrackDRPT->py(), closestTrackDRPT->pz());
    isoPtLepton -= closestTrackVector.perp();
  } else if (closestTrackDR) {
    GlobalVector closestTrackVector(closestTrackDR->px(), closestTrackDR->py(), closestTrackDR->pz());
    isoPtLepton -= closestTrackVector.perp();
  }
// back to normal sum - S.L. 30/10/2007
  if (isoPtLepton<0) isoPtLepton = 0;
//  isoPtLepton <= 0.01 ? isoPtLepton = -1 : isoPtLepton = log(isoPtLepton);
  return isoPtLepton;
}

