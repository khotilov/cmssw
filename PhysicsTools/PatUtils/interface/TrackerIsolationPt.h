//
// $Id: TrackerIsolationPt.h,v 1.2 2008/01/21 16:26:19 lowette Exp $
//

#ifndef PhysicsTools_PatUtils_TrackerIsolationPt_h
#define PhysicsTools_PatUtils_TrackerIsolationPt_h

/**
  \class    TrackerIsolationPt TrackerIsolationPt.h "PhysicsTools/PatUtils/interface/TrackerIsolationPt.h"
  \brief    Calculates a lepton's tracker isolation pt

   TrackerIsolationPt calculates a tracker isolation pt in a cone
   around the lepton's direction, without doing track extrapolation

  \author   Steven Lowette
  \version  $Id: TrackerIsolationPt.h,v 1.2 2008/01/21 16:26:19 lowette Exp $
*/

namespace reco {
  class Track;
}

namespace edm {
  template<typename T> class View;
  class InputTag;
}

namespace pat {
  class Electron;
  class Muon;
  class TrackerIsolationPt {
  public:
    TrackerIsolationPt();
    virtual ~TrackerIsolationPt();
    
    float calculate(const Electron & theElectron, const edm::View<reco::Track> & theTracks, float isoConeElectron = 0.3) const;
    float calculate(const Muon & theMuon, const edm::View<reco::Track> & theTracks, float isoConeMuon = 0.3) const;
    
  private:
    float calculate(const reco::Track & theTrack, const edm::View<reco::Track> & theTracks, float isoCone) const;
  };

}

#endif
