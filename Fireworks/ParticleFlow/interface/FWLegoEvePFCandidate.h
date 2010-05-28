#ifndef Fireworks_ParticleFlow_FWLegoEvePFCandidate_h
#define Fireworks_ParticleFlow_FWLegoEvePFCandidate_h


#include "TEveLine.h"
#include "TEveStraightLineSet.h"

#include "Rtypes.h"



class TEveTrack;

namespace reco {
  class PFCandidate;
}

class FWLegoEvePFCandidate : public TEveStraightLineSet {

 public:
  FWLegoEvePFCandidate(const reco::PFCandidate& pfc);

 private:

};

#endif
