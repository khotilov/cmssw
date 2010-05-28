// -*- C++ -*-
//

// system include files
#include "TEveTrack.h"

// user include files
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "Fireworks/ParticleFlow/interface/setTrackTypePF.h"


//
// constants, enums and typedefs
//
//
// static data member definitions
//
namespace fireworks {

   void setTrackTypePF(const reco::PFCandidate& pfCand,
		       TEveTrack* track ) {

     using namespace reco;

     switch (pfCand.particleId() ) {
     case PFCandidate::e: track->SetLineStyle(5); break;
     case PFCandidate::mu: track->SetLineStyle(6); break;
     case PFCandidate::h0: track->SetLineStyle(3); break;
     case PFCandidate::gamma:  track->SetLineStyle(7); break;
     default: break;
     }
   }
}
