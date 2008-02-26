#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/Wrapper.h"

namespace {
  namespace {
    reco::PFCandidateRef c_r;
    reco::PFCandidateRefProd c_rp;
    reco::PFCandidateRefVector c_rv;
    edm::Wrapper<std::vector<reco::PFCandidate> > w1;
    //Needed since RefToBase is there (from 17x on)
    edm::reftobase::Holder<reco::Candidate, reco::PFCandidateRef> bla1; 
    edm::reftobase::RefHolder<reco::PFCandidateRef> bla2; 
    //This Dummies are needed  
    reco::PFCandidate::ElementInBlock jo1;
    reco::PFCandidate::ElementsInBlocks jo2;  
  }
}
