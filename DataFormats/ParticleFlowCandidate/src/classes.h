#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtra.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"
#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"

#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateElectronExtraFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtraFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/ValueMap.h"

namespace {
  struct dictionary {

    reco::PFCandidateRef c_r;
    reco::PFCandidateRefProd c_rp;
    reco::PFCandidateRefVector c_rv;
    edm::Wrapper<std::vector<reco::PFCandidate> > w1;
    edm::reftobase::Holder<reco::Candidate, reco::PFCandidateRef> bla1; 
    edm::reftobase::RefHolder<reco::PFCandidateRef> bla2; 
    reco::PFCandidatePtr bla3;     
    std::vector<reco::PFCandidatePtr> bla33;     
    edm::Wrapper<edm::ValueMap<edm::Ref<std::vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::PFCandidate>,reco::PFCandidate> > > > bla333;
    edm::ValueMap<edm::Ref<std::vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::PFCandidate>,reco::PFCandidate> > >  bla334;
    std::vector<edm::Ref<std::vector<reco::PFCandidate>,reco::PFCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::PFCandidate>,reco::PFCandidate> > > bla335;

    reco::PFCandidate::ElementInBlock jo1;
    reco::PFCandidate::ElementsInBlocks jo2;  

    reco::IsolatedPFCandidateRef ic_r;
    reco::IsolatedPFCandidateRefProd ic_rp;
    reco::IsolatedPFCandidateRefVector ic_rv;
    edm::Wrapper<std::vector<reco::IsolatedPFCandidate> > iw1;
    reco::IsolatedPFCandidatePtr bla4;     

    reco::PileUpPFCandidateRef puc_r;
    reco::PileUpPFCandidateRefProd puc_rp;
    reco::PileUpPFCandidateRefVector puc_rv;
    edm::Wrapper<std::vector<reco::PileUpPFCandidate> > puw1;
    reco::PileUpPFCandidatePtr bla5;     

    reco::PFCandidateFwdRef c_fr;
    edm::reftobase::Holder<reco::Candidate, reco::PFCandidateFwdRef> bla1f; 
    edm::reftobase::RefHolder<reco::PFCandidateFwdRef> bla2f; 
    reco::PFCandidateFwdPtr fbla3;     
    std::vector<reco::PFCandidateFwdRef> bla33f;     
    std::vector<reco::PFCandidateFwdPtr> bla33fv; 
    
    reco::PFCandidateElectronExtraRef ex_r;
    reco::PFCandidateElectronExtraRefProd ex_rp;
    reco::PFCandidateElectronExtraRefVector ex_rv;
    edm::Wrapper<std::vector<reco::PFCandidateElectronExtra> > ex_w1;

    reco::PFCandidatePhotonExtraRef ph_r;
    reco::PFCandidatePhotonExtraRefProd ph_rp;
    //    reco::PFCandidatePhotonExtraRefVector ph_rv;
    edm::Wrapper<std::vector<reco::PFCandidatePhotonExtra> > ph_w1;

  };
}
