#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositeRefBaseCandidate.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandMatchMapMany.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"

namespace {
  namespace {
    std::vector<reco::Candidate *> v1;
    reco::CandidateCollection o1;
    edm::Wrapper<reco::CandidateCollection> w1;
    std::vector<reco::Particle> v2;
    edm::Wrapper<std::vector<reco::Particle> > w2;
    reco::CandidateRef r1;
    reco::CandidateBaseRef r2;
    std::vector<reco::CandidateBaseRef> rv2;
    edm::reftobase::IndirectHolder<reco::Candidate> rbih1;
    edm::reftobase::RefHolder<reco::CandidateRef> rh1;
    edm::reftobase::IndirectVectorHolder<reco::Candidate> rbih2;
    edm::reftobase::RefVectorHolder<reco::CandidateRefVector> rh2;
    edm::reftobase::Holder<reco::Candidate, reco::CandidateRef> rhcr1;
    edm::reftobase::VectorHolder<reco::Candidate, reco::CandidateRefVector> rhcr2;
    edm::Wrapper<reco::CandidateRefVector> wrv1;
    edm::Wrapper<reco::CandidateBaseRefVector> wrv2;
    reco::CandidateRefProd rp1;
    reco::CandidateBaseRefProd rp2;
    std::vector<edm::RefToBase<reco::Candidate> > vrb1;
    edm::Wrapper<reco::CandFloatAssociations> wav1;
    edm::Wrapper<reco::CandDoubleAssociations> wav2;
    edm::Wrapper<reco::CandIntAssociations> wav3;
    edm::Wrapper<reco::CandUIntAssociations> wav4;
    edm::Wrapper<reco::CandViewFloatAssociations> wav5;
    edm::Wrapper<reco::CandViewDoubleAssociations> wav6;
    edm::Wrapper<reco::CandViewIntAssociations> wav7;
    edm::Wrapper<reco::CandViewUIntAssociations> wav8;
    edm::helpers::KeyVal<reco::CandidateRef,reco::CandidateRef> kv1;
    reco::CandMatchMap cmm1;
    reco::CandMatchMap::const_iterator cmm1it;
    edm::Wrapper<reco::CandMatchMap> wcmm1;
    edm::helpers::KeyVal<reco::CandidateRefProd, reco::CandidateRefProd> kv2;
    reco::CandViewMatchMap cmm2;
    reco::CandViewMatchMap::const_iterator cmm2it;
    edm::Wrapper<reco::CandViewMatchMap> wcmm2;
    edm::helpers::KeyVal<reco::CandidateBaseRefProd, reco::CandidateBaseRefProd> kv3;
    std::map<const reco::Candidate *, const reco::Candidate *> m1;
    std::vector<const reco::Candidate *> vc1;
    reco::CandMatchMapMany cmm3;
    reco::CandMatchMapMany::const_iterator cmm3it;
    edm::Wrapper<reco::CandMatchMapMany> wcmm3;
    edm::Wrapper<std::vector<reco::CandidateBaseRef> > wvrb1;
  }
}
