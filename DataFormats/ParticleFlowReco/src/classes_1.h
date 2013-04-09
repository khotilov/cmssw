#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/Ref.h"
#include <DataFormats/Common/interface/OwnVector.h>
#include <DataFormats/Common/interface/ClonePolicy.h>
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "Math/Cartesian3D.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "Math/Polar3D.h"
#include "Math/CylindricalEta3D.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "Math/GenVector/PositionVector3D.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "Rtypes.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "Math/PxPyPzE4D.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <boost/cstdint.hpp>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h"
/* #include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h" */
/* #include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h" */
#include "DataFormats/ParticleFlowReco/interface/PFSimParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFTrajectoryPoint.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementGsfTrack.h"  //Daniele
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementBrem.h"  //Daniele
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h" //Florian
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedTrackerVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFNuclearInteraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversion.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversionFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFV0.h"
#include "DataFormats/ParticleFlowReco/interface/ConvBremSeed.h"
#include "DataFormats/ParticleFlowReco/interface/ConvBremSeedFwd.h"
//Includes by Jamie
#include "DataFormats/ParticleFlowReco/interface/Calibratable.h"
#include "DataFormats/ParticleFlowReco/interface/CalibrationResultWrapper.h"
#include "DataFormats/ParticleFlowReco/interface/CalibrationProvenance.h"
#include "DataFormats/ParticleFlowReco/interface/CaloWindow.h"
#include "DataFormats/ParticleFlowReco/interface/CaloEllipse.h"
#include "DataFormats/ParticleFlowReco/interface/CaloBox.h"
#include "DataFormats/ParticleFlowReco/interface/ParticleFiltrationDecision.h"

#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexSeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexSeedFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFClusterRefCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/RecoPFClusterRefCandidateFwd.h"

#include <map>

namespace {
  struct dictionary {

    /* NuclearInteraction stuffs  */
    reco::PFNuclearInteraction                                dummy21;
    std::vector<reco::PFNuclearInteraction>                   dummy22;
    edm::Wrapper<std::vector<reco::PFNuclearInteraction> >    dummy23;
    edm::Ref<std::vector<reco::PFNuclearInteraction> >        dummy24;
    edm::RefProd<std::vector<reco::PFNuclearInteraction> >    dummy25;
    edm::RefVector<std::vector<reco::PFNuclearInteraction> >  dummy26;

    reco::PFDisplacedTrackerVertex                                dummy21a;
    std::vector<reco::PFDisplacedTrackerVertex>                   dummy22a;
    edm::Wrapper<std::vector<reco::PFDisplacedTrackerVertex> >    dummy23a;
    edm::Ref<std::vector<reco::PFDisplacedTrackerVertex> >        dummy24a;
    edm::RefProd<std::vector<reco::PFDisplacedTrackerVertex> >    dummy25a;
    edm::RefVector<std::vector<reco::PFDisplacedTrackerVertex> >  dummy26a;
    reco::PFConversionCollection dummy27;
    edm::Wrapper<reco::PFConversionCollection> dummy28;
    edm::Ref<reco::PFConversionCollection> dummy29;
    edm::RefProd<reco::PFConversionCollection> dummy30;
    edm::Wrapper<edm::RefVector<reco::PFConversionCollection> > dummy31;
    std::vector<edm::Ref<std::vector<reco::PFRecTrack>,reco::PFRecTrack,edm::refhelper::FindUsingAdvance<std::vector<reco::PFRecTrack>,reco::PFRecTrack> > > dummy32;

    /* V0 stuffs  */
    reco::PFV0                                dummy33;
    std::vector<reco::PFV0>                   dummy34;
    edm::Wrapper<std::vector<reco::PFV0> >    dummy35;
    edm::Ref<std::vector<reco::PFV0> >        dummy36;
    edm::RefProd<std::vector<reco::PFV0> >    dummy37;
    edm::RefVector<std::vector<reco::PFV0> >  dummy38;
    edm::Ref<std::vector<reco::VertexCompositeCandidate>,reco::VertexCompositeCandidate,edm::refhelper::FindUsingAdvance<std::vector<reco::VertexCompositeCandidate>,reco::VertexCompositeCandidate> > dummy39;

    /* ConvBremSeed stuff */
    reco::ConvBremSeedCollection dummy40;
    edm::Wrapper<reco::ConvBremSeedCollection> dummy41;
    edm::Ref<reco::ConvBremSeedCollection> dummy42;
    edm::RefProd<reco::ConvBremSeedCollection> dummy43;
    edm::Wrapper<edm::RefVector<reco::ConvBremSeedCollection> > dummy44;
    edm::RefToBase<reco::ConvBremSeed> dummy45;
    edm::reftobase::Holder< reco::ConvBremSeed, edm::Ref<reco::ConvBremSeedCollection> > dummy46;
    edm::reftobase::RefHolder< edm::Ref<reco::ConvBremSeedCollection> > dummy47;


  };
}
