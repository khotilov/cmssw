//
// $Id: classes.h,v 1.19 2007/09/23 11:37:00 futyand Exp $
//
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "Rtypes.h" 
#include "Math/Cartesian3D.h" 
#include "Math/Polar3D.h" 
#include "Math/CylindricalEta3D.h" 
#include "Math/PxPyPzE4D.h" 
#include <boost/cstdint.hpp> 
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" 
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h" 
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectron.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/Common/interface/RefProd.h" 
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h" 
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" 
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtra.h" 
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h" 
#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h" 
#include "DataFormats/GeometryVector/interface/GlobalPoint.h" 
#include "DataFormats/EgammaCandidates/interface/ConvertedPhotonFwd.h" 
#include "DataFormats/EgammaCandidates/interface/SiStripElectron.h"
#include "DataFormats/CLHEP/interface/Migration.h" 
#include "DataFormats/GeometryVector/interface/LocalPoint.h" 
#include "DataFormats/EgammaCandidates/interface/SiStripElectronFwd.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h" 
#include "DataFormats/EgammaCandidates/interface/PhotonIsolationAssociation.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoCollection.h"
#include "DataFormats/EgammaCandidates/interface/PMGsfElectronIsoNumCollection.h"
#include "DataFormats/EgammaCandidates/interface/PhotonPi0DiscriminatorAssociation.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/AssociationMap.h"


namespace {
  namespace {
    reco::PhotonCollection v1;
    edm::Wrapper<reco::PhotonCollection> w1;
    edm::Ref<reco::PhotonCollection> r1;
    edm::RefProd<reco::PhotonCollection> rp1;
    edm::Wrapper<edm::RefVector<reco::PhotonCollection> > rv1;

    reco::ElectronCollection v2;
    edm::Wrapper<reco::ElectronCollection> w2;
    edm::Ref<reco::ElectronCollection> r2;
    edm::RefProd<reco::ElectronCollection> rp2;
    edm::Wrapper<edm::RefVector<reco::ElectronCollection> > rv2;

    reco::PixelMatchElectronCollection v3;
    edm::Wrapper<reco::PixelMatchElectronCollection> w3;
    edm::Ref<reco::PixelMatchElectronCollection> r3;
    edm::RefProd<reco::PixelMatchElectronCollection> rp3;
    edm::Wrapper<edm::RefVector<reco::PixelMatchElectronCollection> > rv3;

    reco::PixelMatchGsfElectronCollection v4;
    edm::Wrapper<reco::PixelMatchGsfElectronCollection> w4;
    edm::Ref<reco::PixelMatchGsfElectronCollection> r4;
    edm::RefProd<reco::PixelMatchGsfElectronCollection> rp4;
    edm::Wrapper<edm::RefVector<reco::PixelMatchGsfElectronCollection> > rv4;

    reco::SiStripElectronCollection v5;
    edm::Wrapper<reco::SiStripElectronCollection> w5;
    edm::Ref<reco::SiStripElectronCollection> r5;
    edm::RefProd<reco::SiStripElectronCollection> rp5;
    edm::Wrapper<edm::RefVector<reco::SiStripElectronCollection> > rv5;

    reco::ConvertedPhotonCollection v6;
    edm::Wrapper<reco::ConvertedPhotonCollection> w6;
    edm::Ref<reco::ConvertedPhotonCollection> r6;
    edm::RefProd<reco::ConvertedPhotonCollection> rp6;
    edm::Wrapper<edm::RefVector<reco::ConvertedPhotonCollection> > rv6;

    reco::PhotonIsolationMap v66;
    edm::Wrapper<reco::PhotonIsolationMap> w66;
    edm::helpers::Key<edm::RefProd<reco::PhotonCollection > > h66;

    reco::ElectronIsolationMap v7;
    edm::Wrapper<reco::ElectronIsolationMap> w7;
    edm::helpers::Key<edm::RefProd<reco::ElectronCollection > > h7;
    reco::PMGsfElectronIsoCollectionBase b8;
    reco::PMGsfElectronIsoCollection v8;
    reco::PMGsfElectronIsoCollectionRef r8;
    reco::PMGsfElectronIsoCollectionRefProd rp8;
    reco::PMGsfElectronIsoCollectionRefVector rv8;

    edm::Wrapper<reco::PMGsfElectronIsoCollection> w8;

    reco::PMGsfElectronIsoNumCollectionBase b9;
    reco::PMGsfElectronIsoNumCollection v9;
    reco::PMGsfElectronIsoNumCollectionRef r9;
    reco::PMGsfElectronIsoNumCollectionRefProd rp9;
    reco::PMGsfElectronIsoNumCollectionRefVector rv9;

    edm::Wrapper<reco::PMGsfElectronIsoNumCollection> w9;

    reco::PhotonPi0DiscriminatorAssociationMap v10;
    edm::Wrapper<reco::PhotonPi0DiscriminatorAssociationMap> w10;
    edm::helpers::Key<edm::RefProd<reco::PhotonCollection > > h10;

    edm::reftobase::Holder<reco::Candidate, reco::ElectronRef> rb1;
    edm::reftobase::Holder<reco::Candidate, reco::PhotonRef> rb2;
    edm::reftobase::Holder<reco::Candidate, reco::SiStripElectronRef> rb3;
    edm::reftobase::Holder<reco::Candidate, reco::ConvertedPhotonRef> rb4;

    edm::reftobase::Holder<reco::Candidate, reco::PixelMatchGsfElectronRef> rb11;
    edm::reftobase::RefHolder<reco::PixelMatchGsfElectronRef> rb12;
    edm::reftobase::VectorHolder<reco::Candidate, reco::PixelMatchGsfElectronRefVector> rb13;
    edm::reftobase::RefVectorHolder<reco::PixelMatchGsfElectronRefVector> rb14;

    edm::reftobase::Holder<reco::Candidate, reco::PhotonRef> rb21;
    edm::reftobase::RefHolder<reco::PhotonRef> rb22;
    edm::reftobase::VectorHolder<reco::Candidate, reco::PhotonRefVector> rb23;
    edm::reftobase::RefVectorHolder<reco::PhotonRefVector> rb24;

    edm::reftobase::Holder<reco::Candidate, reco::ConvertedPhotonRef> rb31;
    edm::reftobase::RefHolder<reco::ConvertedPhotonRef> rb32;
    edm::reftobase::VectorHolder<reco::Candidate, reco::ConvertedPhotonRefVector> rb33;
    edm::reftobase::RefVectorHolder<reco::ConvertedPhotonRefVector> rb34;
  }
}
