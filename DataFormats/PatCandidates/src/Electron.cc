//
// $Id: Electron.cc,v 1.13 2008/11/25 08:58:00 fronga Exp $
//

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace pat;


/// default constructor
Electron::Electron() :
    Lepton<reco::GsfElectron>(),
    embeddedGsfTrack_(false),
    embeddedSuperCluster_(false),
    embeddedTrack_(false),
    scSigmaEtaEta_(0),
    scSigmaIEtaIEta_(0), 
    scE1x5_(0),
    scE2x5Max_(0), 
    scE5x5_(0) 
{
}


/// constructor from reco::GsfElectron
Electron::Electron(const reco::GsfElectron & anElectron) :
    Lepton<reco::GsfElectron>(anElectron),
    embeddedGsfTrack_(false),
    embeddedSuperCluster_(false),
    embeddedTrack_(false),
    scSigmaEtaEta_(0),
    scSigmaIEtaIEta_(0), 
    scE1x5_(0),
    scE2x5Max_(0), 
    scE5x5_(0) 
{
}


/// constructor from ref to reco::GsfElectron
Electron::Electron(const edm::RefToBase<reco::GsfElectron> & anElectronRef) :
    Lepton<reco::GsfElectron>(anElectronRef),
    embeddedGsfTrack_(false),
    embeddedSuperCluster_(false),
    embeddedTrack_(false),
    scSigmaEtaEta_(0),
    scSigmaIEtaIEta_(0), 
    scE1x5_(0),
    scE2x5Max_(0), 
    scE5x5_(0) 
{
}

/// constructor from Ptr to reco::GsfElectron
Electron::Electron(const edm::Ptr<reco::GsfElectron> & anElectronRef) :
    Lepton<reco::GsfElectron>(anElectronRef),
    embeddedGsfTrack_(false),
    embeddedSuperCluster_(false),
    embeddedTrack_(false),
    scSigmaEtaEta_(0),
    scSigmaIEtaIEta_(0), 
    scE1x5_(0),
    scE2x5Max_(0), 
    scE5x5_(0) 
{
}


/// destructor
Electron::~Electron() {
}


/// override the reco::GsfElectron::gsfTrack method, to access the internal storage of the supercluster
reco::GsfTrackRef Electron::gsfTrack() const {
  if (embeddedGsfTrack_) {
    return reco::GsfTrackRef(&gsfTrack_, 0);
  } else {
    return reco::GsfElectron::gsfTrack();
  }
}


/// override the reco::GsfElectron::superCluster method, to access the internal storage of the supercluster
reco::SuperClusterRef Electron::superCluster() const {
  if (embeddedSuperCluster_) {
    return reco::SuperClusterRef(&superCluster_, 0);
  } else {
    return reco::GsfElectron::superCluster();
  }
}


/// override the reco::GsfElectron::track method, to access the internal storage of the track
reco::TrackRef Electron::track() const {
  if (embeddedTrack_) {
    return reco::TrackRef(&track_, 0);
  } else {
    return reco::GsfElectron::track();
  }
}

/// method to store the electron's gsfTrack internally
void Electron::embedGsfTrack() {
  gsfTrack_.clear();
  if (reco::GsfElectron::gsfTrack().isNonnull()) {
      gsfTrack_.push_back(*reco::GsfElectron::gsfTrack());
      embeddedGsfTrack_ = true;
  }
}


/// method to store the electron's supercluster internally
void Electron::embedSuperCluster() {
  superCluster_.clear();
  if (reco::GsfElectron::superCluster().isNonnull()) {
      superCluster_.push_back(*reco::GsfElectron::superCluster());
      embeddedSuperCluster_ = true;
  }
}


/// method to store the electron's track internally
void Electron::embedTrack() {
  track_.clear();
  if (reco::GsfElectron::track().isNonnull()) {
      track_.push_back(*reco::GsfElectron::track());
      embeddedTrack_ = true;
  }
}

// method to retrieve a lepton ID (or throw)
float Electron::electronID(const std::string & name) const {
    for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
        if (it->first == name) return it->second;
    }
    cms::Exception ex("Key not found");
    ex << "pat::Electron: the ID " << name << " can't be found in this pat::Electron.\n";
    ex << "The available IDs are: ";
    for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
        ex << "'" << it->first << "' ";
    }
    ex << ".\n";
    throw ex;
}
// check if an ID is there
bool Electron::isElectronIDAvailable(const std::string & name) const {
    for (std::vector<IdPair>::const_iterator it = electronIDs_.begin(), ed = electronIDs_.end(); it != ed; ++it) {
        if (it->first == name) return true;
    }
    return false;
}
/// method to store the electron's cluster shape
void Electron::setClusterShapes (const float& scSigmaEtaEta, const float& scSigmaIEtaIEta, 
                                 const float& scE1x5, const float& scE2x5Max, const float& scE5x5) 
  { 
    scSigmaEtaEta_ = scSigmaEtaEta ; 
    scSigmaIEtaIEta_ = scSigmaIEtaIEta ;
    scE1x5_ = scE1x5 ;
    scE2x5Max_ = scE2x5Max ;
    scE5x5_ = scE5x5 ;
  }


