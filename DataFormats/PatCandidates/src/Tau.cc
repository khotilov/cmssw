//
// $Id: Tau.cc,v 1.7 2008/06/09 09:03:19 gpetrucc Exp $
//

#include "DataFormats/PatCandidates/interface/Tau.h"


using namespace pat;


/// default constructor
Tau::Tau() :
    Lepton<TauType>(),
    embeddedIsolationTracks_(false),
    embeddedLeadTrack_(false),
    embeddedSignalTracks_(false)
{
}


/// constructor from TauType
Tau::Tau(const TauType & aTau) :
    Lepton<TauType>(aTau),
    embeddedIsolationTracks_(false),
    embeddedLeadTrack_(false),
    embeddedSignalTracks_(false)
{
    const reco::PFTau * pfTau = dynamic_cast<const reco::PFTau *>(&aTau);
    if (pfTau != 0) pfSpecific_.push_back(pat::tau::TauPFSpecific(*pfTau));
    const reco::CaloTau * caloTau = dynamic_cast<const reco::CaloTau *>(&aTau);
    if (caloTau != 0) caloSpecific_.push_back(pat::tau::TauCaloSpecific(*caloTau));
}


/// constructor from ref to TauType
Tau::Tau(const edm::RefToBase<TauType> & aTauRef) :
    Lepton<TauType>(aTauRef),
    embeddedIsolationTracks_(false),
    embeddedLeadTrack_(false),
    embeddedSignalTracks_(false)
{
    const reco::PFTau * pfTau = dynamic_cast<const reco::PFTau *>(aTauRef.get());
    if (pfTau != 0) pfSpecific_.push_back(pat::tau::TauPFSpecific(*pfTau));
    const reco::CaloTau * caloTau = dynamic_cast<const reco::CaloTau *>(aTauRef.get());
    if (caloTau != 0) caloSpecific_.push_back(pat::tau::TauCaloSpecific(*caloTau));
}

/// constructor from ref to TauType
Tau::Tau(const edm::Ptr<TauType> & aTauRef) :
    Lepton<TauType>(aTauRef),
    embeddedIsolationTracks_(false),
    embeddedLeadTrack_(false),
    embeddedSignalTracks_(false)
{
    const reco::PFTau * pfTau = dynamic_cast<const reco::PFTau *>(aTauRef.get());
    if (pfTau != 0) pfSpecific_.push_back(pat::tau::TauPFSpecific(*pfTau));
    const reco::CaloTau * caloTau = dynamic_cast<const reco::CaloTau *>(aTauRef.get());
    if (caloTau != 0) caloSpecific_.push_back(pat::tau::TauCaloSpecific(*caloTau));
}



/// destructor
Tau::~Tau() {
}


/// override the TauType::isolationTracks method, to access the internal storage of the track
reco::TrackRefVector Tau::isolationTracks() const {
  if (embeddedIsolationTracks_) {
    reco::TrackRefVector trackRefVec;
    for (unsigned int i = 0; i < isolationTracks_.size(); i++) {
      trackRefVec.push_back(reco::TrackRef(&isolationTracks_, i));
    }
    return trackRefVec;
  } else {
    return TauType::isolationTracks();
  }
}


/// override the TauType::track method, to access the internal storage of the track
reco::TrackRef Tau::leadTrack() const {
  if (embeddedLeadTrack_) {
    return reco::TrackRef(&leadTrack_, 0);
  } else {
    return TauType::leadTrack();
  }
}


/// override the TauType::track method, to access the internal storage of the track
reco::TrackRefVector Tau::signalTracks() const {
  if (embeddedSignalTracks_) {
    reco::TrackRefVector trackRefVec;
    for (unsigned int i = 0; i < signalTracks_.size(); i++) {
      trackRefVec.push_back(reco::TrackRef(&signalTracks_, i));
    }
    return trackRefVec;
  } else {
    return TauType::signalTracks();
  }
}


/// method to store the isolation tracks internally
void Tau::embedIsolationTracks() {
  isolationTracks_.clear();
  reco::TrackRefVector trackRefVec = TauType::isolationTracks();
  for (unsigned int i = 0; i < trackRefVec.size(); i++) {
    isolationTracks_.push_back(*trackRefVec.at(i));
  }
  embeddedIsolationTracks_ = true;
}


/// method to store the isolation tracks internally
void Tau::embedLeadTrack() {
  leadTrack_.clear();
  leadTrack_.push_back(*TauType::leadTrack());
  embeddedLeadTrack_ = true;
}


/// method to store the isolation tracks internally
void Tau::embedSignalTracks(){
  signalTracks_.clear();
  reco::TrackRefVector trackRefVec = TauType::signalTracks();
  for (unsigned int i = 0; i < trackRefVec.size(); i++) {
    signalTracks_.push_back(*trackRefVec.at(i));
  }
  embeddedSignalTracks_ = true;
}

const pat::tau::TauPFSpecific & Tau::pfSpecific() const {
  if (!isPFTau()) throw cms::Exception("Type Error") << "Requesting a PFTau-specific information from a pat::Tau which wasn't made from a PFTau.\n";
  return pfSpecific_[0]; 
}

const pat::tau::TauCaloSpecific & Tau::caloSpecific() const {
  if (!isCaloTau()) throw cms::Exception("Type Error") << "Requesting a CaloTau-specific information from a pat::Tau which wasn't made from a CaloTau.\n";
  return caloSpecific_[0]; 
}
