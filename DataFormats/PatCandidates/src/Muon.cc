//
// $Id: Muon.cc,v 1.11 2008/08/31 20:24:54 hegner Exp $
//

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidate.h"

using namespace pat;


/// default constructor
Muon::Muon() :
    Lepton<MuonType>(),
    embeddedTrack_(false),
    embeddedStandAloneMuon_(false),
    embeddedCombinedMuon_(false),
    embeddedPFCandidate_(false)
{
}


/// constructor from MuonType
Muon::Muon(const MuonType & aMuon) :
    Lepton<MuonType>(aMuon),
    embeddedTrack_(false),
    embeddedStandAloneMuon_(false),
    embeddedCombinedMuon_(false),
    embeddedPFCandidate_(false)
{
}


/// constructor from ref to MuonType
Muon::Muon(const edm::RefToBase<MuonType> & aMuonRef) :
    Lepton<MuonType>(aMuonRef),
    embeddedTrack_(false),
    embeddedStandAloneMuon_(false),
    embeddedCombinedMuon_(false),
    embeddedPFCandidate_(false)
{
}


/// constructor from ref to MuonType
Muon::Muon(const edm::Ptr<MuonType> & aMuonRef) :
    Lepton<MuonType>(aMuonRef),
    embeddedTrack_(false),
    embeddedStandAloneMuon_(false),
    embeddedCombinedMuon_(false),
    embeddedPFCandidate_(false)
{
}


/// destructor
Muon::~Muon() {
}


/// reference to Track reconstructed in the tracker only (reimplemented from reco::Muon)
reco::TrackRef Muon::track() const {
  if (embeddedTrack_) {
    return reco::TrackRef(&track_, 0);
  } else {
    return MuonType::innerTrack();
  }
}


/// reference to Track reconstructed in the muon detector only (reimplemented from reco::Muon)
reco::TrackRef Muon::standAloneMuon() const {
  if (embeddedStandAloneMuon_) {
    return reco::TrackRef(&standAloneMuon_, 0);
  } else {
    return MuonType::outerTrack();
  }
}


/// reference to Track reconstructed in both tracked and muon detector (reimplemented from reco::Muon)
reco::TrackRef Muon::combinedMuon() const {
  if (embeddedCombinedMuon_) {
    return reco::TrackRef(&combinedMuon_, 0);
  } else {
    return MuonType::globalTrack();
  }
}

reco::IsolatedPFCandidateRef Muon::pfCandidateRef() const {
  if (embeddedPFCandidate_) {
    return reco::IsolatedPFCandidateRef(&pfCandidate_, 0);
  } else {
    return pfCandidateRef_;
  }
}


/// return the lepton ID discriminator
float Muon::leptonID() const {
  return leptonID_;
}


/// return the muon segment compatibility -> meant for
float Muon::segmentCompatibility() const {
  return muon::segmentCompatibility(*this);
}


/// return whether it is a good muon
bool Muon::isGoodMuon(const MuonType & muon, reco::Muon::SelectionType type) {
  return muon::isGoodMuon(*this, type);
}


/// embed the Track reconstructed in the tracker only
void Muon::embedTrack() {
  track_.clear();
  if (MuonType::innerTrack().isNonnull()) {
      track_.push_back(*MuonType::innerTrack());
      embeddedTrack_ = true;
  }
}


/// embed the Track reconstructed in the muon detector only
void Muon::embedStandAloneMuon() {
  standAloneMuon_.clear();
  if (MuonType::outerTrack().isNonnull()) {
      standAloneMuon_.push_back(*MuonType::outerTrack());
      embeddedStandAloneMuon_ = true;
  }
}


/// embed the Track reconstructed in both tracked and muon detector
void Muon::embedCombinedMuon() {
  combinedMuon_.clear();
  if (MuonType::globalTrack().isNonnull()) {
      combinedMuon_.push_back(*MuonType::globalTrack());
      embeddedCombinedMuon_ = true;
  }
}

void Muon::embedPFCandidate() {
  pfCandidate_.clear();
  if ( pfCandidateRef_.isAvailable() && pfCandidateRef_.isNonnull()) {
    pfCandidate_.push_back( *pfCandidateRef_ );
    embeddedPFCandidate_ = true;
  }
}

/// method to set the lepton ID discriminator
void Muon::setLeptonID(float id) {
  leptonID_ = id;
}
