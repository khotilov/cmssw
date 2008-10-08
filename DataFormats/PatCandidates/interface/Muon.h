//
// $Id: Muon.h,v 1.17 2008/10/07 18:15:13 lowette Exp $
//

#ifndef DataFormats_PatCandidates_Muon_h
#define DataFormats_PatCandidates_Muon_h

/**
  \class    pat::Muon Muon.h "DataFormats/PatCandidates/interface/Muon.h"
  \brief    Analysis-level muon class

   pat::Muon implements the analysis-level muon class within the 'pat'
   namespace.

   Please post comments and questions to the Physics Tools hypernews:
   https://hypernews.cern.ch/HyperNews/CMS/get/physTools.html

  \author   Steven Lowette, Giovanni Petrucciani, Frederic Ronga, Colin Bernet
  \version  $Id: Muon.h,v 1.17 2008/10/07 18:15:13 lowette Exp $
*/


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/IsolatedPFCandidate.h"

// Define typedefs for convenience
namespace pat {
  class Muon;
  typedef std::vector<Muon>              MuonCollection; 
  typedef edm::Ref<MuonCollection>       MuonRef; 
  typedef edm::RefVector<MuonCollection> MuonRefVector; 
}


// Class definition
namespace pat {


  typedef reco::Muon MuonType;
  typedef reco::MuonCollection MuonTypeCollection;


  class Muon : public Lepton<MuonType> {

    public:

      /// default constructor
      Muon();
      /// constructor from a reco muon
      Muon(const MuonType & aMuon);
      /// constructor from a RefToBase to a reco muon (to be superseded by Ptr counterpart)
      Muon(const edm::RefToBase<MuonType> & aMuonRef);
      /// constructor from a Ptr to a reco muon
      Muon(const edm::Ptr<MuonType> & aMuonRef);
      /// destructor
      virtual ~Muon();

      /// required reimplementation of the Candidate's clone method
      virtual Muon * clone() const { return new Muon(*this); }

      // ---- methods for content embedding ----
      /// reference to Track reconstructed in the tracker only (reimplemented from reco::Muon)
      reco::TrackRef track() const;
      /// reference to Track reconstructed in the tracker only (reimplemented from reco::Muon)
      reco::TrackRef innerTrack() const { return track(); }
      /// reference to Track reconstructed in the muon detector only (reimplemented from reco::Muon)
      reco::TrackRef standAloneMuon() const;
      /// reference to Track reconstructed in the muon detector only (reimplemented from reco::Muon)
      reco::TrackRef outerTrack() const { return standAloneMuon(); }
      /// reference to Track reconstructed in both tracked and muon detector (reimplemented from reco::Muon)
      reco::TrackRef combinedMuon() const;
      /// reference to Track reconstructed in both tracked and muon detector (reimplemented from reco::Muon)
      reco::TrackRef globalTrack() const { return combinedMuon(); }
      /// set reference to Track reconstructed in the tracker only (reimplemented from reco::Muon)
      void embedTrack();
      /// set reference to Track reconstructed in the muon detector only (reimplemented from reco::Muon)
      void embedStandAloneMuon();
      /// set reference to Track reconstructed in both tracked and muon detector (reimplemented from reco::Muon)
      void embedCombinedMuon();

      // ---- PF specific methods ----
      /// reference to the source IsolatedPFCandidates
      /// null if this has been built from a standard muon
      reco::IsolatedPFCandidateRef pfCandidateRef() const;
      /// add a reference to the source IsolatedPFCandidate
      void setPFCandidateRef(const reco::IsolatedPFCandidateRef& ref) {
	pfCandidateRef_ = ref;
      } 
      /// embed the IsolatedPFCandidate pointed to by pfCandidateRef_
      void embedPFCandidate();

      // ---- methods for muon ID ----
      /// return whether it is a good muon
      bool isGoodMuon(const MuonType & muon, reco::Muon::SelectionType type = reco::Muon::TMLastStationLoose);
      /// return the muon segment compatibility
      float segmentCompatibility() const;
      /// return the lepton ID discriminator
      float leptonID() const;
      /// method to set the lepton ID discriminator
      void setLeptonID(float id);

    protected:

      // ---- for content embedding ----
      bool embeddedTrack_;
      std::vector<reco::Track> track_;
      bool embeddedStandAloneMuon_;
      std::vector<reco::Track> standAloneMuon_;
      bool embeddedCombinedMuon_;
      std::vector<reco::Track> combinedMuon_;
      // ---- PF specific members ----
      /// true if the IsolatedPFCandidate is embedded
      bool embeddedPFCandidate_;      
      /// if embeddedPFCandidate_, a copy of the source IsolatedPFCandidate
      /// is stored in this vector
      reco::IsolatedPFCandidateCollection pfCandidate_;
      /// reference to the IsolatedPFCandidate this has been built from
      /// null if this has been built from a standard muon
      reco::IsolatedPFCandidateRef pfCandidateRef_;
      // ---- muon ID's holder ----
      float leptonID_;

  };


}

#endif
