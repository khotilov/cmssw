//
// $Id: Electron.h,v 1.26 2009/10/12 22:18:50 rwolf Exp $
//

#ifndef DataFormats_PatCandidates_Electron_h
#define DataFormats_PatCandidates_Electron_h

/**
  \class    pat::Electron Electron.h "DataFormats/PatCandidates/interface/Electron.h"
  \brief    Analysis-level electron class

   pat::Electron implements the analysis-level electron class within the
   'pat' namespace.

   Please post comments and questions to the Physics Tools hypernews:
   https://hypernews.cern.ch/HyperNews/CMS/get/physTools.html

  \author   Steven Lowette, Giovanni Petrucciani, Frederic Ronga
  \version  $Id: Electron.h,v 1.26 2009/10/12 22:18:50 rwolf Exp $
*/


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// Define typedefs for convenience
namespace pat {
  class Electron;
  typedef std::vector<Electron>              ElectronCollection; 
  typedef edm::Ref<ElectronCollection>       ElectronRef; 
  typedef edm::RefVector<ElectronCollection> ElectronRefVector; 
}


// Class definition
namespace pat {


  class Electron : public Lepton<reco::GsfElectron> {

    public:

      typedef std::pair<std::string,float> IdPair;

      /// default constructor
      Electron();
      /// constructor from a reco electron
      Electron(const reco::GsfElectron & anElectron);
      /// constructor from a RefToBase to a reco electron (to be superseded by Ptr counterpart)
      Electron(const edm::RefToBase<reco::GsfElectron> & anElectronRef);
      /// constructor from a Ptr to a reco electron
      Electron(const edm::Ptr<reco::GsfElectron> & anElectronRef);
      /// destructor
      virtual ~Electron();

      /// required reimplementation of the Candidate's clone method
      virtual Electron * clone() const { return new Electron(*this); }

      // ---- methods for content embedding ----
      /// override the reco::GsfElectron::gsfTrack method, to access the internal storage of the supercluster
      reco::GsfTrackRef gsfTrack() const;
      /// override the reco::GsfElectron::superCluster method, to access the internal storage of the supercluster
      reco::SuperClusterRef superCluster() const;
      /// override the reco::GsfElectron::track method, to access the internal storage of the track
      reco::TrackRef track() const;
      using reco::RecoCandidate::track; // avoid hiding the base implementation
      /// method to store the electron's GsfTrack internally
      void embedGsfTrack();
      /// method to store the electron's SuperCluster internally
      void embedSuperCluster();
      /// method to store the electron's Track internally
      void embedTrack();

      // ---- methods for electron ID ----
      /// Returns a specific electron ID associated to the pat::Electron given its name
      /// For cut-based IDs, the value is 1.0 for good, 0.0 for bad.
      /// Note: an exception is thrown if the specified ID is not available
      float electronID(const std::string & name) const;
      /// Returns true if a specific ID is available in this pat::Electron
      bool isElectronIDAvailable(const std::string & name) const;
      /// Returns all the electron IDs in the form of <name,value> pairs
      /// The 'default' ID is the first in the list
      const std::vector<IdPair> &  electronIDs() const { return electronIDs_; }
      /// Store multiple electron ID values, discarding existing ones
      /// The first one in the list becomes the 'default' electron id 
      void setElectronIDs(const std::vector<IdPair> & ids) { electronIDs_ = ids; }
      
      // ---- overload of isolation functions ----
      /// Overload of pat::Lepton::trackIso(); returns the value of
      /// the summed track pt in a cone of deltaR<0.4
      float trackIso() const { return dr04TkSumPt(); }
      /// Overload of pat::Lepton::trackIso(); returns the value of 
      /// the summed Et of all recHits in the ecal in a cone of 
      /// deltaR<0.4
      float ecalIso()  const { return dr04EcalRecHitSumEt(); }
      /// Overload of pat::Lepton::trackIso(); returns the value of 
      /// the summed Et of all caloTowers in the hcal in a cone of 
      /// deltaR<0.4
      float hcalIso()  const { return dr04HcalTowerSumEt(); }
      /// Overload of pat::Lepton::trackIso(); returns the sum of 
      /// ecalIso() and hcalIso
      float caloIso()  const { return ecalIso()+hcalIso(); }

      // ---- PF specific methods ----
      /// reference to the source PFCandidates
      /// null if this has been built from a standard electron
      reco::PFCandidateRef pfCandidateRef() const;
      /// add a reference to the source IsolatedPFCandidate
      void setPFCandidateRef(const reco::PFCandidateRef& ref) {
	pfCandidateRef_ = ref;
      } 
      /// embed the PFCandidate pointed to by pfCandidateRef_
      void embedPFCandidate();

      /// dB gives the impact parameter wrt the beamline.
      /// If this is not cached it is not meaningful, since
      /// it relies on the distance to the beamline. 
      double dB() const;
      void setDB(double dB) { dB_ = dB; cachedDB_ = true;}

    protected:

      // ---- for content embedding ----
      bool embeddedGsfTrack_;
      std::vector<reco::GsfTrack> gsfTrack_;
      bool embeddedSuperCluster_;
      std::vector<reco::SuperCluster> superCluster_;
      bool embeddedTrack_;
      std::vector<reco::Track> track_;
      // ---- electron ID's holder ----
      std::vector<IdPair> electronIDs_;

      // ---- PF specific members ----
      /// true if the IsolatedPFCandidate is embedded
      bool embeddedPFCandidate_;      
      /// if embeddedPFCandidate_, a copy of the source IsolatedPFCandidate
      /// is stored in this vector
      reco::PFCandidateCollection pfCandidate_;
      /// reference to the IsolatedPFCandidate this has been built from
      /// null if this has been built from a standard electron
      reco::PFCandidateRef pfCandidateRef_;

      // V+Jets group selection variables. 
      bool    cachedDB_;         // have these values been cached? 
      double  dB_;               // track->dxy( beamPoint ) 
  };
}

#endif
