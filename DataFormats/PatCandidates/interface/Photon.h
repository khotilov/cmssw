//
// $Id: Photon.h,v 1.14 2008/07/08 20:56:48 gpetrucc Exp $
//

#ifndef DataFormats_PatCandidates_Photon_h
#define DataFormats_PatCandidates_Photon_h

/**
  \class    pat::Photon Photon.h "DataFormats/PatCandidates/interface/Photon.h"
  \brief    Analysis-level Photon class

   pat::Photon implements the analysis-level photon class within the 'pat'
   namespace.

   Please post comments and questions to the Physics Tools hypernews:
   https://hypernews.cern.ch/HyperNews/CMS/get/physTools.html

  \author   Steven Lowette, Giovanni Petrucciani, Frederic Ronga
  \version  $Id: Photon.h,v 1.14 2008/07/08 20:56:48 gpetrucc Exp $
*/

#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/EgammaCandidates/interface/PhotonID.h"


namespace pat {


  typedef reco::Photon PhotonType;


  class Photon : public PATObject<PhotonType> {

    public:

      /// default constructor
      Photon();
      /// constructor from a reco photon
      Photon(const PhotonType & aPhoton);
      /// constructor from a RefToBase to a reco photon (to be superseded by Ptr counterpart)
      Photon(const edm::RefToBase<PhotonType> & aPhotonRef);
      /// constructor from a Ptr to a reco photon
      Photon(const edm::Ptr<PhotonType> & aPhotonRef);
      /// destructor
      virtual ~Photon();

      /// required reimplementation of the Candidate's clone method
      virtual Photon * clone() const { return new Photon(*this); }

      // ---- methods for content embedding ----
      /// override the superCluster method from CaloJet, to access the internal storage of the supercluster
      reco::SuperClusterRef superCluster() const;
      /// method to store the photon's supercluster internally
      void embedSuperCluster();

      // ---- methods for access the generated photon ----
      /// return the match to the generated photon
      const reco::Particle * genPhoton() const { return genParticle(); }
      /// method to set the generated photon
      void setGenPhoton(const reco::GenParticleRef & gp, bool embed=false) { setGenParticleRef(gp, embed); }

      // ---- methods for photon ID ----
      /// returns the PhotonID object, or a null pointer if no ID is available
      const reco::PhotonID * photonID() const { return photonID_.empty() ? 0 : & photonID_[0]; }
      /// sets the PhotonID object
      void setPhotonID(const reco::PhotonID & photonID) { photonID_.clear(); photonID_.push_back(photonID); }
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isLooseEM() const { return photonIDOrThrow().isLooseEM(); }
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isLoosePhoton() const { return photonIDOrThrow().isLoosePhoton(); }
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isTightPhoton() const { return photonIDOrThrow().isTightPhoton(); }
      /// Returns computed EcalRecHit isolation
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      float isolationEcalRecHit() const { return photonIDOrThrow().isolationEcalRecHit(); }
      /// Returns computed HcalRecHit isolation
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      float isolationHcalRecHit() const { return photonIDOrThrow().isolationHcalRecHit(); }
      /// Returns calculated sum track pT cone of dR
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      float isolationSolidTrkCone() const { return photonIDOrThrow().isolationSolidTrkCone(); }
      /// As above, excluding the core at the center of the cone
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      float isolationHollowTrkCone() const { return photonIDOrThrow().isolationHollowTrkCone(); }
      /// Returns number of tracks in a cone of dR
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      int nTrkSolidCone() const { return photonIDOrThrow().nTrkSolidCone(); }
      /// As above, excluding the core at the center of the cone
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      int nTrkHollowCone() const { return photonIDOrThrow().nTrkHollowCone(); }
      /// Return r9 = e3x3/etotal
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      float r9() const { return photonIDOrThrow().r9(); }
      /// If photon is in ECAL barrel
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isEBPho() const { return photonIDOrThrow().isEBPho(); }
      /// If photon is in ECAL endcap
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isEEPho() const { return photonIDOrThrow().isEEPho(); }
      /// If photon is in EB, and inside the boundaries in super crystals/modules
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isEBGap() const { return photonIDOrThrow().isEBGap(); }
      /// If photon is in EE, and inside the boundaries in supercrystal/D
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isEEGap() const { return photonIDOrThrow().isEEGap(); }
      /// If photon is in boundary between EB and EE
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isEBEEGap() const { return photonIDOrThrow().isEBEEGap(); }
      /// If this is also a GsfElectron
      /// Method from reco::PhotonID, throws exception if there is no photon ID in this pat::Photon
      bool isAlsoElectron() const { return photonIDOrThrow().isAlsoElectron(); }

      // ---- methods for photon isolation ----
      /// Returns the isolation variable for a specifc key (or pseudo-key like CaloIso), or -1.0 if not available
      float isolation(IsolationKeys key) const { 
          if (key >= 0) {
              //if (key >= isolations_.size()) throw cms::Excepton("Missing Data") << "Isolation corresponding to key " << key << " was not stored for this particle.";
              if (size_t(key) >= isolations_.size()) return -1.0;
              return isolations_[key];
          } else switch (key) {
              case CaloIso:  
                  //if (isolations_.size() <= HCalIso) throw cms::Excepton("Missing Data") << "CalIsoo Isolation was not stored for this particle.";
                  if (isolations_.size() <= HCalIso) return -1.0; 
                  return isolations_[ECalIso] + isolations_[HCalIso];
              default:
                  return -1.0;
                  //throw cms::Excepton("Missing Data") << "Isolation corresponding to key " << key << " was not stored for this particle.";
          }
      }
      /// Return the tracker isolation variable that was stored in this object when produced, or -1.0 if there is none
      float trackIso() const { return isolation(TrackerIso); }
      /// Return the sum of ecal and hcal isolation variable that were stored in this object when produced, or -1.0 if at least one is missing
      float caloIso()  const { return isolation(CaloIso); }
      /// Return the ecal isolation variable that was stored in this object when produced, or -1.0 if there is none
      float ecalIso()  const { return isolation(ECalIso); }
      /// Return the hcal isolation variable that was stored in this object when produced, or -1.0 if there is none
      float hcalIso()  const { return isolation(HCalIso); }
      /// Return the user defined isolation variable #index that was stored in this object when produced, or -1.0 if there is none
      float userIso(uint8_t index=0)  const { return isolation(IsolationKeys(UserBaseIso + index)); }
      /// Sets the isolation variable for a specifc key.
      /// Note that you can't set isolation for a pseudo-key like CaloIso
      void setIsolation(IsolationKeys key, float value) {
          if (key >= 0) {
              if (size_t(key) >= isolations_.size()) isolations_.resize(key+1, -1.0);
              isolations_[key] = value;
          } else {
              throw cms::Exception("Illegal Argument") << 
                  "The key for which you're setting isolation does not correspond " <<
                  "to an individual isolation but to the sum of more independent isolations " <<
                  "(e.g. Calo = ECal + HCal), so you can't SET the value, just GET it.\n" <<
                  "Please set up each component independly.\n";
          }
      }
      /// Sets tracker isolation variable
      void setTrackIso(float trackIso) { setIsolation(TrackerIso, trackIso); }
      /// Sets ecal isolation variable
      void setECalIso(float caloIso)   { setIsolation(ECalIso, caloIso);  } 
      /// Sets hcal isolation variable
      void setHCalIso(float caloIso)   { setIsolation(HCalIso, caloIso);  }
      /// Sets user isolation variable #index
      void setUserIso(float value, uint8_t index=0)  { setIsolation(IsolationKeys(UserBaseIso + index), value); }

      // ---- methods for photon isolation deposits ----
      /// Returns the IsoDeposit associated with some key, or a null pointer if it is not available
      const IsoDeposit * isoDeposit(IsolationKeys key) const {
          for (IsoDepositPairs::const_iterator it = isoDeposits_.begin(), ed = isoDeposits_.end(); 
                  it != ed; ++it) 
          {
              if (it->first == key) return & it->second;
          }
          return 0;
      } 
      /// Return the tracker IsoDeposit
      const IsoDeposit * trackerIsoDeposit() const { return isoDeposit(TrackerIso); }
      /// Return the ecal IsoDeposit
      const IsoDeposit * ecalIsoDeposit()    const { return isoDeposit(ECalIso); }
      /// Return the hcal IsoDeposit
      const IsoDeposit * hcalIsoDeposit()    const { return isoDeposit(HCalIso); }
      /// Return a specified user-level IsoDeposit
      const IsoDeposit * userIsoDeposit(uint8_t index=0) const { return isoDeposit(IsolationKeys(UserBaseIso + index)); }
      /// Sets the IsoDeposit associated with some key; if it is already existent, it is overwritten.
      void setIsoDeposit(IsolationKeys key, const IsoDeposit &dep) {
          IsoDepositPairs::iterator it = isoDeposits_.begin(), ed = isoDeposits_.end();
          for (; it != ed; ++it) {
              if (it->first == key) { it->second = dep; return; }
          }
          isoDeposits_.push_back(std::make_pair(key,dep));
      } 
      /// Sets tracker IsoDeposit
      void trackerIsoDeposit(const IsoDeposit &dep) { setIsoDeposit(TrackerIso, dep); }
      /// Sets ecal IsoDeposit
      void ecalIsoDeposit(const IsoDeposit &dep)    { setIsoDeposit(ECalIso, dep); }
      /// Sets hcal IsoDeposit
      void hcalIsoDeposit(const IsoDeposit &dep)    { setIsoDeposit(HCalIso, dep); }
      /// Sets user-level IsoDeposit
      void userIsoDeposit(const IsoDeposit &dep, uint8_t index=0) { setIsoDeposit(IsolationKeys(UserBaseIso + index), dep); }

    protected:

      /// helper method to check photon ID existence
      const reco::PhotonID & photonIDOrThrow() const {
        if (photonID_.empty()) throw cms::Exception("Missing Data") << "This pat::Photon doesn't include a reco::PhotonID.\n";
        return photonID_[0];
      }

    protected:

      // ---- for content embedding ----
      bool embeddedSuperCluster_;
      std::vector<reco::SuperCluster> superCluster_;
      // ---- holder for a reco::PhotonID object ----
      std::vector<reco::PhotonID> photonID_;
      // ---- Isolation and IsoDeposit related datamebers ----
      typedef std::vector<std::pair<IsolationKeys, pat::IsoDeposit> > IsoDepositPairs;
      IsoDepositPairs    isoDeposits_;
      std::vector<float> isolations_;

  };


}

#endif
