//
// $Id: Lepton.h,v 1.16 2009/09/02 13:21:15 veelken Exp $
//

#ifndef DataFormats_PatCandidates_Lepton_h
#define DataFormats_PatCandidates_Lepton_h

/**
  \class    pat::Lepton Lepton.h "DataFormats/PatCandidates/interface/Lepton.h"
  \brief    Analysis-level lepton class

   Lepton implements the analysis-level charged lepton class within the 'pat'
   namespace. It currently provides the link to the generated lepton and
   the isolation information.

   Please post comments and questions to the Physics Tools hypernews:
   https://hypernews.cern.ch/HyperNews/CMS/get/physTools.html

  \author   Steven Lowette, Giovanni Petrucciani, Frederic Ronga
  \version  $Id: Lepton.h,v 1.16 2009/09/02 13:21:15 veelken Exp $
*/

#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"


namespace pat {


  template <class LeptonType>
  class Lepton : public PATObject<LeptonType> {

    public:

      Lepton();
      Lepton(const LeptonType & aLepton);
      Lepton(const edm::RefToBase<LeptonType> & aLeptonRef);
      Lepton(const edm::Ptr<LeptonType> & aLeptonRef);
      virtual ~Lepton();

      virtual Lepton<LeptonType> * clone() const { return new Lepton<LeptonType>(*this); }

      const reco::GenParticle * genLepton() const { return PATObject<LeptonType>::genParticle(); }

      void setGenLepton(const reco::GenParticleRef & gl, bool embed=false) { PATObject<LeptonType>::setGenParticleRef(gl, embed); }

      //============ BEGIN ISOLATION BLOCK =====
      /// Returns the isolation variable for a specifc key (or pseudo-key like CaloIso), or -1.0 if not available
      float userIso(IsolationKeys key) const { 
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

      // ---- specific getters ----
      /// Return the tracker isolation variable that was stored in this object when produced, or -1.0 if there is none
      float trackIso() const { return userIso(TrackerIso); }
      /// Return the sum of ecal and hcal isolation variable that were stored in this object when produced, or -1.0 if at least one is missing
      float caloIso()  const { return userIso(CaloIso); }
      /// Return the ecal isolation variable that was stored in this object when produced, or -1.0 if there is none
      float ecalIso()  const { return userIso(ECalIso); }
      /// Return the hcal isolation variable that was stored in this object when produced, or -1.0 if there is none
      float hcalIso()  const { return userIso(HCalIso); }

      ///PARTICLE FLOW ISOLATION
      ///Return the isolation calculated with all the PFCandidates
      float particleIso() const { return userIso(ParticleIso); }
      ///Return the isolation calculated with only the charged hadron PFCandidates
      float chargedHadronIso() const { return userIso(ChargedHadronIso); }
      ///Return the isolation calculated with only the neutral hadron PFCandidates      
      float neutralHadronIso() const { return userIso(NeutralHadronIso); }	
      ///Return the isolation calculated with only the gamma PFCandidates  
      float photonIso() const { return userIso(PhotonIso); }	

      /// Return the user defined isolation variable #index that was stored in this object when produced, or -1.0 if there is none
      float userIso(uint8_t index=0)  const { return userIso(IsolationKeys(UserBaseIso + index)); }

      // ---- specific setters ----
      /// Sets tracker isolation variable
      void setTrackIso(float trackIso) { setIsolation(TrackerIso, trackIso); }
      /// Sets ecal isolation variable
      void setECalIso(float caloIso)   { setIsolation(ECalIso, caloIso);  } 
      /// Sets hcal isolation variable
      void setHCalIso(float caloIso)   { setIsolation(HCalIso, caloIso);  }
      /// Sets user isolation variable #index
      void setUserIso(float value, uint8_t index=0)  { setIsolation(IsolationKeys(UserBaseIso + index), value); }


      //============ BEGIN ISODEPOSIT BLOCK =====
      /// Returns the IsoDeposit associated with some key, or a null pointer if it is not available
      const IsoDeposit * isoDeposit(IsolationKeys key) const {
          for (IsoDepositPairs::const_iterator it = isoDeposits_.begin(), ed = isoDeposits_.end(); 
                  it != ed; ++it) 
          {
              if (it->first == key) return & it->second;
          }
          return 0;
      } 

      /// Sets the IsoDeposit associated with some key; if it is already existent, it is overwritten.
      void setIsoDeposit(IsolationKeys key, const IsoDeposit &dep) {
          IsoDepositPairs::iterator it = isoDeposits_.begin(), ed = isoDeposits_.end();
          for (; it != ed; ++it) {
              if (it->first == key) { it->second = dep; return; }
          }
          isoDeposits_.push_back(std::make_pair(key,dep));
      } 

      // ---- specific getters ----
      const IsoDeposit * trackerIsoDeposit() const { return isoDeposit(TrackerIso); }
      const IsoDeposit * ecalIsoDeposit()    const { return isoDeposit(ECalIso); }
      const IsoDeposit * hcalIsoDeposit()    const { return isoDeposit(HCalIso); }
      const IsoDeposit * userIsoDeposit(uint8_t index=0) const { return isoDeposit(IsolationKeys(UserBaseIso + index)); }

      // ---- specific setters ----
      void trackerIsoDeposit(const IsoDeposit &dep) { setIsoDeposit(TrackerIso, dep); }
      void ecalIsoDeposit(const IsoDeposit &dep)    { setIsoDeposit(ECalIso, dep); }
      void hcalIsoDeposit(const IsoDeposit &dep)    { setIsoDeposit(HCalIso, dep); }
      void userIsoDeposit(const IsoDeposit &dep, uint8_t index=0) { setIsoDeposit(IsolationKeys(UserBaseIso + index), dep); }


    protected:
      // --- Isolation and IsoDeposit related datamebers ---
      typedef std::vector<std::pair<IsolationKeys, pat::IsoDeposit> > IsoDepositPairs;
      IsoDepositPairs    isoDeposits_;
      std::vector<float> isolations_;
  };


  /// default constructor
  template <class LeptonType>
  Lepton<LeptonType>::Lepton() :
    PATObject<LeptonType>(LeptonType()) {
    // no common constructor, so initialize the candidate manually
    this->setCharge(0);
    this->setP4(reco::Particle::LorentzVector(0, 0, 0, 0));
    this->setVertex(reco::Particle::Point(0, 0, 0));
  }


  /// constructor from LeptonType
  template <class LeptonType>
  Lepton<LeptonType>::Lepton(const LeptonType & aLepton) :
    PATObject<LeptonType>(aLepton) {
  }


  /// constructor from ref to LeptonType
  template <class LeptonType>
  Lepton<LeptonType>::Lepton(const edm::RefToBase<LeptonType> & aLeptonRef) :
    PATObject<LeptonType>(aLeptonRef) {
  }


  /// constructor from ref to LeptonType
  template <class LeptonType>
  Lepton<LeptonType>::Lepton(const edm::Ptr<LeptonType> & aLeptonRef) :
    PATObject<LeptonType>(aLeptonRef) {
  }


  /// destructor
  template <class LeptonType>
  Lepton<LeptonType>::~Lepton() {
  }
}

#endif
