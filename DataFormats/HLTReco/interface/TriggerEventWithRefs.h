#ifndef HLTReco_TriggerEventWithRefs_h
#define HLTReco_TriggerEventWithRefs_h

/** \class trigger::TriggerEventWithRefs
 *
 *  The single EDProduct to be saved for events (RAW case)
 *  describing the details of the (HLT) trigger table
 *
 *  $Date: 2007/11/26 16:55:56 $
 *  $Revision: 1.1 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TauReco/interface/HLTTauFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include <string>
#include <vector>

namespace trigger
{

  typedef uint16_t size_type;

  /// The single EDProduct to be saved in addition for each event
  /// - but only in the RAW case
  class TriggerEventWithRefs {

  private:

    /// Helper class: trigger objects firing a single filter
    class TriggerFilterObject {
    public:
      /// label of filter
      std::string filterLabel_;
      /// 1-after-end (std C++) indices into linearised vector of Refs
      /// (-> first start index is always 0)
      size_type photons_;
      size_type electrons_;
      size_type muons_;
      size_type taus_;
      size_type jets_;
      size_type mets_;
      size_type hts_;
      /// constructor
      TriggerFilterObject() :
	filterLabel_(),
	photons_(0), electrons_(0), muons_(0), taus_(0), jets_(0), mets_(0), hts_(0) { }
      TriggerFilterObject(const std::string& filterLabel,
        size_type np, size_type ne, size_type nm, size_type nt, size_type nj, size_type nM, size_type nH) :
	filterLabel_(filterLabel),
	photons_(np), electrons_(ne), muons_(nm), taus_(nt), jets_(nj), mets_(nM), hts_(nH) { }
    };

  /// data members
  private:
    /// the filters recorded here
    std::vector<TriggerFilterObject> filterObjects_;
    /// non-owning pointers into collections (linearised)
    std::vector<reco::RecoEcalCandidateRef> photons_;
    std::vector<reco::ElectronRef> electrons_;
    std::vector<reco::RecoChargedCandidateRef> muons_;
    std::vector<reco::CaloJetRef> taus_;
    std::vector<reco::CaloJetRef> jets_;
    std::vector<reco::CaloMETRef> mets_;
    std::vector<reco::METRef> hts_;
   
  /// methods
  public:
    /// constructors
    TriggerEventWithRefs(): filterObjects_(),
      photons_(), electrons_(), muons_(), taus_(), jets_(), mets_(), hts_() { }

    /// setters - to build EDProduct
    void addFilterObject(const std::string filterLabel, const TriggerFilterObjectWithRefs& tfowr) {
      photons_.insert(photons_.end(),tfowr.getPhotons().begin(),tfowr.getPhotons().end());
      electrons_.insert(electrons_.end(),tfowr.getElectrons().begin(),tfowr.getElectrons().end());
      muons_.insert(muons_.end(),tfowr.getMuons().begin(),tfowr.getMuons().end());
      taus_.insert(taus_.end(),tfowr.getTaus().begin(),tfowr.getTaus().end());
      jets_.insert(jets_.end(),tfowr.getJets().begin(),tfowr.getJets().end());
      mets_.insert(mets_.end(),tfowr.getMETs().begin(),tfowr.getMETs().end());
      hts_.insert(hts_.end(),tfowr.getHTs().begin(),tfowr.getHTs().end());
      filterObjects_.push_back(
        TriggerFilterObject(filterLabel, 
			    photons_.size(), electrons_.size(), 
			    muons_.size(), taus_.size(), jets_.size(),
			    mets_.size(), hts_.size()
			   )
	);
    }

    /// getters - for user access
    size_type numFilters() const {return filterObjects_.size();}

    const std::string& getFilterLabel(size_type index) const {return filterObjects_.at(index).filterLabel_;}

    size_type find(const std::string& filterLabel) const {
      const size_type n(numFilters());
      for (size_type i=0; i!=n; ++i) {
	if (filterLabel==filterObjects_[i].filterLabel_) {return i;}
      }
      return n;
    }

    size_type numPhotons(size_type index) const {
      return filterObjects_.at(index).photons_ - (index==0? 0 : filterObjects_.at(index-1).photons_);
    }
    size_type numPhotons(const std::string& filterLabel) const {
      return numPhotons(find(filterLabel));
    }

    size_type numElectrons(size_type index) const {
      return filterObjects_.at(index).electrons_ - (index==0? 0 : filterObjects_.at(index-1).electrons_);
    }
    size_type numElectrons(const std::string& filterLabel) const {
      return numElectrons(find(filterLabel));
    }

    size_type numMuons(size_type index) const {
      return filterObjects_.at(index).muons_ - (index==0? 0 : filterObjects_.at(index-1).muons_);
    }
    size_type numMuons(const std::string& filterLabel) const {
      return numMuons(find(filterLabel));
    }

    size_type numTaus(size_type index) const {
      return filterObjects_.at(index).taus_ - (index==0? 0 : filterObjects_.at(index-1).taus_);
    }
    size_type numTaus(const std::string& filterLabel) const {
      return numTaus(find(filterLabel));
    }

    size_type numJets(size_type index) const {
      return filterObjects_.at(index).jets_ - (index==0? 0 : filterObjects_.at(index-1).jets_);
    }
    size_type numJets(const std::string& filterLabel) const {
      return numJets(find(filterLabel));
    }

    size_type numMETs(size_type index) const {
      return filterObjects_.at(index).mets_ - (index==0? 0 : filterObjects_.at(index-1).mets_);
    }
    size_type numMETs(const std::string& filterLabel) const {
      return numMETs(find(filterLabel));
    }

    size_type numHTs(size_type index) const {
      return filterObjects_.at(index).hts_ - (index==0? 0 : filterObjects_.at(index-1).hts_);
    }
    size_type numHTs(const std::string& filterLabel) const {
      return numHTs(find(filterLabel));
    }

    /// iterators

    std::vector<reco::RecoEcalCandidateRef>::const_iterator
      photons_begin(size_type index) const
    { return photons_.begin() + 
      (index==0? 0 : filterObjects_.at(index-1).photons_);
    }
    std::vector<reco::RecoEcalCandidateRef>::const_iterator
      photons_end(size_type index) const
    { return photons_.begin() + filterObjects_.at(index).photons_; }


    std::vector<reco::ElectronRef>::const_iterator
      electrons_begin(size_type index) const
    { return electrons_.begin() + 
      (index==0? 0 : filterObjects_[index-1].electrons_);
    }
    std::vector<reco::ElectronRef>::const_iterator
      electrons_end(size_type index) const
    { return electrons_.begin() + filterObjects_[index].electrons_; }


    std::vector<reco::RecoChargedCandidateRef>::const_iterator
      muons_begin(size_type index) const
    { return muons_.begin() + 
      (index==0? 0 : filterObjects_[index-1].muons_);
    }
    std::vector<reco::RecoChargedCandidateRef>::const_iterator
      muons_end(size_type index) const
    { return muons_.begin() + filterObjects_[index].muons_; }


    std::vector<reco::CaloJetRef>::const_iterator
      taus_begin(size_type index) const
    { return taus_.begin() + 
      (index==0? 0 : filterObjects_[index-1].taus_);
    }
    std::vector<reco::CaloJetRef>::const_iterator
      taus_end(size_type index) const
    { return taus_.begin() + filterObjects_[index].taus_; }


    std::vector<reco::CaloJetRef>::const_iterator
      jets_begin(size_type index) const
    { return jets_.begin() + 
      (index==0? 0 : filterObjects_[index-1].jets_);
    }
    std::vector<reco::CaloJetRef>::const_iterator
      jets_end(size_type index) const
    { return jets_.begin() + filterObjects_[index].jets_; }


    std::vector<reco::CaloMETRef>::const_iterator
      mets_begin(size_type index) const
    { return mets_.begin() + 
      (index==0? 0 : filterObjects_[index-1].mets_);
    }
    std::vector<reco::CaloMETRef>::const_iterator
      mets_end(size_type index) const
    { return mets_.begin() + filterObjects_[index].mets_; }


    std::vector<reco::METRef>::const_iterator
      hts_begin(size_type index) const
    { return hts_.begin() + 
      (index==0? 0 : filterObjects_[index-1].hts_);
    }
    std::vector<reco::METRef>::const_iterator
      hts_end(size_type index) const
    { return hts_.begin() + filterObjects_[index].hts_; }

  };

}

#endif
