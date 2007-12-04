#ifndef HLTReco_TriggerFilterObjectWithRefs_h
#define HLTReco_TriggerFilterObjectWithRefs_h

/** \class trigger::TriggerFilterObjectWithRefs
 *
 *  If HLT cuts of intermediate or final HLT filters are satisfied,
 *  instances of this class hold the combination of reconstructed
 *  physics objects (e/gamma/mu/jet/MMet...) satisfying the cuts.
 *
 *  This implementation is not completely space-efficient as some
 *  physics object containers may stay empty. However, the big
 *  advantage is that the solution is generic, i.e., works for all
 *  possible HLT filters. Hence we accept the reasonably small
 *  overhead of empty containers.
 *
 *  $Date: 2007/12/04 08:35:53 $
 *  $Revision: 1.5 $
 *
 *  \author Martin Grunewald
 *
 */

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Provenance/interface/ProductID.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include <vector>

namespace trigger
{

  /// Transient book-keeping EDProduct filled by HLTFilter module to
  /// record physics objects firing the filter (never persistet in
  /// production; same functionality but different implementation
  /// compared to the old HLT data model's HLTFilterObjectWithRefs
  /// class)
  class TriggerFilterObjectWithRefs {

  /// data members
  private:
    std::vector<reco::RecoEcalCandidateRef> photons_;
    std::vector<reco::ElectronRef> electrons_;
    std::vector<reco::RecoChargedCandidateRef> muons_;
    std::vector<reco::CaloJetRef> jets_;
    std::vector<reco::CompositeCandidateRef> composites_;
    std::vector<reco::CaloMETRef> mets_;
    std::vector<reco::METRef> hts_;
    std::vector<XRef> others_;
    
  /// methods
  public:
    /// constructors
    TriggerFilterObjectWithRefs():
      photons_(), electrons_(), muons_(), jets_(), composites_(), mets_(), hts_(), others_() { }

    /// setters for L3 collections
    void addPhoton(const reco::RecoEcalCandidateRef& photon) {photons_.push_back(photon);}
    void addElectron(const reco::ElectronRef& electron) {electrons_.push_back(electron);}
    void addMuon(const reco::RecoChargedCandidateRef& muon) {muons_.push_back(muon);}
    void addJet(const reco::CaloJetRef& jet) {jets_.push_back(jet);}
    void addComposite(const reco::CompositeCandidateRef& composite) {composites_.push_back(composite);}
    void addMET(const reco::CaloMETRef& met) {mets_.push_back(met);}
    void addHT (const reco::METRef& ht) {hts_.push_back(ht);}

    /// setter for non-L3 collections
    ///   typesafe as original collection is kept and user will match based
    ///   on ProductID - before using key as index into original collection
    template <typename C>
    void addOther(const edm::Ref<C>& ref) {
      addOther(ref.id(),ref.key());
    }
    void addOther(edm::ProductID id, size_type key) {
      others_.push_back(XRef(id,key));
    }

    /// getters
    const std::vector<reco::RecoEcalCandidateRef>& getPhotons() const {return photons_;}
    const std::vector<reco::ElectronRef>& getElectrons() const {return electrons_;}
    const std::vector<reco::RecoChargedCandidateRef>& getMuons() const {return muons_;}
    const std::vector<reco::CaloJetRef>& getJets() const {return jets_;}
    const std::vector<reco::CompositeCandidateRef>& getComposites() const {return composites_;}
    const std::vector<reco::CaloMETRef>& getMETs() const {return mets_;}
    const std::vector<reco::METRef>& getHTs() const {return hts_;}
    const std::vector<XRef>& getOthers() const {return others_;}

    template <typename C>
    void getOthers(const edm::Handle<C>& handle, std::vector<edm::Ref<C> >& vref) const {
      vref.resize(0);
      const size_type n(others_.size());
      for (size_type i=0; i!=n; ++i) {
	if (handle.id()==others_[i].first) {
	  vref.push_back(edm::Ref<C>(handle,others_[i].second));
	}
      }
      return;
    }

  };

}

#endif
