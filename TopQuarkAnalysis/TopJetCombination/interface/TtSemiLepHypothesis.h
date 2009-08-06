#ifndef TtSemiLepHypothesis_h
#define TtSemiLepHypothesis_h

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"

class TtSemiLepHypothesis : public edm::EDProducer {

 public:

  explicit TtSemiLepHypothesis(const edm::ParameterSet&);
  ~TtSemiLepHypothesis();

 protected:
  
  /// produce the event hypothesis as CompositeCandidate and Key
  virtual void produce(edm::Event&, const edm::EventSetup&);
  /// reset candidate pointers before hypo build process
  void resetCandidates();
  /// use one object in a collection to set a ShallowClonePtrCandidate
  template <typename C>
  void setCandidate(const edm::Handle<C>& handle, const int& idx, reco::ShallowClonePtrCandidate*& clone);
  /// return key
  int key() const { return key_; };
  /// return event hypothesis
  reco::CompositeCandidate hypo();
  /// check if index is in valid range of selected jets
  bool isValid(const int& idx, const edm::Handle<std::vector<pat::Jet> >& jets){ return (0<=idx && idx<(int)jets->size()); };
  /// determine lepton type of reco candidate and return a corresponding WDecay::LeptonType; the type is kNone if it is whether a muon nor an electron 
  WDecay::LeptonType leptonType(const reco::RecoCandidate* cand);

  // -----------------------------------------
  // implemet the following two functions
  // for a concrete event hypothesis
  // -----------------------------------------

  /// build the event hypothesis key
  virtual void buildKey() = 0;
  /// build event hypothesis from the reco objects of a semi-leptonic event 
  virtual void buildHypo(edm::Event& event,
			 const edm::Handle<edm::View<reco::RecoCandidate> >& lepton, 
			 const edm::Handle<std::vector<pat::MET> >& neutrino, 
			 const edm::Handle<std::vector<pat::Jet> >& jets, 
			 std::vector<int>& jetPartonAssociation,
			 const unsigned int iComb) = 0;

 protected:

  bool getMatch_;

  edm::InputTag jets_;
  edm::InputTag leps_;
  edm::InputTag mets_;
  edm::InputTag match_;

  int key_;

  reco::ShallowClonePtrCandidate *lightQ_;
  reco::ShallowClonePtrCandidate *lightQBar_;
  reco::ShallowClonePtrCandidate *hadronicB_;
  reco::ShallowClonePtrCandidate *leptonicB_;
  reco::ShallowClonePtrCandidate *neutrino_;
  reco::ShallowClonePtrCandidate *lepton_;
};

// has to be placed in the header since otherwise the function template
// would cause unresolved references in classes derived from this base class
template<typename C>
void
TtSemiLepHypothesis::setCandidate(const edm::Handle<C>& handle, const int& idx, reco::ShallowClonePtrCandidate* &clone) {
  typedef typename C::value_type O;
  edm::Ptr<O> ptr = edm::Ptr<O>(handle, idx);
  clone = new reco::ShallowClonePtrCandidate( ptr, ptr->charge(), ptr->p4(), ptr->vertex() );
}

#endif
