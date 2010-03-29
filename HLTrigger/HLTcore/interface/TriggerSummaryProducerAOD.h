#ifndef HLTcore_TriggerSummaryProducerAOD_h
#define HLTcore_TriggerSummaryProducerAOD_h

/** \class TriggerSummaryProducerAOD
 *
 *  
 *  This class is an EDProducer making the HLT summary object for AOD
 *
 *  $Date: 2010/02/15 17:43:23 $
 *  $Revision: 1.11 $
 *
 *  \author Martin Grunewald
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Selector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"


#include<string>
#include<vector>

//
// class declaration
//

class TriggerSummaryProducerAOD : public edm::EDProducer {
  
  typedef std::set<std::string> InputStringSet;

 public:
  explicit TriggerSummaryProducerAOD(const edm::ParameterSet&);
  ~TriggerSummaryProducerAOD();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  void tokenizeTag(const std::string& tag, std::string& label, std::string& instance, std::string& process) const;

  // additional

  template <typename C>
  void fillTriggerObjectCollections(const edm::Event& );

  template <typename T>
  void fillTriggerObject(const T& );
  void fillTriggerObject(const l1extra::L1HFRings& );
  void fillTriggerObject(const l1extra::L1EtMissParticle& );
  void fillTriggerObject(const reco::CaloMET& );
  void fillTriggerObject(const reco::MET& );

  template <typename C>
    void fillFilterObjectMembers(const edm::Event&, const std::string& tag, const trigger::Vids &, const std::vector<edm::Ref<C> >&);

  template <typename C>
  void fillFilterObjectMember(const int&, const int&, const edm::Ref<C>&);
  void fillFilterObjectMember(const int&, const int&, const edm::Ref<l1extra::L1HFRingsCollection>&);
  void fillFilterObjectMember(const int&, const int&, const edm::Ref<l1extra::L1EtMissParticleCollection>&);
  void fillFilterObjectMember(const int&, const int&, const edm::Ref<reco::CaloMETCollection>&);
  void fillFilterObjectMember(const int&, const int&, const edm::Ref<reco::METCollection>&);

 private:
  /// process name
  std::string pn_;

  /// selector for getMany methods
  edm::ProcessNameSelector selector_;

  /// the pointer to the current TriggerNamesService
  edm::service::TriggerNamesService* tns_;

  /// list of L3 collection labels
  InputStringSet                collectionTagsEvent_;
  InputStringSet                collectionTagsGlobal_;
  /// list of L3 filter labels
  InputStringSet                filterTagsEvent_;
  InputStringSet                filterTagsGlobal_;

  /// trigger object collection
  trigger::TriggerObjectCollection toc_;
  std::vector<std::string> tags_;
  /// global map for indices into toc_: offset per input L3 collection
  std::map<edm::ProductID,int> offset_;

  /// handles to the filter objects
  std::vector<edm::Handle<trigger::TriggerFilterObjectWithRefs> > fobs_;
  /// keys
  trigger::Keys keys_;
  /// ids
  trigger::Vids ids_;

  /// packing decision
  std::vector<bool> maskFilters_;

  /// tokenized version using InputTag
  std::vector<edm::InputTag> collectionTokensEvent_; 

};
#endif
