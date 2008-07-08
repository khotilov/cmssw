#ifndef HLTcore_TriggerSummaryProducerAOD_h
#define HLTcore_TriggerSummaryProducerAOD_h

/** \class TriggerSummaryProducerAOD
 *
 *  
 *  This class is an EDProducer making the HLT summary object for AOD
 *
 *  $Date: 2008/05/19 13:16:46 $
 *  $Revision: 1.8 $
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

#include<string>
#include<vector>

//
// class declaration
//

struct TriggerInputTagComparison {
  bool operator() (const edm::InputTag& lhs, const edm::InputTag& rhs) const {
    return lhs.encode()<rhs.encode();
  }
};

typedef std::set<edm::InputTag,TriggerInputTagComparison> InputTagSet;

class TriggerSummaryProducerAOD : public edm::EDProducer {
  
 public:
  explicit TriggerSummaryProducerAOD(const edm::ParameterSet&);
  ~TriggerSummaryProducerAOD();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob();


  // additional

  template <typename C>
  void fillTriggerObjects(const edm::Event& );

  template <typename C>
  void fillFilterObjects(const edm::Event&, const edm::InputTag& tag, const trigger::Vids &, const std::vector<edm::Ref<C> >&);

 private:
  /// process name
  std::string pn_;

  /// selector for getMany methods
  edm::ProcessNameSelector selector_;

  /// the pointer to the current TriggerNamesService
  edm::service::TriggerNamesService* tns_;

  /// lists of L3 collection labels
  std::vector<edm::InputTag> collectionTags_;
  ///
  InputTagSet                collectionTagsEvent_;
  InputTagSet                collectionTagsGlobal_;
  /// list of L3 filter labels
  InputTagSet                filterTagsEvent_;
  InputTagSet                filterTagsGlobal_;

  /// trigger object collection
  trigger::TriggerObjectCollection toc_;
  std::vector<edm::InputTag> tags_;
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

};
#endif
