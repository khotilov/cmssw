/** \class HLTJetCollectionsFilter
 *
 *
 *  \author Leonard Apanasevich
 *
 */

#include "HLTrigger/JetMET/interface/HLTJetCollectionsFilter.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include<typeinfo>
#include<string>
//
// constructors and destructor
//
template <typename jetType>
HLTJetCollectionsFilter<jetType>::HLTJetCollectionsFilter(const edm::ParameterSet& iConfig):
  HLTFilter(iConfig),
  inputTag_(iConfig.getParameter< edm::InputTag > ("inputTag")),
  originalTag_(iConfig.getParameter< edm::InputTag > ("originalTag")),
  minJetPt_(iConfig.getParameter<double> ("MinJetPt")),
  maxAbsJetEta_(iConfig.getParameter<double> ("MaxAbsJetEta")),
  minNJets_(iConfig.getParameter<unsigned int> ("MinNJets")),
  triggerType_(iConfig.getParameter<int> ("triggerType"))
{
}

template <typename jetType>
HLTJetCollectionsFilter<jetType>::~HLTJetCollectionsFilter(){}

template <typename jetType>
void
HLTJetCollectionsFilter<jetType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  makeHLTFilterDescription(desc);
  desc.add<edm::InputTag>("inputTag",edm::InputTag("hltIterativeCone5CaloJets"));
  desc.add<edm::InputTag>("originalTag",edm::InputTag("hltIterativeCone5CaloJets"));
  desc.add<double>("MinJetPt",30.0);
  desc.add<double>("MaxAbsJetEta",2.6);
  desc.add<unsigned int>("MinNJets",1);
  desc.add<int>("triggerType",0);
  descriptions.add(std::string("hlt")+std::string(typeid(HLTJetCollectionsFilter<jetType>).name()),desc);
}

// ------------ method called to produce the data  ------------
template <typename jetType>
bool
HLTJetCollectionsFilter<jetType>::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  typedef vector<RefVector<vector<jetType>,jetType,refhelper::FindUsingAdvance<vector<jetType>,jetType> > > JetCollectionVector;
  typedef vector<jetType> JetCollection;
  typedef edm::RefVector<JetCollection> JetRefVector;
  typedef edm::Ref<JetCollection> JetRef;

  // The filter object
  if (saveTags()) filterproduct.addCollectionTag(originalTag_);

  Handle < JetCollectionVector > theJetCollectionsHandle;
  iEvent.getByLabel(inputTag_, theJetCollectionsHandle);
  const JetCollectionVector & theJetCollections = *theJetCollectionsHandle;

  // filter decision
  bool accept(false);
  std::vector < JetRef > goodJetRefs;

  for (unsigned int collection = 0; collection < theJetCollections.size(); ++collection) {
    unsigned int numberOfGoodJets(0);
    const JetRefVector & refVector = theJetCollections[collection];

    if (refVector.size() < minNJets_) continue;

    //empty the good jets collection
    goodJetRefs.clear();

    typename JetRefVector::const_iterator jet(refVector.begin());
    for (; jet != refVector.end(); jet++) {
      JetRef jetRef(*jet);
      if (jetRef->pt() >= minJetPt_ && fabs(jetRef->eta()) <= maxAbsJetEta_){
    	  numberOfGoodJets++;
    	  goodJetRefs.push_back(jetRef);
      }
    }
    if (numberOfGoodJets >= minNJets_) {
      accept = true;
      break;
    }
  }

  //fill the filter object
  for (unsigned int refIndex = 0; refIndex < goodJetRefs.size(); ++refIndex) {
    filterproduct.addObject(triggerType_,goodJetRefs.at(refIndex));
  }

  return accept;
}
