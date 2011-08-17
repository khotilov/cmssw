
#include "HLTrigger/JetMET/interface/HLTHtMhtFilter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "DataFormats/METReco/interface/MET.h"


HLTHtMhtFilter::HLTHtMhtFilter(const edm::ParameterSet & iConfig) {

  saveTags_  = iConfig.getParameter<bool>("saveTags");
  htLabels_  = iConfig.getParameter<std::vector<edm::InputTag> >("htLabels");
  mhtLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("mhtLabels");
  minHt_     = iConfig.getParameter<std::vector<double> >("minHt");
  minMht_    = iConfig.getParameter<std::vector<double> >("minMht");
  minMeff_   = iConfig.getParameter<std::vector<double> >("minMeff");
  meffSlope_ = iConfig.getParameter<std::vector<double> >("meffSlope");
  orLabels_  = iConfig.getParameter<std::vector<edm::InputTag> >("toBeORdLabels");

  nOrs_ = htLabels_.size(); // number of settings to .OR.
  if (!( htLabels_.size() == mhtLabels_.size() &&
         htLabels_.size() == minHt_.size() &&
         htLabels_.size() == minMht_.size() &&
         htLabels_.size() == minMeff_.size() &&
         htLabels_.size() == meffSlope_.size() ) ||
	 htLabels_.size() == 0 ) {
    nOrs_ = (mhtLabels_.size() < nOrs_ ? mhtLabels_.size() : nOrs_);
    nOrs_ = (minHt_.size()     < nOrs_ ? minHt_.size()     : nOrs_);
    nOrs_ = (minMht_.size()    < nOrs_ ? minMht_.size()    : nOrs_);
    nOrs_ = (minMeff_.size()   < nOrs_ ? minMeff_.size()   : nOrs_);
    nOrs_ = (meffSlope_.size() < nOrs_ ? meffSlope_.size() : nOrs_);
    edm::LogError("HLTHtMhtFilter") << "inconsistent module configuration!";
  }

  moduleLabel_ = iConfig.getParameter<std::string>("@module_label");
  produces<reco::METCollection>();
  produces<trigger::TriggerFilterObjectWithRefs>();

}


HLTHtMhtFilter::~HLTHtMhtFilter() {
}


void HLTHtMhtFilter::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
  std::vector<edm::InputTag> tmp1; tmp1.push_back(edm::InputTag(""));
  std::vector<double>        tmp2; tmp2.push_back(0.);
  edm::ParameterSetDescription desc;
  desc.add<bool>("saveTags",false);
  desc.add<std::vector<edm::InputTag> >("htLabels",  tmp1);
  desc.add<std::vector<edm::InputTag> >("mhtLabels", tmp1);
  desc.add<std::vector<double> >("minHt",     tmp2);
  desc.add<std::vector<double> >("minMht",    tmp2);
  desc.add<std::vector<double> >("minMeff",   tmp2);
  desc.add<std::vector<double> >("meffSlope", tmp2);
  descriptions.add("hltHtMhtFilter", desc);
}


bool HLTHtMhtFilter::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  // the filter objects to be stored
  std::auto_ptr<reco::METCollection> metobject(new reco::METCollection());
  // the references to the filter objects
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> filterobject (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if (saveTags_) filterobject->addCollectionTag(moduleLabel_);

  // check first the extra input bools
  for (unsigned int i = 0; i < orLabels_.size(); ++i) {
    edm::Handle<bool> orhandle;
    iEvent.getByLabel(orLabels_[i], orhandle);
    if (*orhandle) {
      // put empty objects in the event (so we know this was accepted on an extra bool input)
      iEvent.put(metobject);
      iEvent.put(filterobject);
      // and accept the event
      return true;
    }
  }

  bool accept = false;

  // take the .OR. of all sets of constraints
  for (unsigned int i = 0; i < nOrs_; ++i) {

    // read in the HT and mHT
    edm::Handle<std::vector<reco::MET> > hht;
    iEvent.getByLabel(htLabels_[i], hht);
    double ht = (*hht)[0].sumEt();
    edm::Handle<std::vector<reco::MET> > hmht;
    iEvent.getByLabel(mhtLabels_[i], hmht);
    double mht = (*hmht)[0].pt();

    // check if the event passes this cut set
    accept = accept || (ht > minHt_[i] && mht > minMht_[i] && sqrt(mht + meffSlope_[i]*ht) > minMeff_[i]);
    // in principle we could break if accepted, but in order to save
    // for offline analysis all possible decisions we keep looping here
    // in term of timing this will not matter much; typically 1 or 2 cut-sets
    // will be checked only

    // store the object that was cut on and the ref to it
    metobject->push_back(reco::MET(ht, (*hmht)[0].p4(), reco::MET::Point()));
    edm::Ref<reco::METCollection> metref(iEvent.getRefBeforePut<reco::METCollection>(), i); // point to i'th object
    filterobject->addObject(trigger::TriggerMHT, metref); // save as an MHT

  }

  // put filter object into the Event
  iEvent.put(metobject);
  iEvent.put(filterobject);

  return accept;

}
