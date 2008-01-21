#include "HLTrigger/special/interface/HLTHcalPhiSymFilter.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

HLTHcalPhiSymFilter::HLTHcalPhiSymFilter(const edm::ParameterSet& iConfig)
{
  HBHEHits_ = iConfig.getParameter< edm::InputTag > ("HBHEHitCollection");
  HOHits_ = iConfig.getParameter< edm::InputTag > ("HOHitCollection");
  HFHits_=iConfig.getParameter<edm::InputTag>("HFHitCollection");
  phiSymHBHEHits_ = iConfig.getParameter< std::string > ("phiSymHBHEHitCollection");
  phiSymHOHits_ = iConfig.getParameter< std::string > ("phiSymHOHitCollection");
  phiSymHFHits_ =iConfig.getParameter< std::string > ("phiSymHFHitCollection");

  eCut_HB_ = iConfig.getParameter< double > ("eCut_HB");
  eCut_HE_ = iConfig.getParameter< double > ("eCut_HE");
  eCut_HO_ = iConfig.getParameter<double>("eCut_HO");
  eCut_HF_ = iConfig.getParameter<double>("eCut_HF");

  //register your products
  produces< HBHERecHitCollection >(phiSymHBHEHits_);
  produces< HORecHitCollection >(phiSymHOHits_);
  produces< HFRecHitCollection >(phiSymHFHits_);	

  produces<trigger::TriggerFilterObjectWithRefs>();
}


HLTHcalPhiSymFilter::~HLTHcalPhiSymFilter()
{
 

}


// ------------ method called to produce the data  ------------
bool
HLTHcalPhiSymFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  edm::Handle<HORecHitCollection> HORecHitsH;
  edm::Handle<HFRecHitCollection> HFRecHitsH;
  
  iEvent.getByLabel(HBHEHits_,HBHERecHitsH);
  iEvent.getByLabel(HOHits_,HORecHitsH);
  iEvent.getByLabel(HFHits_,HFRecHitsH);
 
  //Create empty output collections
  std::auto_ptr< HBHERecHitCollection > phiSymHBHERecHitCollection( new HBHERecHitCollection );
  std::auto_ptr< HORecHitCollection > phiSymHORecHitCollection( new HORecHitCollection );
  std::auto_ptr< HFRecHitCollection > phiSymHFRecHitCollection( new HFRecHitCollection );

  // The Filter object. We don't really need to put anything into it, but we 
  // write an empty one for consistency

  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> filterproduct (new trigger::TriggerFilterObjectWithRefs(path(),module()));

  //Select interesting HBHERecHits 
  for (HBHERecHitCollection::const_iterator it=HBHERecHitsH->begin(); it!=HBHERecHitsH->end(); it++) {
    if (it->energy()>eCut_HB_&&it->id().subdet()==1) {
        phiSymHBHERecHitCollection->push_back(*it);
    }
    if (it->energy()>eCut_HE_&&it->id().subdet()==2) {
        phiSymHBHERecHitCollection->push_back(*it);
    }

  }
  
  //Select interesting HORecHits
  for (HORecHitCollection::const_iterator it=HORecHitsH->begin(); it!=HORecHitsH->end(); it++) {
    if (it->energy()>eCut_HO_&&it->id().subdet()==3) {
      phiSymHORecHitCollection->push_back(*it);
    }
  }

  //Select interesting HFRecHits
  for (HFRecHitCollection::const_iterator it=HFRecHitsH->begin(); it!=HFRecHitsH->end(); it++) {
    if (it->energy()>eCut_HF_&&it->id().subdet()==4) {
      phiSymHFRecHitCollection->push_back(*it);
    }
  }

  if ((!phiSymHBHERecHitCollection->size() ) || (!phiSymHORecHitCollection->size()) || (!phiSymHFRecHitCollection->size())) return false;

  //Put selected information in the event
  if (phiSymHBHERecHitCollection->size()>0) iEvent.put( phiSymHBHERecHitCollection, phiSymHBHEHits_);
  if (phiSymHORecHitCollection->size()>0) iEvent.put( phiSymHORecHitCollection, phiSymHOHits_);
  if (phiSymHFRecHitCollection->size()>0) iEvent.put( phiSymHFRecHitCollection, phiSymHFHits_);
  
  iEvent.put(filterproduct);
 
  return true;
}
