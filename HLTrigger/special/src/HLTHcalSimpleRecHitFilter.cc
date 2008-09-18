// -*- C++ -*-
//
// Package:    HLTHcalSimpleRecHitFilter
// Class:      HLTHcalSimpleRecHitFilter
// 
/**\class HLTHcalSimpleRecHitFilter HLTHcalSimpleRecHitFilter.cc Work/HLTHcalSimpleRecHitFilter/src/HLTHcalSimpleRecHitFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Bryan DAHMES
//         Created:  Wed Sep 19 16:21:29 CEST 2007
// $Id: HLTHcalSimpleRecHitFilter.cc,v 1.2 2008/09/17 12:04:35 mzanetti Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

//
// class declaration
//

class HLTHcalSimpleRecHitFilter : public HLTFilter {
public:
    explicit HLTHcalSimpleRecHitFilter(const edm::ParameterSet&);
    ~HLTHcalSimpleRecHitFilter();
    
private:
    virtual void beginJob(const edm::EventSetup&) ;
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    edm::InputTag HcalRecHitCollection_;
    double threshold_;
    std::vector<int> maskedList_;
    
};

//
// constructors and destructor
//
HLTHcalSimpleRecHitFilter::HLTHcalSimpleRecHitFilter(const edm::ParameterSet& iConfig)
{
    //now do what ever initialization is needed
    threshold_     = iConfig.getUntrackedParameter<double>("threshold", 0);
    maskedList_    = iConfig.getUntrackedParameter<std::vector<int> >("maskedChannels", maskedList_); //this is using the ashed index
    HcalRecHitCollection_ = iConfig.getParameter<edm::InputTag>("HFRecHitCollection");
    
}


HLTHcalSimpleRecHitFilter::~HLTHcalSimpleRecHitFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HLTHcalSimpleRecHitFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // using namespace edm;

    // getting very basic uncalRH
    edm::Handle<HFRecHitCollection> crudeHits;
    try {
        iEvent.getByLabel(HcalRecHitCollection_, crudeHits);
    } catch ( std::exception& ex) {
        edm::LogWarning("HLTHcalSimpleRecHitFilter") << HcalRecHitCollection_ << " not available";
    }
    
    bool aboveThreshold = false ; 

    for ( HFRecHitCollection::const_iterator hitItr = crudeHits->begin(); hitItr != crudeHits->end(); ++hitItr ) {     
       HFRecHit hit = (*hitItr);
     
       // masking noisy channels
       std::vector<int>::iterator result;
       result = find( maskedList_.begin(), maskedList_.end(), HcalDetId(hit.id()).hashed_index() );    
       if  (result != maskedList_.end()) 
           continue; 
       
       float ampli_ = hit.energy();
       
       if ( ampli_ >= threshold_ ) {
           aboveThreshold = true;
//            edm::LogInfo("HcalFilter")  << "at evet: " << iEvent.id().event() 
//                                           << " and run: " << iEvent.id().run() 
//                                           << " signal above threshold found: " << ampli_ ; 
           break ;
       }
    }
    return aboveThreshold ; 
}

// ------------ method called once each job just before starting event loop  ------------
void 
HLTHcalSimpleRecHitFilter::beginJob(const edm::EventSetup&) {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTHcalSimpleRecHitFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTHcalSimpleRecHitFilter);
