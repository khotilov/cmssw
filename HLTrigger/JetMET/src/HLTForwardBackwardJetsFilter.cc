/** \class HLTForwardBackwardJetsFilter
 *
 * $Id: HLTForwardBackwardJetsFilter.cc,v 1.5 2008/06/07 19:46:08 apana Exp $
 *
 *
 */

#include "HLTrigger/JetMET/interface/HLTForwardBackwardJetsFilter.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


//
// constructors and destructor
//
HLTForwardBackwardJetsFilter::HLTForwardBackwardJetsFilter(const edm::ParameterSet& iConfig)
{
   inputTag_    = iConfig.getParameter< edm::InputTag > ("inputTag");
   saveTag_     = iConfig.getUntrackedParameter<bool>("saveTag",false);
   minPt_       = iConfig.getParameter<double> ("minPt");
   minEta_ = iConfig.getParameter<double> ("minEta"); 
   maxEta_ = iConfig.getParameter<double> ("maxEta"); 

   //register your products
   produces<trigger::TriggerFilterObjectWithRefs>();
}

HLTForwardBackwardJetsFilter::~HLTForwardBackwardJetsFilter(){}


// ------------ method called to produce the data  ------------
bool
HLTForwardBackwardJetsFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace trigger;
  // The filter object
  std::auto_ptr<trigger::TriggerFilterObjectWithRefs> 
    filterobject (new trigger::TriggerFilterObjectWithRefs(path(),module()));
  if (saveTag_) filterobject->addCollectionTag(inputTag_);

  edm::Handle<reco::CaloJetCollection> recocalojets;
  iEvent.getByLabel(inputTag_,recocalojets);

  // look at all candidates,  check cuts and add to filter object
  int nplusjets(0);
  int nminusjets(0);

  if(recocalojets->size() > 1){
    // events with two or more jets

    // look for jets satifying pt and eta cuts; first on the plus side, then the minus side
    for (reco::CaloJetCollection::const_iterator recocalojet = recocalojets->begin(); 
	 recocalojet!=(recocalojets->end()); recocalojet++) {

      float ptjet=recocalojet->pt();
      float etajet=recocalojet->eta();
      if( ptjet > minPt_ ){
	if ( etajet > minEta_ && etajet < maxEta_ ){
	  nplusjets++;
	  reco::CaloJetRef ref(reco::CaloJetRef(recocalojets,distance(recocalojets->begin(),recocalojet)));
	  filterobject->addObject(TriggerJet,ref);
	}
      }
    }
    if (nplusjets > 0) {   
      for (reco::CaloJetCollection::const_iterator recocalojet = recocalojets->begin(); 
	   recocalojet!=(recocalojets->end()); recocalojet++) {

	float ptjet=recocalojet->pt();
	float etajet=recocalojet->eta();

	if( ptjet > minPt_ ){
	  if ( etajet < -minEta_ && etajet > -maxEta_ ){
	    nminusjets++;
	    reco::CaloJetRef ref(reco::CaloJetRef(recocalojets,distance(recocalojets->begin(),recocalojet)));
	    filterobject->addObject(TriggerJet,ref);
	  }
	}
      }
    }
  } // events with two or more jets
  
  
  
  // filter decision
  bool accept(nplusjets>0 && nminusjets>0);  
  // put filter object into the Event
  iEvent.put(filterobject);
  
  return accept;
}
