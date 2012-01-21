/** \class HLTForwardBackwardJetsFilter
 *
 * $Id: HLTForwardBackwardJetsFilter.cc,v 1.5 2011/10/27 13:41:48 gruen Exp $
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

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// constructors and destructor
//
HLTForwardBackwardJetsFilter::HLTForwardBackwardJetsFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig) 
{
   inputTag_ = iConfig.getParameter< edm::InputTag > ("inputTag");
   minPt_    = iConfig.getParameter<double> ("minPt");
   minEta_   = iConfig.getParameter<double> ("minEta"); 
   maxEta_   = iConfig.getParameter<double> ("maxEta"); 
}

HLTForwardBackwardJetsFilter::~HLTForwardBackwardJetsFilter(){}

void
HLTForwardBackwardJetsFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("inputTag",edm::InputTag("hltIterativeCone5CaloJetsRegional"));
  desc.add<bool>("saveTags",false);
  desc.add<double>("minPt",15.0);
  desc.add<double>("minEta",3.0);
  desc.add<double>("maxEta",5.1);
  descriptions.add("hltForwardBackwardJetsFilter",desc);
}

// ------------ method called to produce the data  ------------
bool
HLTForwardBackwardJetsFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
{
  using namespace trigger;

  // The filter object
  if (saveTags()) filterproduct.addCollectionTag(inputTag_);

  edm::Handle<reco::CaloJetCollection> recocalojets;
  iEvent.getByLabel(inputTag_,recocalojets);

  // look at all candidates,  check cuts and add to filter object
  unsigned int nplusjets(0);
  unsigned int nminusjets(0);

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
	  filterproduct.addObject(TriggerJet,ref);
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
	    filterproduct.addObject(TriggerJet,ref);
	  }
	}
      }
    }
  } // events with two or more jets
  
  // filter decision
  bool accept(nplusjets>0 && nminusjets>0);  
  
  return accept;
}
