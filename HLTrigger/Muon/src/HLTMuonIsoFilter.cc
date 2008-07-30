/** \class HLTMuonIsoFilter
 *
 * See header file for documentation
 *
 *  \author J. Alcaraz
 *
 */

#include "HLTrigger/Muon/interface/HLTMuonIsoFilter.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include <iostream>
//
// constructors and destructor
//
HLTMuonIsoFilter::HLTMuonIsoFilter(const edm::ParameterSet& iConfig) :
   candTag_ (iConfig.getParameter< edm::InputTag > ("CandTag") ),
   previousCandTag_ (iConfig.getParameter<edm::InputTag > ("PreviousCandTag")),
   isoTag_  (iConfig.getParameter< edm::InputTag > ("IsoTag" ) ),
   min_N_   (iConfig.getParameter<int> ("MinN")),
   saveTag_  (iConfig.getUntrackedParameter<bool> ("SaveTag",false)) 
{
   LogDebug("HLTMuonIsoFilter") << " candTag : " << candTag_.encode()
      << "  IsoTag : " << isoTag_.encode()
      << "  MinN : " << min_N_;

   //register your products
   produces<trigger::TriggerFilterObjectWithRefs>();
}

HLTMuonIsoFilter::~HLTMuonIsoFilter()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool
HLTMuonIsoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace trigger;
   using namespace reco;

   // All HLT filters must create and fill an HLT filter object,
   // recording any reconstructed physics objects satisfying (or not)
   // this HLT filter, and place it in the Event.

   // The filter object
   auto_ptr<TriggerFilterObjectWithRefs>
     filterproduct (new TriggerFilterObjectWithRefs(path(),module()));

   // get hold of trks
   Handle<RecoChargedCandidateCollection> mucands;
   if(saveTag_)filterproduct->addCollectionTag(candTag_);
   iEvent.getByLabel (candTag_,mucands);
   Handle<TriggerFilterObjectWithRefs> previousLevelCands;
   iEvent.getByLabel (previousCandTag_,previousLevelCands);
   vector<RecoChargedCandidateRef> vcands;
   previousLevelCands->getObjects(TriggerMuon,vcands);
   
   //get hold of energy deposition
   Handle<edm::ValueMap<bool> > depMap;
   iEvent.getByLabel (isoTag_,depMap);
   
   // look at all mucands,  check cuts and add to filter object
   int n = 0;
   for (unsigned int i=0; i<mucands->size(); i++) {
     RecoChargedCandidateRef candref(mucands,i);
     
     //did this candidate triggered at previous stage.
     if (!triggerdByPreviousLevel(candref,vcands)) continue;
   
     TrackRef tk = candref->get<TrackRef>();
     edm::ValueMap<bool> ::value_type muonIsIsolated = (*depMap)[tk];
     LogDebug("HLTMuonIsoFilter") << " Muon with q*pt= " << tk->charge()*tk->pt() << ", eta= " << tk->eta() << "; Is Muon isolated? " << muonIsIsolated;
     
     if (!muonIsIsolated) continue;
     
     n++;
     filterproduct->addObject(TriggerMuon,candref);
   }

   // filter decision
   const bool accept (n >= min_N_);

   // put filter object into the Event
   iEvent.put(filterproduct);

   LogDebug("HLTMuonIsoFilter") << " >>>>> Result of HLTMuonIsoFilter is " << accept << ", number of muons passing isolation cuts= " << n; 

   return accept;
}

bool HLTMuonIsoFilter::triggerdByPreviousLevel(const reco::RecoChargedCandidateRef & candref, const std::vector<reco::RecoChargedCandidateRef>& vcands){
  bool ok=false;
  uint i=0;
  uint i_max=vcands.size();
  for (;i!=i_max;++i){
    if (candref == vcands[i]) { ok=true; break;}
  }

  return ok;
}
																						       
