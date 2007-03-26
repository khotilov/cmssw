/** \class HLTHighLevel
 *
 * See header file for documentation
 *
 *  $Date: 2007/03/26 11:31:42 $
 *  $Revision: 1.1 $
 *
 *  \author Martin Grunewald
 *
 */

#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/HLTFilterObject.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <cassert>

//
// constructors and destructor
//
HLTHighLevel::HLTHighLevel(const edm::ParameterSet& iConfig) :
  inputTag_ (iConfig.getParameter<edm::InputTag> ("TriggerResultsTag")),
  andOr_    (iConfig.getParameter<bool> ("andOr" )),
  byName_   (iConfig.getParameter<bool> ("byName")),
  n_        (0)
{
  if (byName_) {
    // get names, then derive slot numbers
    HLTPathsByName_= iConfig.getParameter<std::vector<std::string > >("HLTPaths");
    n_=HLTPathsByName_.size();
    HLTPathsByIndex_.resize(n_);
  } else {
    // get slot numbers, then derive names
    HLTPathsByIndex_= iConfig.getParameter<std::vector<unsigned int> >("HLTPaths");
    n_=HLTPathsByIndex_.size();
    HLTPathsByName_.resize(n_);
  }

  // this is a user/analysis filter: it places no product into the event!

}

HLTHighLevel::~HLTHighLevel()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
bool
HLTHighLevel::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;

   const string invalid("@@invalid@@");

   // get hold of TriggerResults Object
   Handle<TriggerResults> trh;
   try {iEvent.getByLabel(inputTag_,trh);} catch(...) {;}
   if (trh.isValid()) {
     LogDebug("") << "TriggerResults found, number of HLT paths: " << trh->size();
   } else {
     LogDebug("") << "TriggerResults product not found - returning result=false!";
     return false;
   }

   // use event data to get the current HLT trigger table
   // this is an ugly hack, the (possibly changing) HLT trigger table 
   // should rather be taken from some runBlock or lumiBlock on file

   unsigned int n(n_);
   if (byName_) {
     for (unsigned int i=0; i!=n; i++) {
       HLTPathsByIndex_[i]=trh->find(HLTPathsByName_[i]);
     }
   } else {
     for (unsigned int i=0; i!=n; i++) {
       if (HLTPathsByIndex_[i]<trh->size()) {
	 HLTPathsByName_[i]=trh->name(HLTPathsByIndex_[i]);
       } else {
	 HLTPathsByName_[i]=invalid;
       }
     }
   }
   
   // for empty input vectors (n==0), default to all HLT trigger paths!
   if (n==0) {
     n=trh->size();
     HLTPathsByName_.resize(n);
     HLTPathsByIndex_.resize(n);
     for (unsigned int i=0; i!=n; i++) {
       HLTPathsByName_[i]=trh->name(i);
       HLTPathsByIndex_[i]=i;
     }
   }

   // report on what is finally used
   LogDebug("") << "HLT trigger paths: " + inputTag_.encode()
		<< " - Number requested: " << n
		<< " - andOr mode: " << andOr_
		<< " - byName: " << byName_;
   if (n>0) {
     LogDebug("") << "  HLT trigger paths requested: index, name and valididty:";
     for (unsigned int i=0; i!=n; i++) {
       LogTrace("") << " " << HLTPathsByIndex_[i]
		    << " " << HLTPathsByName_[i]
		    << " " << ( (HLTPathsByIndex_[i]<trh->size()) && (HLTPathsByName_[i]!=invalid) );
     }
   }

   // count number of requested HLT paths which have fired
   unsigned int fired(0);
   for (unsigned int i=0; i!=n; i++) {
     if (HLTPathsByIndex_[i]<trh->size()) {
       if (trh->accept(HLTPathsByIndex_[i])) {
	 fired++;
       }
     }
   }

   // Boolean filter result
   const bool accept( ((!andOr_) && (fired==n)) ||
		      (( andOr_) && (fired!=0)) );
   LogDebug("") << "Accept = " << accept;

   return accept;

}
