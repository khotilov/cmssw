// File: DataMixingEMWorker.cc
// Description:  see DataMixingEMWorker.h
// Author:  Mike Hildreth, University of Notre Dame
//
//--------------------------------------------

#include <map>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/BranchDescription.h"
//
//
#include "DataMixingEMWorker.h"


using namespace std;

namespace edm
{

  // Virtual constructor

  DataMixingEMWorker::DataMixingEMWorker() { sel_=0;}

  // Constructor 
  DataMixingEMWorker::DataMixingEMWorker(const edm::ParameterSet& ps) : 
							    label_(ps.getParameter<std::string>("Label"))

  {                                                         

    // get the subdetector names
    //    this->getSubdetectorNames();  //something like this may be useful to check what we are supposed to do...

    // create input selector
    if (label_.size()>0){
      sel_=new Selector( ModuleLabelSelector(label_));
    }
    else {
      sel_=new Selector( MatchAllSelector());
    }

    // declare the products to produce, retrieve

    EBProducer_ = ps.getParameter<edm::InputTag>("EBProducer");
    EEProducer_ = ps.getParameter<edm::InputTag>("EEProducer");
    ESProducer_ = ps.getParameter<edm::InputTag>("ESProducer");
    EBrechitCollection_ = ps.getParameter<edm::InputTag>("EBrechitCollection");
    EErechitCollection_ = ps.getParameter<edm::InputTag>("EErechitCollection");
    ESrechitCollection_ = ps.getParameter<edm::InputTag>("ESrechitCollection");
    EBRecHitCollectionDM_        = ps.getParameter<std::string>("EBRecHitCollectionDM");
    EERecHitCollectionDM_        = ps.getParameter<std::string>("EERecHitCollectionDM");
    ESRecHitCollectionDM_        = ps.getParameter<std::string>("ESRecHitCollectionDM");
   //   nMaxPrintout_            = ps.getUntrackedParameter<int>("nMaxPrintout",10);

   //EBalgo_ = new EcalRecHitSimpleAlgo();
   //EEalgo_ = new EcalRecHitSimpleAlgo();

   // don't think I can "produce" in a sub-class...

   //produces< EBRecHitCollection >(EBRecHitCollectionDM_);
   //produces< EERecHitCollection >(EERecHitCollectionDM_);

  }
	       

  // Virtual destructor needed.
  DataMixingEMWorker::~DataMixingEMWorker() { 
    delete sel_;
    sel_=0;
  }  

  void DataMixingEMWorker::addEMSignals(const edm::Event &e) { 
    // fill in maps of hits

    LogInfo("DataMixingEMWorker")<<"===============> adding MC signals for "<<e.id();

    // EB first

   Handle< EBRecHitCollection > pEBRecHits;

   const EBRecHitCollection*  EBRecHits = 0;

   try {
     e.getByLabel(EBProducer_.label(),EBrechitCollection_.label(), pEBRecHits);
     EBRecHits = pEBRecHits.product(); // get a ptr to the product
     LogInfo("DataMixingEMWorker") << "total # EB rechits: " << EBRecHits->size();
   } catch (...) {
   }
 
   if (EBRecHits)
     {
       // loop over rechits, storing them in a map so we can add pileup later
       for(EBRecHitCollection::const_iterator it  = EBRecHits->begin();	
	   it != EBRecHits->end(); ++it) {

	 EBRecHitStorage_.insert(EBRecHitMap::value_type( ( it->id() ), *it ));
	 
         LogDebug("DataMixingEMWorker") << "processed EBRecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();

       }
     }

   // EE next

   Handle< EERecHitCollection > pEERecHits;

   const EERecHitCollection*  EERecHits = 0;

   try {
     e.getByLabel(EEProducer_.label(),EErechitCollection_.label(), pEERecHits);
     EERecHits = pEERecHits.product(); // get a ptr to the product
#ifdef DEBUG
     LogDebug("DataMixingEMWorker") << "total # EE rechits: " << EERecHits->size();
#endif
   } catch (...) {
   }
 
   if (EERecHits)
     {
       // loop over rechits, storing them in a map so we can add pileup later
       for(EERecHitCollection::const_iterator it  = EERecHits->begin();	
	   it != EERecHits->end(); ++it) {

	 EERecHitStorage_.insert(EERecHitMap::value_type( ( it->id() ), *it ));
	 
#ifdef DEBUG	 
         LogDebug("DataMixingEMWorker") << "processed EERecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();
#endif

       }
     }
   // ES next

   Handle< ESRecHitCollection > pESRecHits;

   const ESRecHitCollection*  ESRecHits = 0;

   try {
     e.getByLabel( ESProducer_.label(),ESrechitCollection_.label(), pESRecHits);
     ESRecHits = pESRecHits.product(); // get a ptr to the product
#ifdef DEBUG
     LogDebug("DataMixingEMWorker") << "total # ES rechits: " << ESRecHits->size();
#endif
   } catch (...) {
   }
 
   if (ESRecHits)
     {
       // loop over rechits, storing them in a map so we can add pileup later
       for(ESRecHitCollection::const_iterator it  = ESRecHits->begin();	
	   it != ESRecHits->end(); ++it) {

	 ESRecHitStorage_.insert(ESRecHitMap::value_type( ( it->id() ), *it ));
	 
#ifdef DEBUG	 
         LogDebug("DataMixingEMWorker") << "processed ESRecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();
#endif

       }
     }
    
  } // end of addEMSignals

  void DataMixingEMWorker::addEMPileups(const int bcr, Event *e, unsigned int eventNr) {
  
    LogInfo("DataMixingEMWorker") <<"\n===============> adding pileups from event  "<<e->id()<<" for bunchcrossing "<<bcr;

    // fill in maps of hits; same code as addSignals, except now applied to the pileup events

    // EB first

   Handle< EBRecHitCollection > pEBRecHits;
   const EBRecHitCollection*  EBRecHits = 0;

   try {
     e->getByLabel(EBProducer_.label(),EBrechitCollection_.label(), pEBRecHits);
     EBRecHits = pEBRecHits.product(); // get a ptr to the product
#ifdef DEBUG
     LogDebug("DataMixingEMWorker") << "total # EB rechits: " << EBRecHits->size();
#endif
   } catch (...) {
   }
 
   if (EBRecHits)
     {
       // loop over rechits, adding these to the existing maps
       for(EBRecHitCollection::const_iterator it  = EBRecHits->begin();
	   it != EBRecHits->end(); ++it) {

	 EBRecHitStorage_.insert(EBRecHitMap::value_type( (it->id()), *it ));
	 
#ifdef DEBUG	 
	 LogDebug("DataMixingEMWorker") << "processed EBRecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();
#endif
       }
     }
    // EE Next

   Handle< EERecHitCollection > pEERecHits;
   const EERecHitCollection*  EERecHits = 0;

   try {
     e->getByLabel( EEProducer_.label(),EErechitCollection_.label(), pEERecHits);
     EERecHits = pEERecHits.product(); // get a ptr to the product
#ifdef DEBUG
     LogDebug("DataMixingEMWorker") << "total # EE rechits: " << EERecHits->size();
#endif
   } catch (...) {
   }
 
   if (EERecHits)
     {
       // loop over rechits, adding these to the existing maps
       for(EERecHitCollection::const_iterator it  = EERecHits->begin();
	   it != EERecHits->end(); ++it) {

	 EERecHitStorage_.insert(EERecHitMap::value_type( (it->id()), *it ));
	 
#ifdef DEBUG	 
	 LogDebug("DataMixingEMWorker") << "processed EERecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();
#endif
       }
     }
    // ES Next

   Handle< ESRecHitCollection > pESRecHits;
   const ESRecHitCollection*  ESRecHits = 0;

   try {
     e->getByLabel( ESProducer_.label(),ESrechitCollection_.label(), pESRecHits);
     ESRecHits = pESRecHits.product(); // get a ptr to the product
#ifdef DEBUG
     LogDebug("DataMixingEMWorker") << "total # ES rechits: " << ESRecHits->size();
#endif
   } catch (...) {
   }
 
   if (ESRecHits)
     {
       // loop over rechits, adding these to the existing maps
       for(ESRecHitCollection::const_iterator it  = ESRecHits->begin();
	   it != ESRecHits->end(); ++it) {

	 ESRecHitStorage_.insert(ESRecHitMap::value_type( (it->id()), *it ));
	 
#ifdef DEBUG	 
	 LogDebug("DataMixingEMWorker") << "processed ESRecHit with rawId: "
				      << it->id().rawId() << "\n"
				      << " rechit energy: " << it->energy();
#endif
       }
     }

  }
 
  void DataMixingEMWorker::putEM(edm::Event &e) {

    // collection of rechits to put in the event
    std::auto_ptr< EBRecHitCollection > EBrechits( new EBRecHitCollection );
    std::auto_ptr< EERecHitCollection > EErechits( new EERecHitCollection );
    std::auto_ptr< ESRecHitCollection > ESrechits( new ESRecHitCollection );

    // loop over the maps we have, re-making individual hits or digis if necessary.
    DetId formerID = 0;
    DetId currentID;
    float ESum = 0.;
    float EBTime = 0.;
    EcalRecHit OldHit;
    int nmatch=0;

    // EB first...

    EBRecHitMap::const_iterator iEBchk;

    for(EBRecHitMap::const_iterator iEB  = EBRecHitStorage_.begin();
	iEB != EBRecHitStorage_.end(); ++iEB) {

      currentID = iEB->first; 

      if (currentID == formerID) { // we have to add these rechits together
	nmatch++;                  // use this to avoid using the "count" function
	ESum+=(iEB->second).energy();          // on every element...

	iEBchk = iEB;
	if((iEBchk++) == EBRecHitStorage_.end()) {  //make sure not to lose the last one
	  EcalRecHit aHit(formerID, ESum, EBTime);
	  EBrechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
      }
      else {
	if(nmatch>0) {
	  EcalRecHit aHit(formerID, ESum, EBTime);
	  EBrechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
	else {
	  EBrechits->push_back( OldHit );
	}
	
	iEBchk = iEB;
	if((iEBchk++) == EBRecHitStorage_.end()) {  //make sure not to lose the last one
	  EBrechits->push_back( iEB->second );
	}

	// save pointers for next iteration
	OldHit = iEB->second;
	formerID = currentID;
	ESum = (iEB->second).energy();
	EBTime = (iEB->second).time();  // take time of first hit in sequence - is this ok?
      }
    }

    // EE next...

    // loop over the maps we have, re-making individual hits or digis if necessary.
    formerID = 0;
    ESum = 0.;
    float EETime = 0.;
    
    nmatch=0;

    EERecHitMap::const_iterator iEEchk;

    for(EERecHitMap::const_iterator iEE  = EERecHitStorage_.begin();
	iEE != EERecHitStorage_.end(); ++iEE) {

      currentID = iEE->first; 

      if (currentID == formerID) { // we have to add these rechits together
	nmatch++;                  // use this to avoid using the "count" function
	ESum+=(iEE->second).energy();          // on every element...

	iEEchk = iEE;
	if((iEEchk++) == EERecHitStorage_.end()) {  //make sure not to lose the last one
	  EcalRecHit aHit(formerID, ESum, EETime);
	  EErechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
      }
      else {
	if(nmatch>0) {
	  EcalRecHit aHit(formerID, ESum, EETime);
	  EErechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
	else {
	  EErechits->push_back( OldHit );
	}
	
	iEEchk = iEE;
	if((iEEchk++) == EERecHitStorage_.end()) {  //make sure not to lose the last one
	  EErechits->push_back( iEE->second );
	}

	// save pointers for next iteration
	OldHit = iEE->second;
	formerID = currentID;
	ESum = (iEE->second).energy();
	EETime = (iEE->second).time();  // take time of first hit in sequence - is this ok?
      }
    }

    // ES next...

    // loop over the maps we have, re-making individual hits or digis if necessary.
    formerID = 0;
    ESum = 0.;
    float ESTime = 0.;
    nmatch=0;

    ESRecHitMap::const_iterator iESchk;

    for(ESRecHitMap::const_iterator iES  = ESRecHitStorage_.begin();
	iES != ESRecHitStorage_.end(); ++iES) {

      currentID = iES->first; 

      if (currentID == formerID) { // we have to add these rechits together
	nmatch++;                  // use this to avoid using the "count" function
	ESum+=(iES->second).energy();          // on every element...

	iESchk = iES;
	if((iESchk++) == ESRecHitStorage_.end()) {  //make sure not to lose the last one
	  EcalRecHit aHit(formerID, ESum, ESTime);
	  ESrechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
      }
      else {
	if(nmatch>0) {
	  EcalRecHit aHit(formerID, ESum, ESTime);
	  ESrechits->push_back( aHit );	  
	  // reset energy sum, nmatch
	  ESum = 0 ;
	  nmatch=0;	  
	}
	else {
	  ESrechits->push_back( OldHit );
	}
	
	iESchk = iES;
	if((iESchk++) == ESRecHitStorage_.end()) {  //make sure not to lose the last one
	  ESrechits->push_back( iES->second );
	}

	// save pointers for next iteration
	OldHit = iES->second;
	formerID = currentID;
	ESum = (iES->second).energy();
	ESTime = (iES->second).time();  // take time of first hit in sequence - is this ok?
      }
    }

    // put the collection of reconstructed hits in the event   
    LogInfo("DataMixingEMWorker") << "total # EB Merged rechits: " << EBrechits->size() ;
    LogInfo("DataMixingEMWorker") << "total # EE Merged rechits: " << EErechits->size() ;
    LogInfo("DataMixingEMWorker") << "total # ES Merged rechits: " << ESrechits->size() ;

    e.put( EBrechits, EBRecHitCollectionDM_ );
    e.put( EErechits, EERecHitCollectionDM_ );
    e.put( ESrechits, ESRecHitCollectionDM_ );
    
    // clear local storage after this event

    EBRecHitStorage_.clear();
    EERecHitStorage_.clear();
    ESRecHitStorage_.clear();

  }

} //edm
