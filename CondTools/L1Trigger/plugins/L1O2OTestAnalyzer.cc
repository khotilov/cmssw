// -*- C++ -*-
//
// Package:    L1O2OTestAnalyzer
// Class:      L1O2OTestAnalyzer
// 
/**\class L1O2OTestAnalyzer L1O2OTestAnalyzer.cc CondTools/L1O2OTestAnalyzer/src/L1O2OTestAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Werner Man-Li Sun
//         Created:  Thu Nov  6 23:00:43 CET 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/L1TObjects/interface/L1TriggerKey.h"
#include "CondFormats/L1TObjects/interface/L1TriggerKeyList.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyRcd.h"
#include "CondFormats/DataRecord/interface/L1TriggerKeyListRcd.h"

#include "CondTools/L1Trigger/interface/Exception.h"


//
// class decleration
//

class L1O2OTestAnalyzer : public edm::EDAnalyzer {
   public:
      explicit L1O2OTestAnalyzer(const edm::ParameterSet&);
      ~L1O2OTestAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1O2OTestAnalyzer::L1O2OTestAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


L1O2OTestAnalyzer::~L1O2OTestAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1O2OTestAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   ESHandle< L1TriggerKeyList > pList ;
   iSetup.get< L1TriggerKeyListRcd >().get( pList ) ;

   std::cout << "Found " << pList->tscKeyToTokenMap().size() << " TSC keys:"
	     << std::endl ;

   L1TriggerKeyList::KeyToToken::const_iterator iTSCKey =
     pList->tscKeyToTokenMap().begin() ;
   L1TriggerKeyList::KeyToToken::const_iterator eTSCKey =
     pList->tscKeyToTokenMap().end() ;
   for( ; iTSCKey != eTSCKey ; ++iTSCKey )
     {
       std::cout << iTSCKey->first << " " << iTSCKey->second << std::endl ;
     }
   std::cout << std::endl ;

   L1TriggerKeyList::RecordToKeyToToken::const_iterator iRec =
     pList->recordTypeToKeyToTokenMap().begin() ;
   L1TriggerKeyList::RecordToKeyToToken::const_iterator eRec =
     pList->recordTypeToKeyToTokenMap().end() ;
   for( ; iRec != eRec ; ++iRec )
     {
       const L1TriggerKeyList::KeyToToken& keyTokenMap = iRec->second ;
       std::cout << "For record@type " << iRec->first << ", found "
		 << keyTokenMap.size() << " keys:" << std::endl ;

       L1TriggerKeyList::KeyToToken::const_iterator iKey = keyTokenMap.begin();
       L1TriggerKeyList::KeyToToken::const_iterator eKey = keyTokenMap.end() ;
       for( ; iKey != eKey ; ++iKey )
	 {
	   std::cout << iKey->first << " " << iKey->second << std::endl ;
	 }
       std::cout << std::endl ;
     }

   try
     {
       ESHandle< L1TriggerKey > pKey ;
       iSetup.get< L1TriggerKeyRcd >().get( pKey ) ;

       std::cout << "Current TSC key = " << pKey->getTSCKey() << std::endl ;

       const L1TriggerKey::RecordToKey& recKeyMap = pKey->recordToKeyMap() ;
       L1TriggerKey::RecordToKey::const_iterator iRec = recKeyMap.begin() ;
       L1TriggerKey::RecordToKey::const_iterator eRec = recKeyMap.end() ;
       for( ; iRec != eRec ; ++iRec )
	 {
	   std::cout << iRec->first << " " << iRec->second << std::endl ;
	 }
     }
   catch( cms::Exception& ex )
     {
       std::cout << "No L1TriggerKey found." << std::endl ;
     }

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1O2OTestAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1O2OTestAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1O2OTestAnalyzer);
