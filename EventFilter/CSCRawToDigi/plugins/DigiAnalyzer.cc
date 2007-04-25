/** 
 * Demo analyzer for reading digis
 * author A.Tumanov 2/22/06 
 * ripped from Jeremy's and Rick's analyzers
 *   
 */
#include <iostream>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"
#include "EventFilter/CSCRawToDigi/interface/DigiAnalyzer.h"
#include "DataFormats/CSCDigi/interface/CSCWireDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigi.h"
#include "DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCRPCDigi.h"
#include "DataFormats/CSCDigi/interface/CSCRPCDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

DigiAnalyzer::DigiAnalyzer(edm::ParameterSet const& conf) {

  // If your module takes parameters, here is where you would define
  // their names and types, and access them to initialize internal
  // variables. Example as follows:
  //
  eventNumber = 0;
}

void DigiAnalyzer::analyze(edm::Event const& e, edm::EventSetup const& iSetup) {

  // These declarations create handles to the types of records that you want
  // to retrieve from event "e".
  //
  edm::Handle<CSCWireDigiCollection> wires;
  edm::Handle<CSCStripDigiCollection> strips;
  edm::Handle<CSCComparatorDigiCollection> comparators;
  edm::Handle<CSCALCTDigiCollection> alcts;
  edm::Handle<CSCCLCTDigiCollection> clcts;
  edm::Handle<CSCRPCDigiCollection> rpcs;
  edm::Handle<CSCCorrelatedLCTDigiCollection> correlatedlcts;
  
  // Pass the handle to the method "getByType", which is used to retrieve
  // one and only one instance of the type in question out of event "e". If
  // zero or more than one instance exists in the event an exception is thrown.
  //
  e.getByLabel("cscunpacker","MuonCSCWireDigi",wires);
  e.getByLabel("cscunpacker","MuonCSCStripDigi",strips);
  e.getByLabel("cscunpacker","MuonCSCComparatorDigi",comparators);
  e.getByLabel("cscunpacker","MuonCSCALCTDigi",alcts);
  e.getByLabel("cscunpacker","MuonCSCCLCTDigi",clcts);
  e.getByLabel("cscunpacker","MuonCSCRPCDigi",rpcs);
  e.getByLabel("cscunpacker","MuonCSCCorrelatedLCTDigi",correlatedlcts);
   
  
   
  // read digi collections and print digis
  //
  for (CSCWireDigiCollection::DigiRangeIterator j=wires->begin(); j!=wires->end(); j++) {
    std::vector<CSCWireDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCWireDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }
  
 
  for (CSCStripDigiCollection::DigiRangeIterator j=strips->begin(); j!=strips->end(); j++) {
    std::vector<CSCStripDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCStripDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }


  for (CSCComparatorDigiCollection::DigiRangeIterator j=comparators->begin(); j!=comparators->end(); j++) {
 
    std::vector<CSCComparatorDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCComparatorDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }



  for (CSCALCTDigiCollection::DigiRangeIterator j=alcts->begin(); j!=alcts->end(); j++) {
 
    std::vector<CSCALCTDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCALCTDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }



  for (CSCCLCTDigiCollection::DigiRangeIterator j=clcts->begin(); j!=clcts->end(); j++) {
 
    std::vector<CSCCLCTDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCCLCTDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }



  for (CSCRPCDigiCollection::DigiRangeIterator j=rpcs->begin(); j!=rpcs->end(); j++) {
 
    std::vector<CSCRPCDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCRPCDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }


  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator j=correlatedlcts->begin(); j!=correlatedlcts->end(); j++) {
 
    std::vector<CSCCorrelatedLCTDigi>::const_iterator digiItr = (*j).second.first;
    std::vector<CSCCorrelatedLCTDigi>::const_iterator last = (*j).second.second;
    for( ; digiItr != last; ++digiItr) {
      digiItr->print();
    }
  }



  eventNumber++;
  edm::LogInfo ("DigiAnalyzer")  << "end of event number " << eventNumber;
  
  

}





