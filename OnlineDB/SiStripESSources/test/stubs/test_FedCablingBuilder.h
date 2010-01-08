// Last commit: $Id: test_FedCablingBuilder.h,v 1.6 2008/05/26 13:37:26 giordano Exp $
// Latest tag:  $Name: HEAD $
// Location:    $Source: /cvs_server/repositories/CMSSW/CMSSW/OnlineDB/SiStripESSources/test/stubs/test_FedCablingBuilder.h,v $

#ifndef OnlineDB_SiStripESSources_test_FedCablingBuilder_H
#define OnlineDB_SiStripESSources_test_FedCablingBuilder_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/**
   @class test_FedCablingBuilder 
   @brief Simple class that analyzes Digis produced by RawToDigi unpacker
*/
class test_FedCablingBuilder : public edm::EDAnalyzer {

 public:
  
  test_FedCablingBuilder( const edm::ParameterSet& ) {;}
  virtual ~test_FedCablingBuilder() {;}
  
  void beginJob(){;}
  void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob() {;}
  
};

#endif // OnlineDB_SiStripESSources_test_FedCablingBuilder_H

