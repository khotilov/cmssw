// Last commit: $Id: test_NoiseBuilder.h,v 1.3 2008/05/26 13:37:26 giordano Exp $
// Latest tag:  $Name: HEAD $
// Location:    $Source: /cvs_server/repositories/CMSSW/CMSSW/OnlineDB/SiStripESSources/test/stubs/test_NoiseBuilder.h,v $

#ifndef OnlineDB_SiStripESSources_test_NoiseBuilder_H
#define OnlineDB_SiStripESSources_test_NoiseBuilder_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/**
   @class test_NoiseBuilder 
   @brief Simple class that analyzes Digis produced by RawToDigi unpacker
*/
class test_NoiseBuilder : public edm::EDAnalyzer {

 public:
  
  test_NoiseBuilder( const edm::ParameterSet& ) {;}
  virtual ~test_NoiseBuilder() {;}
  
  void beginJob() {;}
  void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob() {;}
  
};

#endif // OnlineDB_SiStripESSources_test_NoiseBuilder_H

