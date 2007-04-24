// Last commit: $Id: test_SiStripNullKey.h,v 1.3 2007/03/26 10:14:41 bainbrid Exp $

#ifndef DataFormats_SiStripCommon_test_SiStripNullKey_H
#define DataFormats_SiStripCommon_test_SiStripNullKey_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/**
   @class test_SiStripNullKey 
   @author R.Bainbridge
   @brief Simple class that tests SiStripNullKey.
*/
class test_SiStripNullKey : public edm::EDAnalyzer {

 public:
  
  test_SiStripNullKey( const edm::ParameterSet& );
  ~test_SiStripNullKey();
  
  void beginJob( edm::EventSetup const& );
  void analyze( const edm::Event&, const edm::EventSetup& );
  void endJob() {;}
  
};

#endif // DataFormats_SiStripCommon_test_SiStripNullKey_H

