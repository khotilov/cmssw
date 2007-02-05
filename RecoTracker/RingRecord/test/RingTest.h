#ifndef RECOTRACKER_RINGTEST_H
#define RECOTRACKER_RINGTEST_H

//
// Package:         RecoTracker/RingMakerESProducer/test
// Class:           RingTest
// 
// Description:     test rings
//
// Original Author: Oliver Gutsche, gutsche@fnal.gov
// Created:         Fri Dec  8 10:15:02 UTC 2006
//
// $Author: noeding $
// $Date: 2006/08/12 00:23:51 $
// $Revision: 1.6 $
//

#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RingTest : public edm::EDAnalyzer {
 public:
  explicit RingTest( const edm::ParameterSet& );
  ~RingTest();
  
  virtual void analyze( const edm::Event&, const edm::EventSetup& );
 private:
  // ----------member data ---------------------------
  bool dumpRings_;
  std::string fileName_;
  
};

#endif
