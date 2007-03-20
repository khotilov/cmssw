#ifndef L1ScalesProducers_L1ScalesTrivialProducer_h
#define L1ScalesProducers_L1ScalesTrivialProducer_h

// -*- C++ -*-
//
// Package:    L1ScalesProducers
// Class:      L1ScalesTrivialProducer
// 
/**\class L1ScalesTrivialProducer L1ScalesTrivialProducer.h L1TriggerConfig/L1ScalesProducers/interface/L1ScalesTrivialProducer.h

 Description: A Producer for the L1 scales available via EventSetup

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jim Brooke
//         Created:  Wed Sep 27 17:51:32 CEST 2006
// $Id: L1ScalesTrivialProducer.h,v 1.1 2006/12/12 12:37:03 jbrooke Exp $
//
//


// system include files
#include <memory>
#include <boost/shared_ptr.hpp>
#include <vector>

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/DataRecord/interface/L1EmEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"


//
// class declaration
//

class L1ScalesTrivialProducer : public edm::ESProducer {
public:
  L1ScalesTrivialProducer(const edm::ParameterSet&);
  ~L1ScalesTrivialProducer();
  
  std::auto_ptr<L1CaloEtScale> produceEmScale(const L1EmEtScaleRcd&);
  std::auto_ptr<L1CaloEtScale> produceJetScale(const L1JetEtScaleRcd&);

private:
  // ----------member data ---------------------------
  
  double m_emEtScaleInputLsb;
  std::vector<double> m_emEtThresholds;

  double m_jetEtScaleInputLsb;
  std::vector<double> m_jetEtThresholds;
  
};

#endif
