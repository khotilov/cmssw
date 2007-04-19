#ifndef L1GCTEMULATOR_H
#define L1GCTEMULATOR_H

/**\class L1GctEmulator L1GctEmulator.h src/L1Trigger/GlobalCaloTrigger/src/L1GctEmulator.h

 Description:  Framework module that runs the GCT bit-level emulator

 Implementation:
       An EDProducer that contains an instance of L1GlobalCaloTrigger.
*/
//
// Original Author:  Jim Brooke
//         Created:  Thu May 18 15:04:56 CEST 2006
// $Id: L1GctEmulator.h,v 1.7 2007/03/16 14:44:06 heath Exp $
//
//

// system includes


// EDM includes
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Utilities/interface/Exception.h"

// GCT includes
#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetEtCalibrationLut.h"

class L1GctEmulator : public edm::EDProducer {
 public:

  /// constructor
  explicit L1GctEmulator(const edm::ParameterSet& ps);

  /// destructor
  virtual ~L1GctEmulator();

 private:
  virtual void beginJob(const edm::EventSetup& c) ;
  virtual void produce(edm::Event& e, const edm::EventSetup& c);
  virtual void endJob() ;

  void configureGct(const edm::EventSetup& c) ;

  // input label
  std::string m_inputLabel;

  // pointer to the actual emulator
  L1GlobalCaloTrigger* m_gct;

  // the configuration Lut
  L1GctJetEtCalibrationLut* m_jetEtCalibLut;

  // untracked parameters
  bool m_verbose;

  // tracked parameters

};

#endif
