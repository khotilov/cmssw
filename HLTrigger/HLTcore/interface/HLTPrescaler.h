#ifndef HLTPrescaler_H
#define HLTPrescaler_H

/** \class HLTPrescaler
 *
 *  
 *  This class is an EDFilter implementing an HLT
 *  Prescaler module with associated book keeping.
 *
 *  $Date: 2012/01/22 22:32:14 $
 *  $Revision: 1.21 $
 *
 *  \author Martin Grunewald
 *  \author Philipp Schieferdecker
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PrescaleService/interface/PrescaleService.h"


class HLTPrescaler : public edm::EDFilter
{
public:
  //
  // construction/destruction
  //
  explicit HLTPrescaler(edm::ParameterSet const& iConfig);
  virtual ~HLTPrescaler();


  //
  // member functions
  //
  virtual bool beginLuminosityBlock(edm::LuminosityBlock &lb,
				    edm::EventSetup const& iSetup);
  virtual bool filter(edm::Event& iEvent,edm::EventSetup const& iSetup);
  virtual void endJob();
  
  
private:
  //
  //member data
  //

  /// accept one in prescaleFactor_; 0 means never to accept an event
  unsigned int prescaleFactor_;

  /// event counter
  unsigned int eventCount_;

  /// accept counter
  unsigned int acceptCount_;

  /// initial offset
  unsigned int offsetCount_;
  unsigned int offsetPhase_;

  /// prescale service
  edm::service::PrescaleService* prescaleService_;
  
  /// check for (re)initialization of the prescale
  bool newLumi_;

  /// GT payload, to extract the prescale column index
  edm::InputTag gtDigi_;

  /// "seed" used to initialize the prescale counter
  static const
  unsigned int prescaleSeed_;

};

#endif
