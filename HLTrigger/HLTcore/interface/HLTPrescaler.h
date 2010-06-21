#ifndef HLTPrescaler_H
#define HLTPrescaler_H

/** \class HLTPrescaler
 *
 *  
 *  This class is an HLTFilter (-> EDFilter) implementing an HLT
 *  Prescaler module with associated book keeping.
 *
 *  $Date: 2009/12/02 18:47:06 $
 *  $Revision: 1.16 $
 *
 *  \author Martin Grunewald
 *  \author Philipp Schieferdecker
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "FWCore/PrescaleService/interface/PrescaleService.h"


class HLTPrescaler : public HLTFilter
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

  /// prescale service
  edm::service::PrescaleService* prescaleService_;
  
  /// check for (re)initialization of the prescale
  bool newLumi_;

  /// GT payload, to extract the prescale column index
  edm::InputTag gtDigi_;

};

#endif
