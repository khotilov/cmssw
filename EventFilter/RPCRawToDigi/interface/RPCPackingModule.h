#ifndef RPCRawToDigi_RPCPackingModule_H
#define RPCRawToDigi_RPCPackingModule_H

/** \class RPCPackingModule
 *  Driver class for digi to raw data conversions 
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "EventFilter/RPCRawToDigi/interface/EventRecords.h"

#include <vector> 

namespace edm {class ParameterSet;}
namespace edm {class EventSetup; }
namespace edm {class Event; }

class FEDRawData;
class RPCRecordFormatter;

class RPCPackingModule : public edm::EDProducer {
public:

  /// ctor
  explicit RPCPackingModule( const edm::ParameterSet& );

  /// dtor
  virtual ~RPCPackingModule();

  /// get data, convert to raw event, attach again to Event
  virtual void produce( edm::Event&, const edm::EventSetup& );

  static std::vector<rpcrawtodigi::EventRecords> eventRecords(
      int fedId, int trigger_BX, const RPCDigiCollection* , const RPCRecordFormatter& ); 

private:

  FEDRawData * rawData( int fedId, const RPCDigiCollection* , const RPCRecordFormatter& );

private:

//  edm::InputTag digiLabel_;
  unsigned long eventCounter_;
};
#endif
