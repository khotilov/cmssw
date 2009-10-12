// $Id: EventConsumerSelector.cc,v 1.4 2009/08/18 19:26:25 biery Exp $
/// @file: EventConsumerSelector.cc

#include <vector>

#include <boost/lambda/lambda.hpp>

#include "EventFilter/StorageManager/interface/EventConsumerSelector.h"
#include "EventFilter/StorageManager/interface/Exception.h"

using namespace stor;

void EventConsumerSelector::initialize( const InitMsgView& imv )
{

  if( _initialized ) return;

  if( _configInfo.outputModuleLabel() != imv.outputModuleLabel() ) return; 

  _outputModuleId = imv.outputModuleId();

  edm::ParameterSet pset;
  pset.addParameter<Strings>( "SelectEvents", _configInfo.selEvents() );

  Strings tnames;
  imv.hltTriggerNames( tnames );

  std::ostringstream errorMsg;
  errorMsg << "Cannot initialize edm::EventSelector for consumer" <<
    _configInfo.consumerName() << " running on " << _configInfo.remoteHost() <<
    " requesting output module ID" << _outputModuleId <<
    " with label " << _configInfo.outputModuleLabel() <<
    " and HLT trigger names";
  std::for_each(tnames.begin(), tnames.end(), errorMsg << boost::lambda::constant(" ") << _1);
  try
  {
    _eventSelector.reset( new edm::EventSelector( pset, tnames ) );
  }
  catch ( edm::Exception& e )
  {
    errorMsg << e.what();
    
    XCEPT_RAISE(stor::exception::InvalidEventSelection, errorMsg.str());
  }
  catch( std::exception &e )
  {
    errorMsg << e.what();

    XCEPT_RAISE(stor::exception::InvalidEventSelection, errorMsg.str());
  }
  catch(...)
  {
    errorMsg << "Unknown exception";

    XCEPT_RAISE(stor::exception::InvalidEventSelection, errorMsg.str());
  }

  _initialized = true;

}

bool EventConsumerSelector::acceptEvent( const I2OChain& ioc )
{

  if( !_initialized ) return false;
  if( _stale ) return false;

  if( ioc.outputModuleId() != _outputModuleId ) return false;

  std::vector<unsigned char> hlt_out;
  ioc.hltTriggerBits( hlt_out );

  return _eventSelector->wantAll()
    || _eventSelector->acceptEvent( &hlt_out[0], ioc.hltTriggerCount() );

}


/// emacs configuration
/// Local Variables: -
/// mode: c++ -
/// c-basic-offset: 2 -
/// indent-tabs-mode: nil -
/// End: -
