// Last commit: $Id: LatencyHistosUsingDb.cc,v 1.1 2007/12/11 16:09:57 delaer Exp $

#include "DQM/SiStripCommissioningDbClients/interface/LatencyHistosUsingDb.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include <iostream>

using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
LatencyHistosUsingDb::LatencyHistosUsingDb( MonitorUserInterface* mui,
					      const DbParams& params )
  : LatencyHistograms( mui ),
    CommissioningHistosUsingDb( params )
{
  LogTrace(mlDqmClient_) 
    << "[LatencyHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
LatencyHistosUsingDb::LatencyHistosUsingDb( MonitorUserInterface* mui,
					      SiStripConfigDb* const db )
  : LatencyHistograms( mui ),
    CommissioningHistosUsingDb( db )
{
  LogTrace(mlDqmClient_) 
    << "[LatencyHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
LatencyHistosUsingDb::LatencyHistosUsingDb( DaqMonitorBEInterface* bei,
					      SiStripConfigDb* const db ) 
  : LatencyHistograms( bei ),
    CommissioningHistosUsingDb( db )
{
  LogTrace(mlDqmClient_) 
    << "[LatencyHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
LatencyHistosUsingDb::~LatencyHistosUsingDb() {
  LogTrace(mlDqmClient_) 
    << "[LatencyHistosUsingDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
void LatencyHistosUsingDb::uploadToConfigDb() {
  
  if ( !db_ ) {
    edm::LogWarning(mlDqmClient_) 
      << "[LatencyHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripConfigDb interface!"
      << " Aborting upload...";
    return;
  }

  // Update APV descriptions with new Latency settings
  const SiStripConfigDb::DeviceDescriptions& devices = db_->getDeviceDescriptions(); 
  update( const_cast<SiStripConfigDb::DeviceDescriptions&>(devices) );
  if ( !test_ ) { 
    LogTrace(mlDqmClient_) 
      << "[LatencyHistosUsingDb::" << __func__ << "]"
      << " Uploading APV settings to DB...";
    db_->uploadDeviceDescriptions(true); 
    LogTrace(mlDqmClient_) 
      << "[LatencyHistosUsingDb::" << __func__ << "]"
      << " Upload of APV settings to DB finished!";
  } else {
    edm::LogWarning(mlDqmClient_) 
      << "[LatencyHistosUsingDb::" << __func__ << "]"
      << " TEST only! No APV settings will be uploaded to DB...";
  }
  
}

// -----------------------------------------------------------------------------
/** */
void LatencyHistosUsingDb::update( SiStripConfigDb::DeviceDescriptions& devices ) {
  
  // Obtain the latency from the analysis object
  if(!data_.size() || !data_.begin()->second.isValid() ) {
    edm::LogVerbatim(mlDqmClient_) 
      << "[LatencyHistosUsingDb::" << __func__ << "]"
      << " Updated NO Latency settings. No analysis result available !" ;
    return;
  }
  uint16_t latency = uint16_t((data_.begin()->second.maximum()/(-25.))+0.5);
  
  // Iterate through devices and update device descriptions
  uint16_t updated = 0;
  SiStripConfigDb::DeviceDescriptions::iterator idevice;
  for ( idevice = devices.begin(); idevice != devices.end(); idevice++ ) {
    // Check device type
    if ( (*idevice)->getDeviceType() != APV25 ) { continue; }
    // Cast to retrieve appropriate description object
    apvDescription* desc = dynamic_cast<apvDescription*>( *idevice );
    if ( !desc ) { continue; }
    // Retrieve device addresses from device description
    const SiStripConfigDb::DeviceAddress& addr = db_->deviceAddress(*desc);
    // Do it!
    std::stringstream ss;
    ss << "[LatencyHistosUsingDb::" << __func__ << "]"
       << " Updating latency APV settings for crate/FEC/slot/ring/CCU/i2cAddr "
       << addr.fecCrate_ << "/"
       << addr.fecSlot_ << "/"
       << addr.fecRing_ << "/"
       << addr.ccuAddr_ << "/"
       << addr.ccuChan_ << "/"
       << addr.i2cAddr_
       << " from "
       << static_cast<uint16_t>(desc->getLatency());
    desc->setLatency(latency);
    updated++;
    ss << " to "
       << static_cast<uint16_t>(desc->getLatency());
    LogTrace(mlDqmClient_) << ss.str();
  }
  edm::LogVerbatim(mlDqmClient_) 
    << "[LatencyHistosUsingDb::" << __func__ << "] "
    << "Updated Latency settings for " << updated << " devices";
}

