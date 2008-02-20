// Last commit: $Id: FastFedCablingHistosUsingDb.cc,v 1.12 2008/02/19 21:17:18 bainbrid Exp $

#include "DQM/SiStripCommissioningDbClients/interface/FastFedCablingHistosUsingDb.h"
#include "CondFormats/SiStripObjects/interface/FastFedCablingAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include <iostream>

using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
FastFedCablingHistosUsingDb::FastFedCablingHistosUsingDb( MonitorUserInterface* mui,
							  const DbParams& params )
  : CommissioningHistosUsingDb( params ),
    FastFedCablingHistograms( mui )
{
  LogTrace(mlDqmClient_)
    << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
FastFedCablingHistosUsingDb::FastFedCablingHistosUsingDb( MonitorUserInterface* mui,
							  SiStripConfigDb* const db ) 
  : CommissioningHistograms( mui, sistrip::FAST_CABLING ),
    CommissioningHistosUsingDb( db, mui, sistrip::FAST_CABLING ),
    FastFedCablingHistograms( mui )
{
  LogTrace(mlDqmClient_)
    << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
FastFedCablingHistosUsingDb::FastFedCablingHistosUsingDb( DaqMonitorBEInterface* bei,
							  SiStripConfigDb* const db ) 
  : CommissioningHistosUsingDb( db, sistrip::FAST_CABLING ),
    FastFedCablingHistograms( bei )
{
  LogTrace(mlDqmClient_) 
    << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
FastFedCablingHistosUsingDb::~FastFedCablingHistosUsingDb() {
  LogTrace(mlDqmClient_)
    << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
void FastFedCablingHistosUsingDb::uploadConfigurations() {
  LogTrace(mlDqmClient_) 
    << "[FastFedCablingHistosUsingDb::" << __func__ << "]";
  
  if ( !db() ) {
    edm::LogError(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripConfigDb interface!"
      << " Aborting upload...";
    return;
  }

  // Retrieve descriptions for DCU id and DetId  
  SiStripConfigDb::DeviceDescriptions dcus = db()->getDeviceDescriptions( DCU ); 
  SiStripConfigDb::DcuDetIdMap detids = db()->getDcuDetIdMap(); 

  // Update FED connection descriptions
  const SiStripConfigDb::FedConnections& conns = db()->getFedConnections(); 
  update( const_cast<SiStripConfigDb::FedConnections&>(conns), dcus, detids );
  if ( doUploadConf() ) { 
    edm::LogVerbatim(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " Uploading FED connections to DB...";
    db()->uploadFedConnections(true); 
    edm::LogVerbatim(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " Completed database upload of " << conns.size() 
      << " ConnectionDescriptions!";
  } else {
    edm::LogWarning(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " TEST only! No FED connections will be uploaded to DB...";
  }
  
  // Update FED descriptions with enabled/disabled channels
  const SiStripConfigDb::FedDescriptions& feds = db()->getFedDescriptions(); 
  update( const_cast<SiStripConfigDb::FedDescriptions&>(feds) );
  if ( doUploadConf() ) { 
    edm::LogVerbatim(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " Uploading FED descriptions to DB...";
    db()->uploadFedDescriptions(true); 
    edm::LogVerbatim(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " Completed database upload of " << feds.size()
      << " Fed9UDescriptions (with connected channels enabled)!";
  } else {
    edm::LogWarning(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " TEST only! No FED descriptions will be uploaded to DB...";
  }

  // Some debug on good / dirty / missing connections
  connections( dcus, detids );
  
}

// -----------------------------------------------------------------------------
/** */
void FastFedCablingHistosUsingDb::update( SiStripConfigDb::FedConnections& conns,
					  const SiStripConfigDb::DeviceDescriptions& dcus, 
					  const SiStripConfigDb::DcuDetIdMap& detids ) {

  // Clear any connections retrieved from database
  if ( !conns.empty() ) {
    edm::LogWarning(mlDqmClient_)
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " Clearing existing FED channel connections!";
    conns.clear();
  }

  // Update FED-FEC mapping in base class, based on analysis results
  Analyses::iterator ianal = data().begin();
  Analyses::iterator janal = data().end();
  for ( ; ianal != janal; ++ianal ) {
    
    FastFedCablingAnalysis* anal = dynamic_cast<FastFedCablingAnalysis*>( ianal->second );
    if ( !anal ) { 
      edm::LogError(mlDqmClient_)
	<< "[FastFedCablingHistosUsingDb::" << __func__ << "]"
	<< " NULL pointer to analysis object!";
      continue; 
    }

    if ( !anal->isValid() ) { continue; }
    
    SiStripFecKey fec_key( anal->fecKey() );
    SiStripFedKey fed_key( anal->fedKey() );
    
#ifdef USING_NEW_DATABASE_MODEL	
    ConnectionDescription* conn = new ConnectionDescription();
    conn->setFedId( fed_key.fedId() );
    conn->setFedChannel( fed_key.fedChannel() );
    conn->setFecHardwareId( "" ); //@@
    conn->setFecCrateSlot( fec_key.fecCrate() );
    conn->setFecSlot( fec_key.fecSlot() );
    conn->setRingSlot( fec_key.fecRing() );
    conn->setCcuAddress( fec_key.ccuAddr() );
    conn->setI2cChannel( fec_key.ccuChan() );
    conn->setApvAddress( SiStripFecKey::i2cAddr(anal->lldCh(),true) );
    conn->setDcuHardId( anal->dcuHardId() );
#else
    FedChannelConnectionDescription* conn = new FedChannelConnectionDescription(); 
    conn->setFedId( fed_key.fedId() );
    conn->setFedChannel( fed_key.fedChannel() );
    conn->setFecSupervisor( "" );
    conn->setFecSupervisorIP( "" );
    conn->setFecInstance( fec_key.fecCrate() );
    conn->setSlot( fec_key.fecSlot() );
    conn->setRing( fec_key.fecRing() );
    conn->setCcu( fec_key.ccuAddr() );
    conn->setI2c( fec_key.ccuChan() );
    conn->setApv( SiStripFecKey::i2cAddr(anal->lldCh(),true) );
    conn->setDcuHardId( anal->dcuHardId() );
    SiStripConfigDb::DcuDetIdMap::const_iterator idet = detids.find( anal->dcuHardId() );
    if ( idet != detids.end() ) { 
      conn->setDetId( idet->second->getDetId() );
      conn->setApvPairs( idet->second->getApvNumber()/2 );
      conn->setFiberLength( static_cast<uint32_t>( idet->second->getFibreLength() ) );
    } else {
      conn->setDetId( sistrip::invalid32_ );
      conn->setApvPairs( sistrip::invalid_ ); 
      conn->setFiberLength( sistrip::invalid_ );
    }
#endif
    
    conns.push_back(conn);

  }  
  
#ifndef USING_NEW_DATABASE_MODEL	
  std::stringstream ss; 
  ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
     << " Dump of " << conns.size() 
     << " FedChannelConnection descriptions: "
     << std::endl;
  SiStripConfigDb::FedConnections::iterator ifed = conns.begin();
  for ( ; ifed != conns.end(); ifed++ ) { (*ifed)->toXML(ss); }
  LogTrace(mlTest_) << ss.str();
#endif
  
}

// -----------------------------------------------------------------------------
/** */
void FastFedCablingHistosUsingDb::update( SiStripConfigDb::FedDescriptions& feds ) {

  // Iterate through feds and disable all channels 
  SiStripConfigDb::FedDescriptions::iterator ifed;
  try {
    for ( ifed = feds.begin(); ifed != feds.end(); ifed++ ) {
      for ( uint16_t ichan = 0; ichan < sistrip::FEDCH_PER_FED; ichan++ ) {
	Fed9U::Fed9UAddress addr( ichan );
	Fed9U::Fed9UAddress addr0( ichan, static_cast<Fed9U::u8>(0) );
	Fed9U::Fed9UAddress addr1( ichan, static_cast<Fed9U::u8>(1) );
	(*ifed)->setFedFeUnitDisable( addr, true );
	(*ifed)->setApvDisable( addr0, true );
	(*ifed)->setApvDisable( addr1, true );
      }
    }
  } catch( ICUtils::ICException& e ) {
    edm::LogWarning(mlDqmClient_) << e.what();
  }
  
  // Counters for number of connected / enabled channels
  uint16_t connected = 0;
  std::map< uint16_t, std::vector<uint16_t> > enabled;

  // Iterate through feds and enable connected channels
  for ( ifed = feds.begin(); ifed != feds.end(); ifed++ ) {
    for ( uint16_t ichan = 0; ichan < sistrip::FEDCH_PER_FED; ichan++ ) {
      
      // Retrieve FEC key from FED-FEC map
      SiStripFedKey fed( static_cast<uint16_t>( (*ifed)->getFedId() ),
			 SiStripFedKey::feUnit(ichan),
			 SiStripFedKey::feChan(ichan) );
      uint32_t fed_key = fed.key();
      
      // Retrieve analysis for given FED id and channel
      Analyses::const_iterator iter = data().find( fed_key );
      if ( iter == data().end() ) { continue; }
      
      if ( !iter->second->isValid() ) { continue; }
      
      FastFedCablingAnalysis* anal = dynamic_cast<FastFedCablingAnalysis*>( iter->second );
      if ( !anal ) { 
	edm::LogError(mlDqmClient_)
	  << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
	  << " NULL pointer to OptoScanAnalysis object!";
	continue; 
      }
      
      // Retrieve FED id and channel 
      SiStripFedKey key( anal->fedKey() );
      uint16_t fed_id = key.fedId();
      uint16_t fed_ch = key.fedChannel();
      
      // Enable front-end unit and channel
      Fed9U::Fed9UAddress addr( fed_ch );
      Fed9U::Fed9UAddress addr0( fed_ch, static_cast<Fed9U::u8>(0) );
      Fed9U::Fed9UAddress addr1( fed_ch, static_cast<Fed9U::u8>(1) );
      (*ifed)->setFedFeUnitDisable( addr, false );
      (*ifed)->setApvDisable( addr0, false );
      (*ifed)->setApvDisable( addr1, false );
      connected++;
      enabled[fed_id].push_back(fed_ch);
      
    }
  }
      
  // Some debug 
  std::stringstream sss;
  if ( !feds.empty() ) {
    sss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
	<< " Enabled a total of " << connected 
	<< " FED channels and disabled " << feds.size() * 96 - connected 
	<< " FED channels (" << 100 * connected / ( feds.size() * 96 )
	<< "% of total)";
    edm::LogVerbatim(mlDqmClient_) << sss.str();
  } else {
    sss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
	<< " Found no FEDs! (and therefore no connected channels)";
    edm::LogWarning(mlDqmClient_) << sss.str();
  }
  
  // Some debug 
  std::stringstream ss;
  ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
     << " Dump of enabled FED channels:" 
     << std::endl;
  std::map< uint16_t, std::vector<uint16_t> >::const_iterator fed = enabled.begin();
  for ( ; fed != enabled.end(); fed++ ) {
    ss << " Enabled " << fed->second.size()
       << " channels for FED id " 
       << std::setw(3) << fed->first << ": ";
    if ( !fed->second.empty() ) { 
      uint16_t first = fed->second.front();
      uint16_t last = fed->second.front();
      std::vector<uint16_t>::const_iterator chan = fed->second.begin();
      for ( ; chan != fed->second.end(); chan++ ) { 
	if ( chan != fed->second.begin() ) {
	  if ( *chan != last+1 ) { 
	    ss << std::setw(2) << first << "->" << std::setw(2) << last << ", ";
	    if ( chan != fed->second.end() ) { first = *(chan+1); }
	  } 
	}
	last = *chan;
      }
      if ( first != last ) { ss << std::setw(2) << first << "->" << std::setw(2) << last; }
      ss << std::endl;
    }
  }
  LogTrace(mlDqmClient_) << ss.str();
  
}

// -----------------------------------------------------------------------------
// 
void FastFedCablingHistosUsingDb::addDcuDetIds() {

  if ( !cabling() ) {
    edm::LogError(mlDqmClient_) 
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripFedCabling object!";
    return;
  }

  // retrieve descriptions for dcu id and det id 
  SiStripConfigDb::DeviceDescriptions dcus = db()->getDeviceDescriptions( DCU ); 
  SiStripConfigDb::DcuDetIdMap detids = db()->getDcuDetIdMap(); 
  
  if ( dcus.empty() ) { 
    edm::LogError(mlCabling_)
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " No DCU descriptions found!";
    return;
  }
  
  if ( detids.empty() ) { 
    edm::LogWarning(mlCabling_)
      << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
      << " DCU-DetId map is empty!";
  }
  
  Analyses::iterator ianal = data().begin();
  Analyses::iterator janal = data().end();
  for ( ; ianal != janal; ++ianal ) { 

    // check if dcu id and lld channel have been identified 
    if ( !ianal->second->isValid() ) { continue; }
    
    // retrieve analysis object
    FastFedCablingAnalysis* anal = dynamic_cast<FastFedCablingAnalysis*>( ianal->second );
    
    if ( !anal ) {
      edm::LogError(mlDqmClient_) 
	<< "[FastFedCablingHistosUsingDb::" << __func__ << "]"
	<< " NULL pointer to FastFedCablingAnalysis object!";
      return;
    }
    
    // find dcu that matches analysis result 
    bool found = false;
    SiStripConfigDb::DeviceDescriptions::const_iterator idcu = dcus.begin();
    SiStripConfigDb::DeviceDescriptions::const_iterator jdcu = dcus.end();
    while ( !found && idcu != jdcu ) {
      dcuDescription* dcu = dynamic_cast<dcuDescription*>( *idcu );
      if ( dcu ) { 
	if ( dcu->getDcuType() == "FEH" ) { 
	  if ( dcu->getDcuHardId() == anal->dcuHardId() ) {
	    found = true; 
	    anal->dcuId( dcu->getDcuHardId() ); 
	    const SiStripConfigDb::DeviceAddress& addr = db()->deviceAddress(*dcu);
	    uint32_t fec_key = SiStripFecKey( addr.fecCrate_,
					      addr.fecSlot_,
					      addr.fecRing_,
					      addr.ccuAddr_,
					      addr.ccuChan_,
					      anal->lldCh() ).key();
	    anal->fecKey( fec_key );
	    SiStripConfigDb::DcuDetIdMap::const_iterator idet = detids.find( dcu->getDcuHardId() );
	    if ( idet != detids.end() ) { anal->detId( idet->second->getDetId() ); }
	  }
	}
      }
      idcu++;
    }

  }

}

// -----------------------------------------------------------------------------
/** */
void FastFedCablingHistosUsingDb::create( SiStripConfigDb::AnalysisDescriptions& desc,
					  Analysis analysis ) {

#ifdef USING_NEW_DATABASE_MODEL
  
  FastFedCablingAnalysis* anal = dynamic_cast<FastFedCablingAnalysis*>( analysis->second );
  if ( !anal ) { return; }

  if ( !anal->isValid() ) { return; } //@@ only store valid descriptions for now!
  
  SiStripFecKey fec_key( anal->fecKey() );
  SiStripFedKey fed_key( anal->fedKey() );
    
  for ( uint16_t iapv = 0; iapv < 2; ++iapv ) {
    
    // Create description
    FastFedCablingAnalysisDescription* tmp;
    tmp = new FastFedCablingAnalysisDescription( anal->highLevel(),
						 anal->highRms(),
						 anal->lowLevel(),
						 anal->lowRms(),
						 anal->max(),
						 anal->min(),
						 anal->dcuId(),
						 anal->lldCh(),
						 anal->isDirty(),
						 FastFedCablingAnalysis::threshold_,
						 FastFedCablingAnalysis::dirtyThreshold_,
						 fec_key.fecCrate(),
						 fec_key.fecSlot(),
						 fec_key.fecRing(),
						 fec_key.ccuAddr(),
						 fec_key.ccuChan(),
						 SiStripFecKey::i2cAddr( fec_key.lldChan(), !iapv ), 
						 db()->dbParams().partition_,
						 db()->dbParams().runNumber_,
						 anal->isValid(),
						 "",
						 fed_key.fedId(),
						 fed_key.feUnit(),
						 fed_key.feChan(),
						 fed_key.fedApv() );
      
    // Add comments
    typedef std::vector<std::string> Strings;
    Strings errors = anal->getErrorCodes();
    Strings::const_iterator istr = errors.begin();
    Strings::const_iterator jstr = errors.end();
    for ( ; istr != jstr; ++istr ) { tmp->addComments( *istr ); }

    // Store description
    desc.push_back( tmp );
      
  }

#endif
  
}

// -----------------------------------------------------------------------------
// prints debug info on good, dirty, missing connections, and missing devices
void FastFedCablingHistosUsingDb::connections( const SiStripConfigDb::DeviceDescriptions& dcus, 
					       const SiStripConfigDb::DcuDetIdMap& detids ) {
  
  // strings
  std::vector<std::string> valid;
  std::vector<std::string> dirty;
  std::vector<std::string> trimdac;
  std::vector<std::string> missing;
  std::vector<std::string> devices;
  uint32_t missing_pairs = 0;

  // iterate through analyses
  std::vector<uint32_t> found_dcus;
  Analyses::iterator ianal = data().begin();
  Analyses::iterator janal = data().end();
  for ( ; ianal != janal; ++ianal ) {

    // extract fast fed cabling object
    FastFedCablingAnalysis* anal = dynamic_cast<FastFedCablingAnalysis*>( ianal->second );
    if ( !anal ) { continue; }
    
    // construct string for valid and dirty connections
    std::stringstream ss;
    SiStripFecKey( anal->fecKey() ).terse(ss); ss << " "; 
    SiStripFedKey( anal->fedKey() ).terse(ss); ss << " ";
    ss << "DcuId=" << std::hex << std::setw(8) << std::setfill('0') << anal->dcuId() << std::dec << " ";
    ss << "DetId=" << std::hex << std::setw(8) << std::setfill('0') << anal->detId() << std::dec;
    if ( anal->isValid() ) { valid.push_back( ss.str() ); }
    if ( anal->isDirty() ) { dirty.push_back( ss.str() ); }
    if ( anal->badTrimDac() ) { trimdac.push_back( ss.str() ); }
    
    // record "found" dcus
    found_dcus.push_back( anal->dcuHardId() );

  }
  
  // iterate through dcu devices
  SiStripConfigDb::DeviceDescriptions::const_iterator idcu = dcus.begin();
  SiStripConfigDb::DeviceDescriptions::const_iterator jdcu = dcus.end();
  for ( ; idcu != jdcu; ++idcu ) {
    
    // extract dcu description
    dcuDescription* dcu = dynamic_cast<dcuDescription*>( *idcu );
    if ( !dcu ) { continue; }
    if ( dcu->getDcuType() != "FEH" ) { continue; }
    SiStripConfigDb::DeviceAddress dcu_addr = db()->deviceAddress( *dcu );
    
    // continue if dcu has been "found"
    std::vector<uint32_t>::const_iterator iter = find ( found_dcus.begin(), found_dcus.end(), dcu->getDcuHardId() );
    if ( iter != found_dcus.end() ) { continue; }
    
    // find detid for "missing" dcu
    SiStripConfigDb::DcuDetIdMap::const_iterator idet = detids.find( dcu->getDcuHardId() );
    if ( idet == detids.end() ) { continue; }
    if ( idet->second ) { continue; }
    
    // retrieve number of apv pairs
    uint16_t npairs = idet->second->getApvNumber()/2; 
  
    // retrieve apvs for given dcu
    vector<bool> addrs; 
    addrs.resize(6,false);
    SiStripConfigDb::DeviceDescriptions apvs = db()->getDeviceDescriptions( APV25 );
    SiStripConfigDb::DeviceDescriptions::const_iterator iapv = apvs.begin();
    SiStripConfigDb::DeviceDescriptions::const_iterator japv = apvs.end();
    for ( ; iapv != japv; ++iapv ) {
      apvDescription* apv = dynamic_cast<apvDescription*>( *iapv );
      if ( !apv ) { continue; }
      SiStripConfigDb::DeviceAddress apv_addr = db()->deviceAddress( *apv );
      if ( apv_addr.fecCrate_ == dcu_addr.fecCrate_ &&
	   apv_addr.fecSlot_  == dcu_addr.fecSlot_ &&
	   apv_addr.fecRing_  == dcu_addr.fecRing_ &&
	   apv_addr.ccuAddr_  == dcu_addr.ccuAddr_ &&
	   apv_addr.ccuChan_  == dcu_addr.ccuChan_ ) {
	uint16_t pos = apv_addr.i2cAddr_ - 32;
	if ( pos < 6 ) { addrs[pos] = true; }
      }
    }
    
    // construct strings for missing fibres
    uint16_t pairs = 0;
    if ( addrs[0] || addrs[1] ) {
      pairs++;
      std::stringstream ss;
      SiStripFecKey( dcu_addr.fecCrate_,
		     dcu_addr.fecSlot_,
		     dcu_addr.fecRing_,
		     dcu_addr.ccuAddr_,
		     dcu_addr.ccuChan_,
		     1 ).terse(ss); 
      ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
      ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
      missing.push_back( ss.str() ); 
    } 
    if ( addrs[2] || addrs[3] ) {
      pairs++;
      std::stringstream ss;
      SiStripFecKey( dcu_addr.fecCrate_,
		     dcu_addr.fecSlot_,
		     dcu_addr.fecRing_,
		     dcu_addr.ccuAddr_,
		     dcu_addr.ccuChan_,
		     2 ).terse(ss); 
      ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
      ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
      missing.push_back( ss.str() ); 
    } 
    if ( addrs[4] || addrs[5] ) {
      pairs++;
      std::stringstream ss;
      SiStripFecKey( dcu_addr.fecCrate_,
		     dcu_addr.fecSlot_,
		     dcu_addr.fecRing_,
		     dcu_addr.ccuAddr_,
		     dcu_addr.ccuChan_,
		     3 ).terse(ss); 
      ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
      ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
      missing.push_back( ss.str() ); 
    }
    
    if ( pairs != npairs ) {
      
      missing_pairs = npairs - pairs;
      
      if ( !addrs[0] ) { 
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       1, 32 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
      if ( !addrs[1] ) { 
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       1, 33 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
      if ( !addrs[2] && npairs == 3 ) {
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       2, 34 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
      if ( !addrs[3] && npairs == 3 ) {
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       2, 35 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
      if ( !addrs[4] ) { 
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       3, 36 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
      if ( !addrs[5] ) { 
	std::stringstream ss;
	SiStripFecKey( dcu_addr.fecCrate_,
		       dcu_addr.fecSlot_,
		       dcu_addr.fecRing_,
		       dcu_addr.ccuAddr_,
		       dcu_addr.ccuChan_,
		       3, 37 ).terse(ss); 
	ss << " DcuId=" << std::hex << std::setw(8) << std::setfill('0') << dcu->getDcuHardId() << std::dec;
	ss << " DetId=" << std::hex << std::setw(8) << std::setfill('0') << idet->first << std::dec;
	devices.push_back( ss.str() ); 
      }
      
    }

  }

  // summary
  { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " Summary of connections: " << std::endl
       << " \"Good\" connections     : " << valid.size() << std::endl
       << " \"Dirty\" connections    : " << dirty.size() << std::endl
       << " \"Bad\" TrimDAQ settings : " << trimdac.size() << std::endl
       << " \"Missing\" connections  : " << missing.size() << std::endl
       << " \"Missing\" APV pairs    : " << missing_pairs << std::endl
       << " \"Missing\" APVs         : " << devices.size() << std::endl;
    edm::LogVerbatim(mlCabling_) << ss.str();
  }

  // good connections
  if ( !valid.empty() ) { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " List of \"good\" connections: " << std::endl;
    std::vector<std::string>::const_iterator istr = valid.begin();
    std::vector<std::string>::const_iterator jstr = valid.end();
    for ( ; istr != jstr; ++istr ) { ss << *istr << std::endl; }
    LogTrace(mlCabling_) << ss.str();
  }

  // dirty connections
  if ( !dirty.empty() ) { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " List of \"dirty\" connections: " << std::endl;
    std::vector<std::string>::const_iterator istr = dirty.begin();
    std::vector<std::string>::const_iterator jstr = dirty.end();
    for ( ; istr != jstr; ++istr ) { ss << *istr << std::endl; }
    edm::LogWarning(mlCabling_) << ss.str();
  }

  // TrimDAC connections
  if ( !trimdac.empty() ) { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " List of \"bad\" TrimDAC settings: " << std::endl;
    std::vector<std::string>::const_iterator istr = trimdac.begin();
    std::vector<std::string>::const_iterator jstr = trimdac.end();
    for ( ; istr != jstr; ++istr ) { ss << *istr << std::endl; }
    edm::LogWarning(mlCabling_) << ss.str();
  }

  // missing connections
  if ( !missing.empty() ) { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " List of \"missing\" connections: " << std::endl;
    std::vector<std::string>::const_iterator istr = missing.begin();
    std::vector<std::string>::const_iterator jstr = missing.end();
    for ( ; istr != jstr; ++istr ) { ss << *istr << std::endl; }
    edm::LogError(mlCabling_) << ss.str();
  }

  // missing devices
  if ( !devices.empty() ) { 
    std::stringstream ss;
    ss << "[FastFedCablingHistosUsingDb::" << __func__ << "]"
       << " List of \"missing\" APVs: " << std::endl;
    std::vector<std::string>::const_iterator istr = devices.begin();
    std::vector<std::string>::const_iterator jstr = devices.end();
    for ( ; istr != jstr; ++istr ) { ss << *istr << std::endl; }
    edm::LogWarning(mlCabling_) << ss.str();
  }

}
