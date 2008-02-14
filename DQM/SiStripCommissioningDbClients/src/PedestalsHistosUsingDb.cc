// Last commit: $Id: PedestalsHistosUsingDb.cc,v 1.10 2008/02/07 17:02:58 bainbrid Exp $

#include "DQM/SiStripCommissioningDbClients/interface/PedestalsHistosUsingDb.h"
#include "CondFormats/SiStripObjects/interface/PedestalsAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include <iostream>

using namespace sistrip;

// -----------------------------------------------------------------------------
/** */
PedestalsHistosUsingDb::PedestalsHistosUsingDb( MonitorUserInterface* mui,
						const DbParams& params )
  : CommissioningHistosUsingDb( params ),
    PedestalsHistograms( mui )
{
  LogTrace(mlDqmClient_) 
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
PedestalsHistosUsingDb::PedestalsHistosUsingDb( MonitorUserInterface* mui,
						SiStripConfigDb* const db )
  : CommissioningHistograms( mui, sistrip::PEDESTALS ),
    CommissioningHistosUsingDb( db, mui, sistrip::PEDESTALS ),
    PedestalsHistograms( mui )
{
  LogTrace(mlDqmClient_) 
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
PedestalsHistosUsingDb::PedestalsHistosUsingDb( DaqMonitorBEInterface* bei,
						SiStripConfigDb* const db ) 
  : CommissioningHistosUsingDb( db, sistrip::PEDESTALS ),
    PedestalsHistograms( bei )
{
  LogTrace(mlDqmClient_) 
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
/** */
PedestalsHistosUsingDb::~PedestalsHistosUsingDb() {
  LogTrace(mlDqmClient_) 
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
/** */
void PedestalsHistosUsingDb::uploadConfigurations() {

  if ( !db() ) {
    edm::LogError(mlDqmClient_) 
      << "[PedestalsHistosUsingDb::" << __func__ << "]"
      << " NULL pointer to SiStripConfigDb interface!"
      << " Aborting upload...";
    return;
  }
  
  // Update FED descriptions with new peds/noise values
  const SiStripConfigDb::FedDescriptions& feds = db()->getFedDescriptions(); 
  update( const_cast<SiStripConfigDb::FedDescriptions&>(feds) );
  if ( doUploadConf() ) { 
    edm::LogVerbatim(mlDqmClient_) 
      << "[PedestalsHistosUsingDb::" << __func__ << "]"
      << " Uploading pedestals/noise to DB...";
    db()->uploadFedDescriptions(true); // always major version 
    edm::LogVerbatim(mlDqmClient_) 
      << "[PedestalsHistosUsingDb::" << __func__ << "]"
      << " Completed database upload of " << feds.size() 
      << " FED descriptions!";
  } else {
    edm::LogWarning(mlDqmClient_) 
      << "[PedestalsHistosUsingDb::" << __func__ << "]"
      << " TEST only! No pedestals/noise values will be uploaded to DB...";
  }
  
}

// -----------------------------------------------------------------------------
/** */
void PedestalsHistosUsingDb::update( SiStripConfigDb::FedDescriptions& feds ) {
 
  // Iterate through feds and update fed descriptions
  uint16_t updated = 0;
  SiStripConfigDb::FedDescriptions::iterator ifed;
  for ( ifed = feds.begin(); ifed != feds.end(); ifed++ ) {
    
    for ( uint16_t ichan = 0; ichan < sistrip::FEDCH_PER_FED; ichan++ ) {

      // Build FED and FEC keys
      const FedChannelConnection& conn = cabling()->connection( (*ifed)->getFedId(), ichan );
      if ( conn.fecCrate() == sistrip::invalid_ ||
	   conn.fecSlot() == sistrip::invalid_ ||
	   conn.fecRing() == sistrip::invalid_ ||
	   conn.ccuAddr() == sistrip::invalid_ ||
	   conn.ccuChan() == sistrip::invalid_ ||
	   conn.lldChannel() == sistrip::invalid_ ) { continue; }
      SiStripFedKey fed_key( conn.fedId(), 
			     SiStripFedKey::feUnit( conn.fedCh() ),
			     SiStripFedKey::feChan( conn.fedCh() ) );
      SiStripFecKey fec_key( conn.fecCrate(), 
			     conn.fecSlot(), 
			     conn.fecRing(), 
			     conn.ccuAddr(), 
			     conn.ccuChan(), 
			     conn.lldChannel() );

      // Locate appropriate analysis object 
      Analyses::const_iterator iter = data().find( fec_key.key() );
      if ( iter != data().end() ) {

	// Check if analysis is valid
	if ( !iter->second->isValid() ) { continue; }

	PedestalsAnalysis* anal = dynamic_cast<PedestalsAnalysis*>( iter->second );
	if ( !anal ) { continue; }
	
	// Iterate through APVs and strips
	for ( uint16_t iapv = 0; iapv < sistrip::APVS_PER_FEDCH; iapv++ ) {
	  for ( uint16_t istr = 0; istr < anal->peds()[iapv].size(); istr++ ) { 
	    
	    static float high_threshold = 5.;
	    static float low_threshold  = 2.;
	    static bool  disable_strip  = false;
	    Fed9U::Fed9UStripDescription data( static_cast<uint32_t>( anal->peds()[iapv][istr] ), 
					       high_threshold, 
					       low_threshold, 
					       anal->noise()[iapv][istr],
					       disable_strip );
	    Fed9U::Fed9UAddress addr( ichan, iapv, istr );
	    (*ifed)->getFedStrips().setStrip( addr, data );
	    
	  }
	}
	updated++;
      
      } else {
	edm::LogWarning(mlDqmClient_) 
	  << "[PedestalsHistosUsingDb::" << __func__ << "]"
	  << " Unable to find pedestals/noise for FedKey/Id/Ch: " 
	  << hex << setw(8) << setfill('0') << fed_key.key() << dec << "/"
	  << (*ifed)->getFedId() << "/"
	  << ichan
	  << " and device with FEC/slot/ring/CCU/LLD " 
	  << fec_key.fecCrate() << "/"
	  << fec_key.fecSlot() << "/"
	  << fec_key.fecRing() << "/"
	  << fec_key.ccuAddr() << "/"
	  << fec_key.ccuChan() << "/"
	  << fec_key.channel();
      }
    }
  }

  edm::LogVerbatim(mlDqmClient_) 
    << "[PedestalsHistosUsingDb::" << __func__ << "]"
    << " Updated FED pedestals/noise for " 
    << updated << " channels";

}

// -----------------------------------------------------------------------------
/** */
void PedestalsHistosUsingDb::create( SiStripConfigDb::AnalysisDescriptions& desc,
					  Analysis analysis ) {

  PedestalsAnalysis* anal = dynamic_cast<PedestalsAnalysis*>( analysis->second );
  if ( !anal ) { return; }
  
  SiStripFecKey fec_key( anal->fecKey() );
  SiStripFedKey fed_key( anal->fedKey() );
  
  for ( uint16_t iapv = 0; iapv < 2; ++iapv ) {
    
    // Create description
    PedestalsAnalysisDescription* tmp;
    tmp = new PedestalsAnalysisDescription( anal->dead()[iapv],
					    anal->noisy()[iapv],
					    anal->pedsMean()[iapv],
					    anal->pedsSpread()[iapv],
					    anal->noiseMean()[iapv],
					    anal->noiseSpread()[iapv],
					    anal->rawMean()[iapv],
					    anal->rawSpread()[iapv],
					    anal->pedsMax()[iapv], 
					    anal->pedsMin()[iapv], 
					    anal->noiseMax()[iapv],
					    anal->noiseMin()[iapv],
					    anal->rawMax()[iapv],
					    anal->rawMin()[iapv],
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

}

