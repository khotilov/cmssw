// Last commit: $Id: test_SiStripFecKey.cc,v 1.5 2008/01/15 16:27:57 bainbrid Exp $

#include "DataFormats/SiStripCommon/test/plugins/test_SiStripFecKey.h"
#include "FWCore/Framework/interface/Event.h" 
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/Constants.h" 
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include <algorithm>
//#include <functional>

using namespace sistrip;

// -----------------------------------------------------------------------------
// 
testSiStripFecKey::testSiStripFecKey( const edm::ParameterSet& pset ) 
{
  LogTrace(mlDqmCommon_)
    << "[testSiStripFecKey::" << __func__ << "]"
    << " Constructing object...";
}

// -----------------------------------------------------------------------------
// 
testSiStripFecKey::~testSiStripFecKey() {
  LogTrace(mlDqmCommon_)
    << "[testSiStripFecKey::" << __func__ << "]"
    << " Destructing object...";
}

// -----------------------------------------------------------------------------
// 
void testSiStripFecKey::beginJob( const edm::EventSetup& setup ) {
  
  uint32_t cntr = 0;
  uint32_t start = time(NULL);

  edm::LogInfo(mlDqmCommon_)
    << "[SiStripFecKey::" << __func__ << "]"
    << " Tests the generation of keys...";

  std::vector<uint32_t> keys;
  
  // FEC crates
  for ( uint16_t icrate = 0; icrate <= sistrip::FEC_CRATE_MAX+1; icrate++ ) {
    if ( icrate > 1 && icrate < sistrip::FEC_CRATE_MAX ) { continue; }

    // FEC slots
    for ( uint16_t ifec = 0; ifec <= sistrip::SLOTS_PER_CRATE+1; ifec++ ) {
      if ( ifec > 1 && ifec < sistrip::SLOTS_PER_CRATE ) { continue; }

      // FEC rings
      for ( uint16_t iring = 0; iring <= sistrip::FEC_RING_MAX+1; iring++ ) {
	if ( iring > 1 && iring < sistrip::FEC_RING_MAX ) { continue; }

	// CCU addr
	for ( uint16_t iccu = 0; iccu <= sistrip::CCU_ADDR_MAX+1; iccu++ ) {
	  if ( iccu > 1 && iccu < sistrip::CCU_ADDR_MAX ) { continue; }
	  
	  // CCU channel
	  for ( uint16_t ichan = 0; ichan <= sistrip::CCU_CHAN_MAX+1; ichan++ ) {
 	    if ( ichan > 1 && 
 		 ichan != sistrip::CCU_CHAN_MIN &&
 		 ichan < sistrip::CCU_CHAN_MAX-1 ) { continue; }
	    
	    // LLD channels
	    for ( uint16_t illd = 0; illd <= sistrip::LLD_CHAN_MAX+1; illd++ ) {
	      if ( illd > 1 && illd < sistrip::LLD_CHAN_MAX ) { continue; }
	      
	      // APV
	      for ( uint16_t iapv = 0; iapv <= sistrip::APV_I2C_MAX+1; iapv++ ) {
		if ( iapv > 1 && 
		     iapv != sistrip::APV_I2C_MIN &&
		     iapv < sistrip::APV_I2C_MAX ) { continue; }
		
		// Some debug
		if ( !(cntr%1000) ) {
		  LogTrace(mlDqmCommon_)
		    << "[SiStripFecKey::" << __func__ << "]"
		    << " Cntr: " << cntr;
		}
		cntr++;
		
		// Print out FEC
		std::stringstream ss;
		ss << std::endl
		   << "[SiStripFecKey::" << __func__ << "]"
		   << " crate/FEC/ring/CCU/module/LLD/I2C: "
		   << icrate << "/"
		   << ifec << "/"
		   << iring << "/"
		   << iccu << "/"
		   << ichan << "/"
		   << illd << "/"
		   << iapv << std::endl << std::endl;
		
		SiStripFecKey tmp1( icrate, ifec, iring, iccu, ichan, illd, iapv );
		SiStripFecKey tmp2 = SiStripFecKey( tmp1.key() );
		SiStripFecKey tmp3 = SiStripFecKey( tmp1.path() );
		SiStripFecKey tmp4 = SiStripFecKey( tmp1 );
		SiStripFecKey tmp5; tmp5 = tmp1;

		keys.push_back(tmp1.key());
		
		ss << ">>> original       : "; tmp1.terse(ss); ss << std::endl;
		ss << ">>> from FEC key   : "; tmp1.terse(ss); ss << std::endl;
		ss << ">>> from directory : "; tmp1.terse(ss); ss << std::endl;
		ss << ">>> directory      : " << tmp1.path() << std::endl;
		ss << ">>> isValid        : " << tmp1.isValid()
		   << " " << tmp1.isValid( tmp1.granularity() )
		   << " " << tmp1.isValid( sistrip::APV ) << std::endl
		   << ">>> isInvalid      : " << tmp1.isInvalid()
		   << " " << tmp1.isInvalid( tmp1.granularity() )
		   << " " << tmp1.isInvalid( sistrip::APV );

// 		ss << ">>> original:" << std::endl << tmp1 << std::endl
// 		   << ">>> from FEC key:" << std::endl << tmp2 << std::endl
// 		   << ">>> from directory:" << std::endl << tmp3 << std::endl
// 		   << ">>> isValid:   " << tmp1.isValid()
// 		   << " " << tmp1.isValid( tmp1.granularity() )
// 		   << " " << tmp1.isValid( sistrip::APV ) << std::endl
// 		   << ">>> isInvalid: " << tmp1.isInvalid()
// 		   << " " << tmp1.isInvalid( tmp1.granularity() )
// 		   << " " << tmp1.isInvalid( sistrip::APV );

		LogTrace(mlDqmCommon_) << ss.str();
		
	      }
	    }
	  }
	}
      }
    }
  }

  std::sort( keys.begin(), keys.end() );
  typedef std::vector<uint32_t>::iterator iter;
  SiStripFecKey value( static_cast<uint16_t>(4),
		       static_cast<uint16_t>(21),
		       static_cast<uint16_t>(8),
		       static_cast<uint16_t>(127) );
  std::pair<iter,iter> temp = 
    std::equal_range( keys.begin(), 
 		      keys.end(),
		      value.key(),
 		      ConsistentWithKey(value) );
  edm::LogVerbatim(mlDqmCommon_)
    << "[SiStripFecKey::" << __func__ << "]"
    << " number of keys = " << keys.size()
    << " number of matching = " << temp.second - temp.first;

  if ( temp.first != temp.second ) {
    std::stringstream ss;
    ss << std::endl;
    for ( iter ii = temp.first; ii != temp.second; ++ii ) {
      SiStripFecKey(*ii).terse(ss); ss << std::endl;
    }
    LogTrace(mlDqmCommon_)
      << "[SiStripFecKey::" << __func__ << "] begin"
      << ss.str()
      << "[SiStripFecKey::" << __func__ << "] end";
  }

  if ( find( keys.begin(), keys.end(), value.key() ) != keys.end() ) {
    edm::LogVerbatim(mlDqmCommon_)
      << "[SiStripFecKey::" << __func__ << "]"
      << " found!!! ";
  }
  
  edm::LogVerbatim(mlDqmCommon_)
    << "[SiStripFecKey::" << __func__ << "]"
    << " Processed " << cntr
    << " FecKeys in " << (time(NULL)-start)
    << " seconds at an average rate of " << (cntr*1.) / ((time(NULL)-start)*1.)
    << " per second...";

  // Tests for utility methods

  SiStripFecKey invalid;
  SiStripFecKey inv(sistrip::invalid_,
		    sistrip::invalid_,
		    sistrip::invalid_,
		    sistrip::invalid_,
		    sistrip::invalid_,
		    sistrip::invalid_,
		    sistrip::invalid_);
  SiStripFecKey valid(1,2,1,1,16,1,32);
  SiStripFecKey all(0,0,0,0,0,0,0);
  SiStripFecKey same(valid);
  SiStripFecKey equal = valid;
  SiStripFecKey equals; 
  equals = valid;
  SiStripFecKey to_gran(valid,sistrip::CCU_CHAN); 

  std::stringstream ss;

  ss << "[SiStripFecKey::" << __func__ << "]"
     << " Tests for utility methods..." << std::endl;

  ss << ">>>> invalid.path: " << invalid << std::endl
     << ">>>> inv.path:     " << inv << std::endl
     << ">>>> valid.path:   " << valid << std::endl
     << ">>>> all.path:     " << all << std::endl
     << ">>>> same.path:    " << same << std::endl
     << ">>>> equal.path:   " << equal << std::endl
     << ">>>> equals.path:  " << equals << std::endl
     << ">>>> to_gran.path:  " << to_gran << std::endl;
  
  ss << std::hex
     << ">>>> invalid.key:  " << invalid.key() << std::endl
     << ">>>> valid.key:    " << valid.key() << std::endl
     << ">>>> all.key:      " << all.key() << std::endl
     << std::dec;
  
  ss << ">>>> invalid.isInvalid: " << invalid.isInvalid() << std::endl
     << ">>>> invalid.isValid:   " << invalid.isValid() << std::endl
     << ">>>> valid.isInvalid:   " << valid.isInvalid() << std::endl
     << ">>>> valid.isValid:     " << valid.isValid() << std::endl
     << ">>>> all.isInvalid:     " << all.isInvalid() << std::endl
     << ">>>> all.isValid:       " << all.isValid() << std::endl;

  ss << ">>>> valid.isEqual(valid):        " << valid.isEqual(valid) << std::endl
     << ">>>> valid.isConsistent(valid):   " << valid.isConsistent(valid) << std::endl
     << ">>>> valid.isEqual(invalid):      " << valid.isEqual(invalid) << std::endl
     << ">>>> valid.isConsistent(invalid): " << valid.isConsistent(invalid) << std::endl
     << ">>>> valid.isEqual(all):          " << valid.isEqual(all) << std::endl
     << ">>>> valid.isConsistent(all):     " << valid.isConsistent(all) << std::endl
     << ">>>> valid.isEqual(same):         " << valid.isEqual(same) << std::endl
     << ">>>> valid.isEqual(equal):        " << valid.isEqual(equal) << std::endl
     << ">>>> valid.isEqual(equals):       " << valid.isEqual(equals) << std::endl;
  
  LogTrace(mlDqmCommon_) << ss.str();
  
}

// -----------------------------------------------------------------------------
// 
void testSiStripFecKey::analyze( const edm::Event& event, 
				  const edm::EventSetup& setup ) {
  LogTrace(mlDqmCommon_) 
    << "[SiStripFecKey::" << __func__ << "]"
    << " Analyzing run/event "
    << event.id().run() << "/"
    << event.id().event();
}


