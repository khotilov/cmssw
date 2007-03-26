// Last commit: $Id: SiStripFecKey.cc,v 1.4 2007/03/21 08:22:59 bainbrid Exp $

#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForHardwareSystems.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForDqm.h"
#include "DataFormats/SiStripCommon/interface/ConstantsForView.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include <iomanip>
#include <sstream>

// -----------------------------------------------------------------------------
//
SiStripFecKey::SiStripFecKey( const uint16_t& fec_crate, 
			      const uint16_t& fec_slot, 
			      const uint16_t& fec_ring, 
			      const uint16_t& ccu_addr, 
			      const uint16_t& ccu_chan,
			      const uint16_t& lld_chan,
			      const uint16_t& i2c_addr ) :
  SiStripKey(),
  fecCrate_(fec_crate), 
  fecSlot_(fec_slot),
  fecRing_(fec_ring), 
  ccuAddr_(ccu_addr),
  ccuChan_(ccu_chan),
  lldChan_(lld_chan),
  i2cAddr_(i2c_addr)
{
  // order is important!
  initFromValue();
  initFromKey();
  initFromPath();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripFecKey::SiStripFecKey( const uint32_t& fec_key ) :
  SiStripKey(fec_key),
  fecCrate_(sistrip::invalid_), 
  fecSlot_(sistrip::invalid_),
  fecRing_(sistrip::invalid_), 
  ccuAddr_(sistrip::invalid_),
  ccuChan_(sistrip::invalid_), 
  lldChan_(sistrip::invalid_),
  i2cAddr_(sistrip::invalid_)
{
  // order is important!
  initFromKey(); 
  initFromValue();
  initFromPath();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripFecKey::SiStripFecKey( const std::string& path ) :
  SiStripKey(path),
  fecCrate_(sistrip::invalid_), 
  fecSlot_(sistrip::invalid_),
  fecRing_(sistrip::invalid_), 
  ccuAddr_(sistrip::invalid_),
  ccuChan_(sistrip::invalid_), 
  lldChan_(sistrip::invalid_),
  i2cAddr_(sistrip::invalid_)
{
  // order is important!
  initFromPath();
  initFromValue();
  initFromKey(); 
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripFecKey::SiStripFecKey( const SiStripFecKey& input ) :
  SiStripKey(),
  fecCrate_(input.fecCrate()), 
  fecSlot_(input.fecSlot()),
  fecRing_(input.fecRing()), 
  ccuAddr_(input.ccuAddr()),
  ccuChan_(input.ccuChan()), 
  lldChan_(input.lldChan()), 
  i2cAddr_(input.i2cAddr())
{
  key(input.key());
  path(input.path());
  granularity(input.granularity());
}

// -----------------------------------------------------------------------------
// 
SiStripFecKey::SiStripFecKey( const SiStripKey& input ) :
  SiStripKey(),
  fecCrate_(sistrip::invalid_), 
  fecSlot_(sistrip::invalid_),
  fecRing_(sistrip::invalid_), 
  ccuAddr_(sistrip::invalid_),
  ccuChan_(sistrip::invalid_), 
  lldChan_(sistrip::invalid_),
  i2cAddr_(sistrip::invalid_)
{
  SiStripKey& temp = const_cast<SiStripKey&>(input);
  SiStripFecKey& fec_key = dynamic_cast<SiStripFecKey&>(temp);
  if ( (&fec_key) ) {
    key(fec_key.key());
    path(fec_key.path());
    granularity(fec_key.granularity());
    fecCrate_ = fec_key.fecCrate(); 
    fecSlot_ = fec_key.fecSlot();
    fecRing_ = fec_key.fecRing(); 
    ccuAddr_ = fec_key.ccuAddr();
    ccuChan_ = fec_key.ccuChan(); 
    lldChan_ = fec_key.lldChan();
    i2cAddr_ = fec_key.i2cAddr();
  }
}

// -----------------------------------------------------------------------------
// 
SiStripFecKey::SiStripFecKey() :
  SiStripKey(),
  fecCrate_(sistrip::invalid_), 
  fecSlot_(sistrip::invalid_),
  fecRing_(sistrip::invalid_), 
  ccuAddr_(sistrip::invalid_),
  ccuChan_(sistrip::invalid_), 
  lldChan_(sistrip::invalid_),
  i2cAddr_(sistrip::invalid_)
{;}

// -----------------------------------------------------------------------------
// 
uint16_t SiStripFecKey::hybridPos( const uint16_t& i2c_addr ) {
  if ( i2c_addr < sistrip::APV_I2C_MIN ||
       i2c_addr > sistrip::APV_I2C_MAX ) {
    return sistrip::invalid_;
  }
  return ( i2c_addr - sistrip::APV_I2C_MIN + 1 );
}

// -----------------------------------------------------------------------------
// 
uint16_t SiStripFecKey::i2cAddr( const uint16_t& hybrid_pos ) {
  if ( !hybrid_pos ||
       hybrid_pos > 
       ( sistrip::APV_I2C_MAX - 
	 sistrip::APV_I2C_MIN + 1 ) ) {
    return sistrip::invalid_;
  }
  return ( hybrid_pos + sistrip::APV_I2C_MIN - 1 );
}

// -----------------------------------------------------------------------------
// 
uint16_t SiStripFecKey::i2cAddr( const uint16_t& lld_chan,
				 const bool& first_apv ) {
  if ( lld_chan < sistrip::LLD_CHAN_MIN ||
       lld_chan > sistrip::LLD_CHAN_MAX ) {
    return sistrip::invalid_; 
  }
  return ( lld_chan * sistrip::APVS_PER_CHAN + (first_apv?1:2) );
}

// -----------------------------------------------------------------------------
// 
uint16_t SiStripFecKey::lldChan( const uint16_t& i2c_addr ) {
  if ( i2c_addr < sistrip::APV_I2C_MIN ||
       i2c_addr > sistrip::APV_I2C_MAX ) {
    return sistrip::invalid_;
  }
  return ( ( i2c_addr - sistrip::APV_I2C_MIN ) / 2 + 1 );
}

// -----------------------------------------------------------------------------
// 
bool SiStripFecKey::firstApvOfPair( const uint16_t& i2c_addr ) {
  if ( i2c_addr < sistrip::APV_I2C_MIN ||
       i2c_addr > sistrip::APV_I2C_MAX ) {
    return sistrip::invalid_;
  }
  return ( ( ( i2c_addr - sistrip::APV_I2C_MIN ) % 2 ) == 0 );
}

// -----------------------------------------------------------------------------
// 
bool SiStripFecKey::isEqual( const SiStripKey& key ) const {
  SiStripKey& temp = const_cast<SiStripKey&>(key);
  SiStripFecKey& input = dynamic_cast<SiStripFecKey&>(temp);
  if ( !(&input) ) { return false; }
  if ( fecCrate_ == input.fecCrate() &&
       fecSlot_ == input.fecSlot() &&
       fecRing_ == input.fecRing() &&
       ccuAddr_ == input.ccuAddr() &&
       ccuChan_ == input.ccuChan() &&
       lldChan_ == input.lldChan() &&
       i2cAddr_ == input.i2cAddr() ) { 
    return true;
  } else { return false; }
}

// -----------------------------------------------------------------------------
// 
bool SiStripFecKey::isConsistent( const SiStripKey& key ) const {
  SiStripKey& temp = const_cast<SiStripKey&>(key);
  SiStripFecKey& input = dynamic_cast<SiStripFecKey&>(temp);
  if ( !(&input) ) { return false; }
  if ( isEqual(input) ) { return true; }
  else if ( ( fecCrate_ == 0 || input.fecCrate() == 0 ) &&
	    ( fecSlot_ == 0 || input.fecSlot() == 0 ) &&
	    ( fecRing_ == 0 || input.fecRing() == 0 ) &&
	    ( ccuAddr_ == 0 || input.ccuAddr() == 0 ) &&
	    ( lldChan_ == 0 || input.lldChan() == 0 ) &&
	    ( i2cAddr_ == 0 || input.i2cAddr() == 0 ) ) {
    return true;
  } else { return false; }
}

// -----------------------------------------------------------------------------
//
bool SiStripFecKey::isValid() const { 
  return isValid(sistrip::APV); 
}

// -----------------------------------------------------------------------------
//
bool SiStripFecKey::isValid( const sistrip::Granularity& gran ) const {
  if ( gran == sistrip::FEC_SYSTEM ) { return true; }
  else if ( gran == sistrip::UNDEFINED_GRAN ) { return false; }
  
  if ( fecCrate_ != sistrip::invalid_ ) {
    if ( gran == sistrip::FEC_CRATE ) { return true; }
    if ( fecSlot_ != sistrip::invalid_ ) {
      if ( gran == sistrip::FEC_RING ) { return true; }
      if ( fecRing_ != sistrip::invalid_ ) {
	if ( gran == sistrip::FEC_RING ) { return true; }
	if ( ccuAddr_ != sistrip::invalid_ ) {
	  if ( gran == sistrip::CCU_ADDR ) { return true; }
	  if ( ccuChan_ != sistrip::invalid_ ) {
	    if ( gran == sistrip::CCU_CHAN ) { return true; }
	    if ( lldChan_ != sistrip::invalid_ ) {
	      if ( gran == sistrip::LLD_CHAN ) { return true; }
	      if ( i2cAddr_ != sistrip::invalid_ ) {
		if ( gran == sistrip::APV ) { return true; }
	      }
	    }
	  }
	}
      }
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
//
bool SiStripFecKey::isInvalid() const { 
  return isValid(sistrip::APV); 
}

// -----------------------------------------------------------------------------
//
bool SiStripFecKey::isInvalid( const sistrip::Granularity& gran ) const {
  if ( gran == sistrip::FEC_SYSTEM ) { return false; }
  else if ( gran == sistrip::UNDEFINED_GRAN ) { return true; }

  if ( fecCrate_ == sistrip::invalid_ ) {
    if ( gran == sistrip::FEC_CRATE ) { return true; }
    if ( fecSlot_ == sistrip::invalid_ ) {
      if ( gran == sistrip::FEC_RING ) { return true; }
      if ( fecRing_ == sistrip::invalid_ ) {
	if ( gran == sistrip::FEC_RING ) { return true; }
	if ( ccuAddr_ == sistrip::invalid_ ) {
	  if ( gran == sistrip::CCU_ADDR ) { return true; }
	  if ( ccuChan_ == sistrip::invalid_ ) {
	    if ( gran == sistrip::CCU_CHAN ) { return true; }
	    if ( lldChan_ == sistrip::invalid_ ) {
	      if ( gran == sistrip::LLD_CHAN  ) { return true; }
	      if ( i2cAddr_ == sistrip::invalid_ ) {
		if ( gran == sistrip::APV  ) { return true; }
	      }
	    }
	  }
	}
      }
    }
  }
  return false;
}

// -----------------------------------------------------------------------------
//
void SiStripFecKey::initFromValue() {

  // FEC crate  
  if ( fecCrate_ >= sistrip::FEC_CRATE_MIN &&
       fecCrate_ <= sistrip::FEC_CRATE_MAX ) {
    fecCrate_ = fecCrate_;
  } else if ( fecCrate_ == 0 ) { 
    fecCrate_ = 0;
  } else { fecCrate_ = sistrip::invalid_; }

  // FEC slot
  if ( fecSlot_ >= sistrip::CRATE_SLOT_MIN &&
       fecSlot_ <= sistrip::CRATE_SLOT_MAX ) {
    fecSlot_ = fecSlot_;
  } else if ( fecSlot_ == 0 ) { 
    fecSlot_ = 0;
  } else { fecSlot_ = sistrip::invalid_; }

  // FEC ring
  if ( fecRing_ >= sistrip::FEC_RING_MIN &&
       fecRing_ <= sistrip::FEC_RING_MAX ) {
    fecRing_ = fecRing_;
  } else if ( fecRing_ == 0 ) { 
    fecRing_ = 0;
  } else { fecRing_ = sistrip::invalid_; }

  // CCU addr
  if ( ccuAddr_ >= sistrip::CCU_ADDR_MIN &&
       ccuAddr_ <= sistrip::CCU_ADDR_MAX ) {
    ccuAddr_ = ccuAddr_;
  } else if ( ccuAddr_ == 0 ) { 
    ccuAddr_ = 0;
  } else { ccuAddr_ = sistrip::invalid_; }

  // CCU chan
  if ( ccuChan_ >= sistrip::CCU_CHAN_MIN &&
       ccuChan_ <= sistrip::CCU_CHAN_MAX ) {
    ccuChan_ = ccuChan_;
  } else if ( ccuChan_ == 0 ) { 
    ccuChan_ = 0;
  } else { ccuChan_ = sistrip::invalid_; }
  
  // LLD channel
  if ( lldChan_ >= sistrip::LLD_CHAN_MIN &&
       lldChan_ <= sistrip::LLD_CHAN_MAX ) {
    lldChan_ = lldChan_;
  } else if ( lldChan_ == 0 ) { 
    lldChan_ = 0;
  } else { lldChan_ = sistrip::invalid_; }
  
  // APV I2C address
  if ( i2cAddr_ >= sistrip::APV_I2C_MIN &&
       i2cAddr_ <= sistrip::APV_I2C_MAX ) { 
    i2cAddr_ = i2cAddr_;
//     // Consistency check wrt LLD channel
//     if ( lldChan(i2cAddr_) && lldChan_ &&
// 	 lldChan(i2cAddr_) != lldChan_ ) {
//       i2cAddr_ = sistrip::invalid_;
//     }
  } else if ( i2cAddr_ == 0 ) { 
    i2cAddr_ = 0;
  } else { i2cAddr_ = sistrip::invalid_; }
  
}

// -----------------------------------------------------------------------------
//
void SiStripFecKey::initFromKey() {
  
  if ( key() == sistrip::invalid32_ ) { 

    // ---------- Set FecKey based on member data ----------
    
    // Initialise to null value
    key(0);
    
    // Extract FEC crate  
    if ( fecCrate_ >= sistrip::FEC_CRATE_MIN &&
	 fecCrate_ <= sistrip::FEC_CRATE_MAX ) {
      key( key() | (fecCrate_<<fecCrateOffset_) );
    } else if ( fecCrate_ == 0 ) { 
      key( key() | (fecCrate_<<fecCrateOffset_) );
    } else { 
      key( key() | (fecCrateMask_<<fecCrateOffset_) ); 
    }

    // Extract FEC slot
    if ( fecSlot_ >= sistrip::CRATE_SLOT_MIN &&
	 fecSlot_ <= sistrip::CRATE_SLOT_MAX ) {
      key( key() | (fecSlot_<<fecSlotOffset_) );
    } else if ( fecSlot_ == 0 ) { 
      key( key() | (fecSlot_<<fecSlotOffset_) );
    } else { 
      key( key() | (fecSlotMask_<<fecSlotOffset_) ); 
    }

    // Extract FEC ring
    if ( fecRing_ >= sistrip::FEC_RING_MIN &&
	 fecRing_ <= sistrip::FEC_RING_MAX ) {
      key( key() | (fecRing_<<fecRingOffset_) );
    } else if ( fecRing_ == 0 ) { 
      key( key() | (fecRing_<<fecRingOffset_) );
    } else { 
      key( key() | (fecRingMask_<<fecRingOffset_) ); 
    }

    // Extract CCU addr
    if ( ccuAddr_ >= sistrip::CCU_ADDR_MIN &&
	 ccuAddr_ <= sistrip::CCU_ADDR_MAX ) {
      key( key() | (ccuAddr_<<ccuAddrOffset_) );
    } else if ( ccuAddr_ == 0 ) { 
      key( key() | (ccuAddr_<<ccuAddrOffset_) );
    } else { 
      key( key() | (ccuAddrMask_<<ccuAddrOffset_) ); 
    }

    // Extract CCU chan
    if ( ccuChan_ >= sistrip::CCU_CHAN_MIN &&
	 ccuChan_ <= sistrip::CCU_CHAN_MAX ) {
      key( key() | ( (ccuChan_-(sistrip::CCU_CHAN_MIN-1)) << ccuChanOffset_ ) ); 
    } else if ( ccuChan_ == 0 ) { 
      key( key() | (ccuChan_<<ccuChanOffset_) );
    } else { 
      key( key() | (ccuChanMask_<<ccuChanOffset_) ); 
    }
    
    // Extract LLD channel
    if ( lldChan_ >= sistrip::LLD_CHAN_MIN &&
	 lldChan_ <= sistrip::LLD_CHAN_MAX ) {
      key( key() | (lldChan_<<lldChanOffset_) ); 
    } else if ( lldChan_ == 0 ) { 
      key( key() | (lldChan_<<lldChanOffset_) );
    } else { 
      key( key() | (lldChanMask_<<lldChanOffset_) ); 
    }
    
    // Extract APV I2C address
    if ( i2cAddr_ >= sistrip::APV_I2C_MIN &&
	 i2cAddr_ <= sistrip::APV_I2C_MAX ) {
      key( key() | ( hybridPos(i2cAddr_) << i2cAddrOffset_ ) ); 
//       // Consistency check wrt LLD channel
//       if ( lldChan(i2cAddr_) && lldChan_ &&
// 	   lldChan(i2cAddr_) != lldChan_ ) {
// 	i2cAddr_ |= (i2cAddrMask_<<i2cAddrOffset_) ); 
//       }
    } else if ( i2cAddr_ == 0 ) { 
      key( key() | (i2cAddr_<<i2cAddrOffset_) );
    } else { 
      key( key() | (i2cAddrMask_<<i2cAddrOffset_) ); 
    }
    
  } else {
    
    // ---------- Set member data based on FEC key ----------

    fecCrate_ = ( key()>>fecCrateOffset_ ) & fecCrateMask_;
    fecSlot_  = ( key()>>fecSlotOffset_ )  & fecSlotMask_;
    fecRing_  = ( key()>>fecRingOffset_ )  & fecRingMask_;
    ccuAddr_  = ( key()>>ccuAddrOffset_ )  & ccuAddrMask_;
    ccuChan_  = ( key()>>ccuChanOffset_ )  & ccuChanMask_;
    lldChan_  = ( key()>>lldChanOffset_ )  & lldChanMask_;
    i2cAddr_  = ( key()>>i2cAddrOffset_ )  & i2cAddrMask_;

    if ( fecCrate_ == fecCrateMask_ ) { fecCrate_ = sistrip::invalid_; } 
    if ( fecSlot_ == fecSlotMask_ )   { fecSlot_ = sistrip::invalid_; } 
    if ( fecRing_ == fecRingMask_ )   { fecRing_ = sistrip::invalid_; } 
    if ( ccuAddr_ == ccuAddrMask_ )   { ccuAddr_ = sistrip::invalid_; } 
    if ( ccuChan_ == ccuChanMask_ )   { ccuChan_ = sistrip::invalid_; }
    else if ( ccuChan_ )              { ccuChan_ += (sistrip::CCU_CHAN_MIN-1); }
    if ( lldChan_ == lldChanMask_ )   { lldChan_ = sistrip::invalid_; }
    if ( i2cAddr_ == i2cAddrMask_ )   { i2cAddr_ = sistrip::invalid_; }
    else if ( i2cAddr_ )              { i2cAddr_ = i2cAddr(i2cAddr_); }
    
//     // Consistency check wrt LLD channel
//     if ( lldChan(i2cAddr_) && lldChan_ &&
// 	 lldChan(i2cAddr_) != lldChan_ ) {
//       i2cAddr_ = sistrip::invalid_;
//     }
    
  }
  
}

// -----------------------------------------------------------------------------
// 
void SiStripFecKey::initFromPath() {
  
  if ( path() == sistrip::null_ ) {
    
    // ---------- Set directory path based on member data ----------

    std::stringstream dir;
    
    dir << sistrip::root_ << sistrip::dir_ 
	<< sistrip::controlView_ << sistrip::dir_;

    // Add FEC crate
    if ( fecCrate_ ) {
      dir << sistrip::fecCrate_ << fecCrate_ << sistrip::dir_;
      
      // Add FEC slot
      if ( fecSlot_ ) {
	dir << sistrip::fecSlot_ << fecSlot_ << sistrip::dir_;
	
	// Add FEC ring
	if ( fecRing_ ) {
	  dir << sistrip::fecRing_ << fecRing_ << sistrip::dir_;
	  
	  // Add CCU address
	  if ( ccuAddr_ ) {
	    dir << sistrip::ccuAddr_ << ccuAddr_ << sistrip::dir_;
	    
	    // Add CCU channel
	    if ( ccuChan_ ) {
	      dir << sistrip::ccuChan_ << ccuChan_ << sistrip::dir_;

	      // Add LLD channel
	      if ( lldChan_ ) {
		dir << sistrip::lldChan_ << lldChan_ << sistrip::dir_;

		// Add APV I2C address
		if ( i2cAddr_ ) {
		  dir << sistrip::apv_ << i2cAddr_ << sistrip::dir_;
		}
	      }
	    }
	  }
	}
      }
    }
    
    std::string temp( dir.str() );
    path( temp );

  } else {
    
    // ---------- Set member data based on directory path ----------
    
    fecCrate_ = 0;
    fecSlot_  = 0;
    fecRing_  = 0;
    ccuAddr_  = 0;
    ccuChan_  = 0;
    lldChan_  = 0;
    i2cAddr_  = 0;
    
    uint32_t curr = 0; // current string position
    uint32_t next = 0; // next string position
    next = path().find( sistrip::controlView_, curr );

    // Extract view 
    curr = next;
    if ( curr != std::string::npos ) { 
      next = path().find( sistrip::fecCrate_, curr );
      std::string control_view( path(), 
				curr+sistrip::controlView_.size(), 
				(next-sistrip::dir_.size())-curr );
      
      // Extract FEC crate
      curr = next;
      if ( curr != std::string::npos ) { 
	next = path().find( sistrip::fecSlot_, curr );
	std::string fec_crate( path(), 
			       curr+sistrip::fecCrate_.size(), 
			       (next-sistrip::dir_.size())-curr );
	fecCrate_ = std::atoi( fec_crate.c_str() );

	// Extract FEC slot
	curr = next;
	if ( curr != std::string::npos ) { 
	  next = path().find( sistrip::fecRing_, curr );
	  std::string fec_slot( path(), 
				curr+sistrip::fecSlot_.size(), 
				(next-sistrip::dir_.size())-curr );
	  fecSlot_ = std::atoi( fec_slot.c_str() );

	  // Extract FEC ring
	  curr = next;
	  if ( curr != std::string::npos ) { 
	    next = path().find( sistrip::ccuAddr_, curr );
	    std::string fec_ring( path(), 
				  curr+sistrip::fecRing_.size(),
				  (next-sistrip::dir_.size())-curr );
	    fecRing_ = std::atoi( fec_ring.c_str() );

	    // Extract CCU address
	    curr = next;
	    if ( curr != std::string::npos ) { 
	      next = path().find( sistrip::ccuChan_, curr );
	      std::string ccu_addr( path(), 
				    curr+sistrip::ccuAddr_.size(), 
				    (next-sistrip::dir_.size())-curr );
	      ccuAddr_ = std::atoi( ccu_addr.c_str() );

	      // Extract CCU channel
	      curr = next;
	      if ( curr != std::string::npos ) { 
		next = path().find( sistrip::lldChan_, curr );
		std::string ccu_chan( path(), 
				      curr+sistrip::ccuChan_.size(), 
				      (next-sistrip::dir_.size())-curr );
		ccuChan_ = std::atoi( ccu_chan.c_str() );
		
		// Extract LLD channel
		curr = next;
		if ( curr != std::string::npos ) { 
		  next = path().find( sistrip::apv_, curr );
		  std::string lld_chan( path(), 
					curr+sistrip::lldChan_.size(), 
					(next-sistrip::dir_.size())-curr );
		  lldChan_ = std::atoi( lld_chan.c_str() );
		  
		  // Extract I2C address
		  curr = next;
		  if ( curr != std::string::npos ) { 
		    next = std::string::npos;
		    std::string i2c_addr( path(), 
					  curr+sistrip::apv_.size(),
					  next-curr );
		    i2cAddr_ = std::atoi( i2c_addr.c_str() );
		  }
		}
	      }
	    }
	  }
	}
      }
    } else {
      std::stringstream ss;
      ss << sistrip::root_ << sistrip::dir_ 
	 << sistrip::unknownView_ << sistrip::dir_;
      std::string temp( ss.str() );
      path( temp );
    }
    
  }
  
}

// -----------------------------------------------------------------------------
// 
void SiStripFecKey::initGranularity() {
  
  granularity( sistrip::FEC_SYSTEM );
  channel(0);
  if ( fecCrate_ && fecCrate_ != sistrip::invalid_ ) {
    granularity( sistrip::FEC_CRATE ); 
    channel(fecCrate_);
    if ( fecSlot_ && fecSlot_ != sistrip::invalid_ ) {
      granularity( sistrip::FEC_SLOT );
      channel(fecSlot_);
      if ( fecRing_ && fecRing_ != sistrip::invalid_ ) {
	granularity( sistrip::FEC_RING );
	channel(fecRing_);
	if ( ccuAddr_ && ccuAddr_ != sistrip::invalid_ ) {
	  granularity( sistrip::CCU_ADDR );
	  channel(ccuAddr_);
	  if ( ccuChan_ && ccuChan_ != sistrip::invalid_ ) {
	    granularity( sistrip::CCU_CHAN );
	    channel(ccuChan_);
	    if ( lldChan_ && lldChan_ != sistrip::invalid_ ) {
	      granularity( sistrip::LLD_CHAN );
	      channel(lldChan_);
	      if ( i2cAddr_ && i2cAddr_ != sistrip::invalid_ ) {
		granularity( sistrip::APV );
		channel(i2cAddr_);
	      } else if ( i2cAddr_ == sistrip::invalid_ ) { 
		granularity( sistrip::UNDEFINED_GRAN ); 
		channel(sistrip::invalid_);
	      }
	    } else if ( lldChan_ == sistrip::invalid_ ) { 
	      granularity( sistrip::UNDEFINED_GRAN ); 
	      channel(sistrip::invalid_);
	    }
	  } else if ( ccuChan_ == sistrip::invalid_ ) { 
	    granularity( sistrip::UNDEFINED_GRAN ); 
	    channel(sistrip::invalid_);
	  }
	} else if ( ccuAddr_ == sistrip::invalid_ ) { 
	  granularity( sistrip::UNDEFINED_GRAN ); 
	  channel(sistrip::invalid_);
	}
      } else if ( fecRing_ == sistrip::invalid_ ) { 
	granularity( sistrip::UNDEFINED_GRAN ); 
	channel(sistrip::invalid_);
      }
    } else if ( fecSlot_ == sistrip::invalid_ ) { 
      granularity( sistrip::UNDEFINED_GRAN ); 
      channel(sistrip::invalid_);
    }
  } else if ( fecCrate_ == sistrip::invalid_ ) { 
    granularity( sistrip::UNDEFINED_GRAN ); 
    channel(sistrip::invalid_);
  }

}

// -----------------------------------------------------------------------------
//
std::ostream& operator<< ( std::ostream& os, const SiStripFecKey& input ) {
  return os << std::endl
	    << " [SiStripFecKey::print]" << std::endl
	    << std::hex
	    << " FEC key              : 0x" 
	    << std::setfill('0') 
	    << std::setw(8) << input.key() << std::endl
	    << std::setfill(' ') 
	    << std::dec
	    << " FEC VME crate        : " << input.fecCrate() << std::endl
	    << " FEC VME slot         : " << input.fecSlot() << std::endl 
	    << " FEC control ring     : " << input.fecRing() << std::endl
	    << " CCU I2C address      : " << input.ccuAddr() << std::endl
	    << " CCU chan (FE module) : " << input.ccuChan() << std::endl
	    << " LaserDriver channel  : " << input.lldChan() << std::endl 
	    << " APV I2C address      : " << input.i2cAddr() << std::endl 
	    << " Directory            : " << input.path() << std::endl
	    << " Granularity          : "
	    << SiStripEnumsAndStrings::granularity( input.granularity() ) << std::endl
 	    << " Channel              : " << input.channel() << std::endl
	    << " isValid              : " << input.isValid();
}

