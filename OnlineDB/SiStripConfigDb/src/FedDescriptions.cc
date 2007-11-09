// Last commit: $Id: FedDescriptions.cc,v 1.12 2007/11/07 15:55:33 bainbrid Exp $
// Latest tag:  $Name:  $
// Location:    $Source: /cvs_server/repositories/CMSSW/CMSSW/OnlineDB/SiStripConfigDb/src/FedDescriptions.cc,v $

#include "OnlineDB/SiStripConfigDb/interface/SiStripConfigDb.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "CondFormats/SiStripObjects/interface/FedChannelConnection.h"

using namespace std;
using namespace sistrip;

// -----------------------------------------------------------------------------
// 
const SiStripConfigDb::FedDescriptions& SiStripConfigDb::getFedDescriptions() {

  if ( !deviceFactory(__func__) ) { return feds_; }
  if ( !resetFeds_ ) { return feds_; }

  try {
    deviceFactory(__func__)->setUsingStrips( usingStrips_ );
    
    //@@ handle fact that "latest" version for FED is signified by "-1", not "0"!
    int16_t major = dbParams_.fedMajor_; 
    int16_t minor = dbParams_.fedMinor_; 
    if ( dbParams_.fedMajor_ == 0 && 
	 dbParams_.fedMinor_ == 0 ) {
      major = -1;
      minor = -1;
    }
    feds_ = *( deviceFactory(__func__)->getFed9UDescriptions( dbParams_.partition_, 
							      major, 
							      minor ) );
    resetFeds_ = false;
  }
  catch (... ) {
    handleException( __func__ );
  }
  
  // Debug 
  stringstream ss; 
  ss << "[SiStripConfigDb::" << __func__ << "]";
  if ( feds_.empty() ) { ss << " Found no FED descriptions"; }
  else { ss << " Found " << feds_.size() << " FED descriptions"; }
  if ( !dbParams_.usingDb_ ) { ss << " in " << dbParams_.inputFedXml_.size() << " 'fed.xml' file(s)"; }
  else { ss << " in database partition '" << dbParams_.partition_ << "'"; }
  if ( feds_.empty() ) { edm::LogWarning(mlConfigDb_) << ss; }
  else { LogTrace(mlConfigDb_) << ss; }
  
  return feds_;
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::resetFedDescriptions() {
  feds_.clear();
  resetFeds_ = true;
}

// -----------------------------------------------------------------------------
// 
void SiStripConfigDb::uploadFedDescriptions( bool new_major_version ) { //@@ this ok???

  if ( !deviceFactory(__func__) ) { return; }
  
  try { 
    SiStripConfigDb::FedDescriptions::iterator ifed = feds_.begin();
    for ( ; ifed != feds_.end(); ifed++ ) {
      deviceFactory(__func__)->setFed9UDescription( **ifed, 
 						    (uint16_t*)(&dbParams_.fedMajor_), 
 						    (uint16_t*)(&dbParams_.fedMinor_),
 						    false ); // (new_major_version?1:0) ); //@@ should be true???
    }
  }
  catch (...) { 
    handleException( __func__ ); 
  }
  
}

// -----------------------------------------------------------------------------
// 
const SiStripConfigDb::FedDescriptions& SiStripConfigDb::createFedDescriptions( const SiStripFecCabling& fec_cabling ) {
  
  // Static container
  static FedDescriptions static_fed_descriptions;
  static_fed_descriptions.clear();
  
  // Create FED cabling from FEC cabling
  vector<FedChannelConnection> conns;
  fec_cabling.connections( conns );
  SiStripFedCabling* fed_cabling = new SiStripFedCabling( conns );

  // Retrieve and iterate through FED ids
  vector<uint16_t>::const_iterator ifed = fed_cabling->feds().begin();
  for ( ; ifed != fed_cabling->feds().end(); ifed++ ) {
    // Create FED description
    try {
      Fed9U::Fed9UDescription* f = new Fed9U::Fed9UDescription();
      f->setFedId( *ifed );
      f->setFedHardwareId( *ifed );
      Fed9U::Fed9UAddress addr;
      for ( uint32_t i = 0; i < Fed9U::APVS_PER_FED; i++ ) {
	addr.setFedApv(i);
	vector<Fed9U::u32> pedestals(128,100);
	vector<bool> disableStrips(128,false);
	vector<Fed9U::u32> highThresholds(128,50);
	vector<Fed9U::u32> lowThresholds(128,20);
	vector<Fed9U::Fed9UStripDescription> apvStripDescription(128);
	for ( uint32_t j = 0; j < Fed9U::STRIPS_PER_APV; j++) {
	  apvStripDescription[j].setPedestal(pedestals[j]);
	  apvStripDescription[j].setDisable(disableStrips[j]);
	  apvStripDescription[j].setLowThreshold(lowThresholds[j]);
	  apvStripDescription[j].setHighThreshold(highThresholds[j]);
	}
	f->getFedStrips().setApvStrips (addr, apvStripDescription);
      }
      static_fed_descriptions.push_back( f );
    } catch(...) {
      stringstream ss; 
      ss << "Problems creating description for FED id " << *ifed;
      handleException( __func__, ss.str() );
    }
  } 
  
  if ( static_fed_descriptions.empty() ) {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " No FED connections created!";
  }
  
  return static_fed_descriptions;

}

// -----------------------------------------------------------------------------
/** */ 
const vector<uint16_t>& SiStripConfigDb::getFedIds() {
  
  static vector<uint16_t> fed_ids;
  fed_ids.clear();
  
  bool using_strips = usingStrips_;
  deviceFactory(__func__)->setUsingStrips( false );
  getFedDescriptions();
  deviceFactory(__func__)->setUsingStrips( using_strips );
  FedDescriptions::iterator ifed = feds_.begin();
  for ( ; ifed != feds_.end(); ifed++ ) { 
    fed_ids.push_back( (*ifed)->getFedId() );
  }
  if ( fed_ids.empty() ) {
    edm::LogWarning(mlConfigDb_)
      << "[SiStripConfigDb::" << __func__ << "]"
      << " No FED ids found!"; 
  }
  return fed_ids;
}
