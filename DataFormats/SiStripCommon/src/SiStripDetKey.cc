// Last commit: $Id: SiStripDetKey.cc,v 1.7 2007/07/31 15:20:25 ratnik Exp $

#include "DataFormats/SiStripCommon/interface/SiStripDetKey.h"
#include "DataFormats/SiStripCommon/interface/Constants.h" 
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include <iomanip>

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const DetId& det_id,
			      const uint16_t& apv_pair_number,
			      const uint16_t& apv_within_pair ) :
  SiStripKey(),
  apvPairNumber_(apv_pair_number), 
  apvWithinPair_(apv_within_pair)
{
  // order is important!
  initFromValue();
  initFromKey();
  initFromPath();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const SiStripDetId& det_id ) :
  SiStripKey(),
  apvPairNumber_(sistrip::invalid_), 
  apvWithinPair_(sistrip::invalid_)
{
  // order is important!
  initFromValue();
  initFromKey();
  initFromPath();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const uint32_t& det_key ) :
  SiStripKey(det_key),
  apvPairNumber_(sistrip::invalid_), 
  apvWithinPair_(sistrip::invalid_)
{
  // order is important!
  initFromKey(); 
  initFromValue();
  initFromPath();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const std::string& path ) :
  SiStripKey(path),
  apvPairNumber_(sistrip::invalid_), 
  apvWithinPair_(sistrip::invalid_)
{
  // order is important!
  initFromPath();
  initFromValue();
  initFromKey();
  initGranularity();
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const SiStripDetKey& input ) :
  SiStripKey(),
  apvPairNumber_(input.apvPairNumber()), 
  apvWithinPair_(input.apvWithinPair())
{
  key(input.key());
  path(input.path());
  granularity(input.granularity());
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey( const SiStripKey& input ) :
  SiStripKey(),
  apvPairNumber_(sistrip::invalid_), 
  apvWithinPair_(sistrip::invalid_)
{
  SiStripKey& temp = const_cast<SiStripKey&>(input);
  SiStripDetKey& det_key = dynamic_cast<SiStripDetKey&>(temp);
  if ( (&det_key) ) {
    key(det_key.key());
    path(det_key.path());
    granularity(det_key.granularity());
    apvPairNumber_ = det_key.apvPairNumber();
    apvWithinPair_ = det_key.apvWithinPair();
  }
}

// -----------------------------------------------------------------------------
// 
SiStripDetKey::SiStripDetKey() : 
  SiStripKey(),
  apvPairNumber_(sistrip::invalid_), 
  apvWithinPair_(sistrip::invalid_)
{;}

// -----------------------------------------------------------------------------
// 
bool SiStripDetKey::isEqual( const SiStripKey& input ) const {
  SiStripKey& temp = const_cast<SiStripKey&>(input);
  if ( &dynamic_cast<SiStripDetKey&>(temp) ) { return true; }
  else { return false; }
}

// -----------------------------------------------------------------------------
// 
bool SiStripDetKey::isConsistent( const SiStripKey& input ) const {
  return isEqual(input);
}

// -----------------------------------------------------------------------------
//
bool SiStripDetKey::isValid() const { 
  return false;
}

// -----------------------------------------------------------------------------
//
bool SiStripDetKey::isValid( const sistrip::Granularity& gran ) const {
  return false; 
}

// -----------------------------------------------------------------------------
//
bool SiStripDetKey::isInvalid() const { 
  return true;
}

// -----------------------------------------------------------------------------
//
bool SiStripDetKey::isInvalid( const sistrip::Granularity& gran ) const {
  return true;
}

// -----------------------------------------------------------------------------
// 
void SiStripDetKey::initFromValue() {;}

// -----------------------------------------------------------------------------
//
void SiStripDetKey::initFromKey() {;}

// -----------------------------------------------------------------------------
// 
void SiStripDetKey::initFromPath() {;}

// -----------------------------------------------------------------------------
// 
void SiStripDetKey::initGranularity() {;}

// -----------------------------------------------------------------------------
//
void SiStripDetKey::print( std::stringstream& ss ) const {
  ss << " [SiStripDetKey::print]" << std::endl
     << std::hex
     << " 32-bit key  : 0x" 
     << std::setfill('0') 
     << std::setw(8) << key() << std::endl
     << std::setfill(' ') 
     << std::dec
     << " Directory   : " << path() << std::endl
     << " Granularity : "
     << SiStripEnumsAndStrings::granularity( granularity() ) << std::endl
     << " Channel     : " << channel() << std::endl
     << " isValid    : " << isValid();
}

// -----------------------------------------------------------------------------
//
std::ostream& operator<< ( std::ostream& os, const SiStripDetKey& input ) {
  std::stringstream ss;
  input.print(ss);
  os << ss.str();
  return os;
}
