// Last commit: $Id: SiStripHistoTitle.cc,v 1.2 2007/03/21 08:22:59 bainbrid Exp $

#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include "DataFormats/SiStripCommon/interface/SiStripKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripFedKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripFecKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripNullKey.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include <iostream>
#include <iomanip>

// -----------------------------------------------------------------------------
//
SiStripHistoTitle::SiStripHistoTitle( const sistrip::RunType& run_type, 
				      const SiStripKey&       key_object,
				      const std::string&      extra_info ) 
  : title_(""),
    runType_(run_type),
    keyType_(sistrip::UNKNOWN_KEY),
    keyValue_(sistrip::invalid32_),
    granularity_(sistrip::UNKNOWN_GRAN),
    channel_(sistrip::invalid_),
    extraInfo_(extra_info)
{
  SiStripKey& temp = const_cast<SiStripKey&>(key_object);
  if ( &dynamic_cast<SiStripFedKey&>(temp) ) {
    keyType_ = sistrip::FED_KEY;
  } else if ( &dynamic_cast<SiStripFecKey&>(temp) ) {
    keyType_ = sistrip::FEC_KEY;
  } else {
    keyType_ = sistrip::UNKNOWN_KEY;
  }
  if ( &key_object ) {
    keyValue_ = key_object.key();
    granularity_ = key_object.granularity();
    channel_ = key_object.channel();
  }
  setTitle();
}

// -----------------------------------------------------------------------------
//
SiStripHistoTitle::SiStripHistoTitle( const sistrip::RunType&     run_type, 
				      const sistrip::KeyType&     key_type,
				      const uint32_t&             key_value,
				      const sistrip::Granularity& gran,
				      const uint16_t&             channel,
				      const std::string&          extra_info ) 
  : title_(""),
    runType_(run_type),
    keyType_(key_type),
    keyValue_(key_value),
    granularity_(gran),
    channel_(channel),
    extraInfo_(extra_info)
{
  setTitle();
}

// -----------------------------------------------------------------------------
//
SiStripHistoTitle::SiStripHistoTitle( const std::string& histo_title ) 
  : title_(histo_title),
    runType_(sistrip::UNDEFINED_RUN_TYPE),
    keyType_(sistrip::UNDEFINED_KEY),
    keyValue_(sistrip::invalid32_),
    granularity_(sistrip::UNDEFINED_GRAN),
    channel_(sistrip::invalid_),
    extraInfo_("")
{
  extractTitle();
}

// -----------------------------------------------------------------------------
//
void SiStripHistoTitle::setTitle() {

  std::stringstream title;

  // Append RunType, KeyType and KeyValue
  title << SiStripEnumsAndStrings::runType( runType_ )
	<< sistrip::sep_
	<< SiStripEnumsAndStrings::keyType( keyType_ )
	<< sistrip::hex_ 
	<< std::setfill('0') << std::setw(8) << std::hex << keyValue_ << std::dec
	<< sistrip::sep_;
  
  // Append Granularity and channel number
  title << SiStripEnumsAndStrings::granularity( granularity_ );
  if ( channel_ ) { title << channel_; }
  
  // Append extra info
  if ( extraInfo_ != "" ) { 
    title << sistrip::sep_ << extraInfo_; 
  }
  
  title_ = title.str();

}

// -----------------------------------------------------------------------------
// 
void SiStripHistoTitle::extractTitle() {
  
  std::string::size_type length = title_.length();
  std::string::size_type position = 0;
  std::string::size_type pos = 0;
  std::string::size_type siz = 0;
  
  // Extract RunType
  siz = title_.find(sistrip::sep_,position) - position;
  runType_ = SiStripEnumsAndStrings::runType( title_.substr(position,siz) );
  std::string run_type = SiStripEnumsAndStrings::runType( runType_ );
  position += title_.substr(position).find( run_type ) + run_type.size() + sistrip::sep_.size();
  if ( position >= length ) { return; }
  
  // Extract KeyType
  siz = title_.find(sistrip::sep_,position) - position;
  keyType_ = SiStripEnumsAndStrings::keyType( title_.substr(position,siz) );
  std::string key_type = SiStripEnumsAndStrings::keyType( keyType_ );
  position += title_.substr(position).find( key_type ) + key_type.size() + sistrip::hex_.size();
  if ( position >= length ) { return; }
  
  // Extract KeyValue
  siz = 8;
  std::stringstream key; 
  key << title_.substr(position,siz);
  key >> std::hex >> keyValue_;
  position += siz + sistrip::sep_.size();
  if ( position >= length ) { return; }
  
  // Extract Granularity
  pos = title_.find(sistrip::sep_,position);
  if ( pos == std::string::npos || pos < position ) { siz = std::string::npos - position; }
  else { siz = pos - position; }
  granularity_ = SiStripEnumsAndStrings::granularity( title_.substr(position,siz) );
  std::string gran = SiStripEnumsAndStrings::granularity( granularity_ );
  position += title_.substr(position).find( gran ) + gran.size();
  if ( position > length ) { return; }

  // Extract Channel 
  pos = title_.find(sistrip::sep_,position);
  if ( pos == std::string::npos || pos < position ) { siz = std::string::npos - position; }
  else { siz = pos - position; }
  if ( position == length || !siz ) {
    if ( granularity_ != sistrip::UNDEFINED_GRAN ) { channel_ = 0; }
    else if ( granularity_ == sistrip::UNKNOWN_GRAN ) { channel_ = sistrip::invalid_; }
  } else {
    std::stringstream chan; 
    chan << title_.substr(position,siz);
    chan >> std::dec >> channel_;
  }
  position += siz + sistrip::sep_.size();
  if ( position >= length ) { return; }
  
  // Extract ExtraInfo
  extraInfo_ = title_.substr( position, std::string::npos - position ); 
  
}

// -----------------------------------------------------------------------------
//
std::ostream& operator<< ( std::ostream& os, const SiStripHistoTitle& title ) {
  std::stringstream ss;
  ss << "[SiStripHistoTitle::print]" << std::endl
     << " Title          : " << title.title() << std::endl
     << " RunType        : " << SiStripEnumsAndStrings::runType( title.runType() ) << std::endl
     << " KeyType        : " << SiStripEnumsAndStrings::keyType( title.keyType() ) << std::endl
     << " KeyValue (hex) : " << std::hex << std::setfill('0') << std::setw(8) << title.keyValue() << std::dec << std::endl
     << " Granularity    : " << SiStripEnumsAndStrings::granularity( title.granularity() ) << std::endl
     << " Channel        : " << title.channel() << std::endl
     << " ExtraInfo      : ";
  if ( title.extraInfo() != "" ) { ss << "\"" << title.extraInfo() << "\""; }
  else { ss << "(none)"; }
  os << ss.str();
  return os;
}

