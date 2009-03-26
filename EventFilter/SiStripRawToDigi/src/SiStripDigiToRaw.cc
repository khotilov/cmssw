// Last commit: $Id: OldSiStripDigiToRaw.cc,v 1.32 2009/03/04 13:37:50 bainbrid Exp $

#include "EventFilter/SiStripRawToDigi/interface/SiStripDigiToRaw.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DataFormats/SiStripCommon/interface/SiStripConstants.h"
#include "FWCore/Utilities/interface/CRC16.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Fed9UUtils.hh"
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

// -----------------------------------------------------------------------------
/** */
OldSiStripDigiToRaw::OldSiStripDigiToRaw( std::string mode, 
				    int16_t nbytes,
				    bool use_fed_key ) : 
  readoutMode_(mode),
  nAppendedBytes_(nbytes),
  useFedKey_(use_fed_key)
{
  if ( edm::isDebugEnabled() ) {
    LogDebug("DigiToRaw")
      << "[OldSiStripDigiToRaw::OldSiStripDigiToRaw]"
      << " Constructing object...";
  }
}

// -----------------------------------------------------------------------------
/** */
OldSiStripDigiToRaw::~OldSiStripDigiToRaw() {
  if ( edm::isDebugEnabled() ) {
    LogDebug("DigiToRaw")
      << "[OldSiStripDigiToRaw::~OldSiStripDigiToRaw]"
      << " Destructing object...";
  }
}

// -----------------------------------------------------------------------------
/** 
    Input: DetSetVector of SiStripDigis. Output: FEDRawDataCollection.
//     Retrieves and iterates through FED buffers, extract FEDRawData
//     from collection and (optionally) dumps raw data to stdout, locates
//     start of FED buffer by identifying DAQ header, creates new
//     Fed9UEvent object using current FEDRawData buffer, dumps FED
//     buffer to stdout, retrieves data from various header fields
*/
void OldSiStripDigiToRaw::createFedBuffers( edm::Event& event,
					 edm::ESHandle<SiStripFedCabling>& cabling,
					 edm::Handle< edm::DetSetVector<SiStripDigi> >& collection,
					 auto_ptr<FEDRawDataCollection>& buffers ) {

  try {

    const uint16_t strips_per_fed = 96 * 256; 
    vector<uint16_t> raw_data; 
    raw_data.reserve(strips_per_fed);
    
    const vector<uint16_t>& fed_ids = cabling->feds();
    vector<uint16_t>::const_iterator ifed;

    for ( ifed = fed_ids.begin(); ifed != fed_ids.end(); ifed++ ) {

      raw_data.clear(); raw_data.resize( strips_per_fed, 0 );
      
      const std::vector<FedChannelConnection>& conns = cabling->connections(*ifed);
      std::vector<FedChannelConnection>::const_iterator iconn = conns.begin();
      for ( ; iconn != conns.end(); iconn++ ) {

	// Determine FED key from cabling
	uint32_t fed_key = ( ( iconn->fedId() & sistrip::invalid_ ) << 16 ) | ( iconn->fedCh() & sistrip::invalid_ );
	
	// Determine whether DetId or FED key should be used to index digi containers
	uint32_t key = ( useFedKey_ || readoutMode_ == "SCOPE_MODE" ) ? fed_key : iconn->detId();

	// Determine APV pair number (needed only when using DetId)
	uint16_t ipair = ( useFedKey_ || readoutMode_ == "SCOPE_MODE" ) ? 0 : iconn->apvPairNumber();
	
	// Check key is non-zero and valid
	if ( !key || ( key == sistrip::invalid32_ ) ) { continue; }

	// Find digis for DetID in collection
	vector< edm::DetSet<SiStripDigi> >::const_iterator digis = collection->find( key );
	if (digis == collection->end()) { continue; } 

	edm::DetSet<SiStripDigi>::const_iterator idigi;
	for ( idigi = digis->data.begin(); idigi != digis->data.end(); idigi++ ) {
	  if ( (*idigi).strip() < ipair*256 ||
	       (*idigi).strip() > ipair*256+255 ) { continue; }
	  unsigned short strip = iconn->fedCh()*256 + (*idigi).strip() % 256;

	  if ( strip >= strips_per_fed ) {
	    if ( edm::isDebugEnabled() ) {
	      std::stringstream ss;
	      ss << "[OldSiStripDigiToRaw::createFedBuffers]"
		 << " strip >= strips_per_fed";
	      edm::LogWarning("DigiToRaw") << ss.str();
	    }
	    continue;
	  }
	  
	  // check if value already exists
	  if ( edm::isDebugEnabled() ) {
	    if ( raw_data[strip] && raw_data[strip] != (*idigi).adc() ) {
	      std::stringstream ss; 
	      ss << "[OldSiStripDigiToRaw::createFedBuffers]" 
		 << " Incompatible ADC values in buffer!"
		 << "  FedId/FedCh: " << *ifed << "/" << iconn->fedCh()
		 << "  DetStrip: " << (*idigi).strip() 
		 << "  FedStrip: " << strip
		 << "  AdcValue: " << (*idigi).adc()
	       << "  RawData[" << strip << "]: " << raw_data[strip];
	      edm::LogWarning("DigiToRaw") << ss.str();
	    }
	  }

	  // Add digi to buffer
	  raw_data[strip] = (*idigi).adc();

	}
	// if ((*idigi).strip() >= (ipair+1)*256) break;
      }

      // instantiate appropriate buffer creator object depending on readout mode
      Fed9U::Fed9UBufferCreator* creator = 0;
      if ( readoutMode_ == "SCOPE_MODE" ) {
	if ( edm::isDebugEnabled() ) {
	  edm::LogWarning("DigiToRaw")
 	    << "[OldSiStripDigiToRaw::createFedBuffers]" 
	    << " Fed9UBufferCreatorScopeMode not implemented yet!";
	}
      } else if ( readoutMode_ == "VIRGIN_RAW" ) {
	creator = new Fed9U::Fed9UBufferCreatorRaw();
      } else if ( readoutMode_ == "PROCESSED_RAW" ) {
	creator = new Fed9U::Fed9UBufferCreatorProcRaw();
      } else if ( readoutMode_ == "ZERO_SUPPRESSED" ) {
	creator = new Fed9U::Fed9UBufferCreatorZS();
      } else {
	if ( edm::isDebugEnabled() ) {
	  edm::LogWarning("DigiToRaw")
 	    << "[OldSiStripDigiToRaw::createFedBuffers]" 
	    << " UNKNOWN readout mode";
	}
      }
  
      if ( !creator ) { 
	if ( edm::isDebugEnabled() ) {
	  edm::LogWarning("DigiToRaw")
 	    << "[OldSiStripDigiToRaw::createFedBuffers]" 
	    << " NULL pointer to Fed9UBufferCreator";
	}
	return; 
      }

      // generate FED buffer and pass to Daq
      Fed9U::Fed9UBufferGenerator generator( creator );
      generator.generateFed9UBuffer( raw_data );
      //generator.setSlink64();
      FEDRawData& fedrawdata = buffers->FEDData( *ifed ); 
      // calculate size of FED buffer in units of bytes (unsigned char)
      int nbytes = generator.getBufferSize() * 4;
      // resize (public) "data_" member of struct FEDRawData
      fedrawdata.resize( nbytes );
      // copy FED buffer to struct FEDRawData using Fed9UBufferGenerator
      unsigned char* chars = const_cast<unsigned char*>( fedrawdata.data() );
      unsigned int* ints = reinterpret_cast<unsigned int*>( chars );
      generator.getBuffer( ints );
      if ( creator ) { delete creator; }

      //@@ THIS IS TEMPORARY FIX SO THAT DAQ WORKS. CONCERNS 32-BIT WORD SWAPPING
      FEDRawData temp( fedrawdata );
      uint32_t* temp_u32 = reinterpret_cast<uint32_t*>( const_cast<unsigned char*>( temp.data() ) );
      uint32_t* fedrawdata_u32 = reinterpret_cast<uint32_t*>( const_cast<unsigned char*>( fedrawdata.data() ) );
      uint16_t iter = 0; 
      while ( iter < fedrawdata.size() / sizeof(uint32_t) ) {
	fedrawdata_u32[iter] = temp_u32[iter+1];
	fedrawdata_u32[iter+1] = temp_u32[iter];
	iter+=2;
      }
      
      //@@ OVERWRITE HEADER AND TRAILER FOR DAQ
      FEDHeader header( fedrawdata.data() );
      header.set( fedrawdata.data(), 0, ( 0xFFFFFF && event.id().event() ), 0, *ifed ); //@@ LIMIT LV1 TO 24-BITS!!!
      
      FEDTrailer trailer( fedrawdata.data() + ( fedrawdata.size() - 8 )  );
      trailer.set( fedrawdata.data() + ( fedrawdata.size() - 8 ), 
		   fedrawdata.size() / 8,
		   evf::compute_crc( fedrawdata.data(), fedrawdata.size() ), 0, 0 );
      
    }
  }
  catch ( std::string err ) {
    if ( edm::isDebugEnabled() ) {
      edm::LogWarning("DigiToRaw") 
	<< "OldSiStripDigiToRaw::createFedBuffers] " 
	<< "Exception caught : " << err;
    }
  }
  
}


