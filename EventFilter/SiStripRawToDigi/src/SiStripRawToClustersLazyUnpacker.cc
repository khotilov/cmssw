#include "EventFilter/SiStripRawToDigi/interface/SiStripRawToClustersLazyUnpacker.h"
#include <sstream>
#include <iostream>

using namespace sistrip;

SiStripRawToClustersLazyUnpacker::SiStripRawToClustersLazyUnpacker(const SiStripRegionCabling& regioncabling, const SiStripClusterizerFactory& clustfact, const FEDRawDataCollection& data, bool dump) :

  raw_(&data),
  regions_(&(regioncabling.getRegionCabling())),
  clusterizer_(&clustfact),
  buffers_(),
  rawToDigi_(0,0,0,0,0),
  dump_(dump),
  mode_(sistrip::READOUT_MODE_INVALID),
  fedRawData_()
{
  buffers_.assign(1024,static_cast<sistrip::FEDBuffer*>(0));
}

SiStripRawToClustersLazyUnpacker::~SiStripRawToClustersLazyUnpacker() {

  std::vector< sistrip::FEDBuffer*>::iterator ibuffer = buffers_.begin();
  for (; ibuffer!=buffers_.end(); ibuffer++) {
    if (*ibuffer) delete *ibuffer;
  }
}

void SiStripRawToClustersLazyUnpacker::fill(const uint32_t& index, record_type& record) {

  // Get region, subdet and layer from element-index
  uint32_t region = SiStripRegionCabling::region(index);
  uint32_t subdet = static_cast<uint32_t>(SiStripRegionCabling::subdet(index));
  uint32_t layer = SiStripRegionCabling::layer(index);
 
  // Retrieve cabling for element
  const SiStripRegionCabling::ElementCabling& element = (*regions_)[region][subdet][layer];
  
  // Loop dets
  SiStripRegionCabling::ElementCabling::const_iterator idet = element.begin();
  for (;idet!=element.end();idet++) {
    
    // If det id is null or invalid continue.
    if ( !(idet->first) || (idet->first == sistrip::invalid32_) ) { continue; }
    
    // Loop over apv-pairs of det
    std::vector<FedChannelConnection>::const_iterator iconn = idet->second.begin();
    for (;iconn!=idet->second.end();iconn++) {
      
      // If fed id is null or connection is invalid continue
      if ( !iconn->fedId() || !iconn->isConnected() ) { continue; }    
      
      // If Fed hasnt already been initialised, extract data and initialise
      if (!buffers_[iconn->fedId()]) {
	
	// Retrieve FED raw data for given FED
	const FEDRawData& input = raw_->FEDData( static_cast<int>(iconn->fedId()) );
	
	// Cache new correctly-ordered FEDRawData object (to maintain scope for Fed9UEvent)
	fedRawData_.push_back( FEDRawData() );
	rawToDigi_.locateStartOfFedBuffer( iconn->fedId(), input, fedRawData_.back() );
	
	// Check on FEDRawData pointer
	if ( !fedRawData_.back().data() ) {
	  edm::LogWarning(mlRawToCluster_)
	    << "[SiStripRawToClustersLazyGetter::" 
	    << __func__ 
	    << "]"
	    << " NULL pointer to FEDRawData for FED id " 
	    << iconn->fedId();
	  continue;
	}	
	
	// Check on FEDRawData size
	if ( !fedRawData_.back().size() ) {
	  edm::LogWarning(mlRawToCluster_)
	    << "[SiStripRawToClustersLazyGetter::" 
	    << __func__ << "]"
	    << " FEDRawData has zero size for FED id " 
	    << iconn->fedId();
	  continue;
	}
	
	// construct FEDBuffer
	try {buffers_[iconn->fedId()] = new sistrip::FEDBuffer(fedRawData_.back().data(),fedRawData_.back().size());}
	catch (const cms::Exception& e) { 
	  edm::LogWarning(mlRawToCluster_) 
	    << e.what();
	  if ( buffers_[iconn->fedId()] ) { delete buffers_[iconn->fedId()]; }
	  buffers_[iconn->fedId()] = 0;
	  continue;
	}

	// dump of FEDRawData to stdout
	if ( dump_ ) {
	  std::stringstream ss;
	  rawToDigi_.dumpRawData( iconn->fedId(), input, ss );
	  LogTrace(mlRawToDigi_) 
	    << ss.str();
	}

	// record readout mode
	mode_ = buffers_[iconn->fedId()]->readoutMode();
      }

      // check channel
      if (!buffers_[iconn->fedId()]->channelGood(iconn->fedCh())) continue; 
      
      // Determine key that should be used to index digi containers
      uint32_t key = iconn->detId();
      
      // Determine APV std::pair number
      uint16_t ipair = iconn->apvPairNumber();


      if (mode_ == sistrip::READOUT_MODE_ZERO_SUPPRESSED ) { 
	
	// create unpacker
	sistrip::FEDZSChannelUnpacker unpacker = sistrip::FEDZSChannelUnpacker::zeroSuppressedModeUnpacker(buffers_[iconn->fedId()]->channel(iconn->fedCh()));
	
	// unpack
	while (unpacker.hasData()) {
	  clusterizer_->algorithm()->add(record,key,unpacker.strip()+ipair*256,unpacker.adc());
	  unpacker++;
	}
      }

      else if (mode_ == sistrip::READOUT_MODE_ZERO_SUPPRESSED_LITE ) { 
	
	// create unpacker
	sistrip::FEDZSChannelUnpacker unpacker = sistrip::FEDZSChannelUnpacker::zeroSuppressedLiteModeUnpacker(buffers_[iconn->fedId()]->channel(iconn->fedCh()));
	
	// unpack
	while (unpacker.hasData()) {
	  clusterizer_->algorithm()->add(record,key,unpacker.strip()+ipair*256,unpacker.adc());
	  unpacker++;
	}
      }

      else {
	  edm::LogWarning(mlRawToCluster_)
	    << "[SiStripRawToClustersLazyGetter::" 
	    << __func__ << "]"
	    << " FEDRawData readout mode "
	    << mode_
	    << " from FED id "
	    << iconn->fedId() 
	    << " not supported."; 
	  continue;
	}
    }
    clusterizer_->algorithm()->endDet(record,idet->first);
  }
}

sistrip::RawToClustersLazyUnpacker::RawToClustersLazyUnpacker(const SiStripRegionCabling& regioncabling, const SiStripClusterizerFactory& clustfact, const FEDRawDataCollection& data) :

  raw_(&data),
  regions_(&(regioncabling.getRegionCabling())),
  clusterizer_(&clustfact),
  fedEvents_(),
  rawToDigi_(0,0,0,0,0),
  fedRawData_()
{
  fedEvents_.assign(1024,static_cast<Fed9U::Fed9UEvent*>(0));
}

sistrip::RawToClustersLazyUnpacker::~RawToClustersLazyUnpacker() {
  std::vector< Fed9U::Fed9UEvent*>::iterator ifedevent = fedEvents_.begin();
  for (; ifedevent!=fedEvents_.end(); ifedevent++) {
    if (*ifedevent) {
      delete (*ifedevent);
      *ifedevent = 0;
    }
  }
}

void sistrip::RawToClustersLazyUnpacker::fill(const uint32_t& index, record_type& record) {

  //Get region, subdet and layer from element-index
  uint32_t region = SiStripRegionCabling::region(index);
  uint32_t subdet = static_cast<uint32_t>(SiStripRegionCabling::subdet(index));
  uint32_t layer = SiStripRegionCabling::layer(index);
 
  //Retrieve cabling for element
  const SiStripRegionCabling::ElementCabling& element = (*regions_)[region][subdet][layer];
  
  //Loop dets
  SiStripRegionCabling::ElementCabling::const_iterator idet = element.begin();
  for (;idet!=element.end();idet++) {
    
    //If det id is null or invalid continue.
    if ( !(idet->first) || (idet->first == sistrip::invalid32_) ) { continue; }
    
    //Loop over apv-pairs of det
    std::vector<FedChannelConnection>::const_iterator iconn = idet->second.begin();
    for (;iconn!=idet->second.end();iconn++) {
      
      //If fed id is null or connection is invalid continue
      if ( !iconn->fedId() || !iconn->isConnected() ) { continue; }    
      
      //If Fed hasnt already been initialised, extract data and initialise
      if (!fedEvents_[iconn->fedId()]) {
	
	// Retrieve FED raw data for given FED
	const FEDRawData& input = raw_->FEDData( static_cast<int>(iconn->fedId()) );
	
	// Cache new correctly-ordered FEDRawData object (to maintain scope for Fed9UEvent)
	fedRawData_.push_back( FEDRawData() );
	rawToDigi_.locateStartOfFedBuffer( iconn->fedId(), input, fedRawData_.back() );

	// Recast data to suit Fed9UEvent
	Fed9U::u32* data_u32 = reinterpret_cast<Fed9U::u32*>( const_cast<unsigned char*>( fedRawData_.back().data() ) );
	Fed9U::u32  size_u32 = static_cast<Fed9U::u32>( fedRawData_.back().size() / 4 ); 
	
	// Check on FEDRawData pointer
	if ( !data_u32 ) {
	  if ( edm::isDebugEnabled() ) {
	    edm::LogWarning(mlRawToCluster_)
	      << "[SiStripRawToClustersLazyGetter::" 
	      << __func__ 
	      << "]"
	      << " NULL pointer to FEDRawData for FED id " 
	      << iconn->fedId();
	  }
	  continue;
	}	
	
	// Check on FEDRawData size
	if ( !size_u32 ) {
	  if ( edm::isDebugEnabled() ) {
	    edm::LogWarning(mlRawToCluster_)
	      << "[SiStripRawToClustersLazyGetter::" 
	      << __func__ << "]"
	      << " FEDRawData has zero size for FED id " 
	      << iconn->fedId();
	  }
	  continue;
	}
	
	// Construct Fed9UEvent using present FED buffer
	try {fedEvents_[iconn->fedId()] = new Fed9U::Fed9UEvent(data_u32,0,size_u32);} 
	catch(...) { 
	  rawToDigi_.handleException( __func__, "Problem when constructing Fed9UEvent" ); 
	  if ( fedEvents_[iconn->fedId()] ) { delete fedEvents_[iconn->fedId()]; }
	  fedEvents_[iconn->fedId()] = 0;
	  continue;
	}
      }
      
      //Calculate corresponding FED unit, channel
      uint16_t iunit = 0;
      uint16_t ichan = 0;
      uint16_t chan = 0;
      try {
	Fed9U::Fed9UAddress addr;
	addr.setFedChannel( static_cast<unsigned char>( iconn->fedCh() ) );
        //0-7 (internal)
	iunit = addr.getFedFeUnit();
	//0-11 (internal)
	ichan = addr.getFeUnitChannel();
	//0-95 (internal)
	chan = 12*( iunit ) + ichan;
      } catch(...) { 
	rawToDigi_.handleException( __func__, "Problem using Fed9UAddress" ); 
      } 
      
      try {
	// temporary check to find corrupted FED data
        #ifdef USE_PATCH_TO_CATCH_CORRUPT_FED_DATA
	uint16_t last_strip = 0;
	uint16_t strips = 256 * iconn->nApvPairs();
        #endif
	Fed9U::Fed9UEventIterator fed_iter = const_cast<Fed9U::Fed9UEventChannel&>(fedEvents_[iconn->fedId()]->channel( iunit, ichan )).getIterator();
	for (Fed9U::Fed9UEventIterator i = fed_iter+7; i.size() > 0;) {
	  uint16_t first_strip = iconn->apvPairNumber()*256 + *i++;
	  unsigned char width = *i++; 
	  for ( uint16_t istr = 0; istr < ((uint16_t)width); istr++) {
	    uint16_t strip = first_strip + istr;
	    // temporary check to find corrupted FED data
            #ifdef USE_PATCH_TO_CATCH_CORRUPT_FED_DATA
	    if ( !( strip < strips && ( !strip || strip > last_strip ) ) ) { 
	      edm::LogWarning(mlRawToDigi_)
		<< "[SiStripRawToClustersLazyUnpacker::" 
		<< __func__ << "]"
		<< " Corrupt FED data found for FED id " 
		<< iconn->fedId()
		<< " and channel " 
		<< iconn->fedCh()
		<< "!  present strip: " 
		<< strip
		<< "  last strip: " 
		<< last_strip
		<< "  detector strips: " 
		<< strips;
	      continue; 
	    }
	    last_strip = strip;
            #endif	 
	    clusterizer_->algorithm()->add(record,idet->first,strip,(uint16_t)(*i++));
	  }
	}
      } catch(...) { 
	std::stringstream sss;
	sss << "Problem accessing ZERO_SUPPR data for FED id/ch: " 
	    << iconn->fedId() 
	    << "/" 
	    << chan;
	rawToDigi_.handleException( __func__, sss.str() ); 
      } 
    }
    clusterizer_->algorithm()->endDet(record,idet->first);
  }
}
  
