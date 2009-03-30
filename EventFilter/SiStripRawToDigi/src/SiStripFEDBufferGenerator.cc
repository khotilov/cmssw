#include "EventFilter/SiStripRawToDigi/interface/SiStripFEDBufferGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <cstring>
#include <stdexcept>

namespace sistrip {
  
  //FEDStripData
  
  FEDStripData::FEDStripData(bool dataIsAlreadyConvertedTo8Bit, const size_t samplesPerChannel)
    : data_(FEDCH_PER_FED,ChannelData(dataIsAlreadyConvertedTo8Bit,samplesPerChannel))
  {
    if (samplesPerChannel > SCOPE_MODE_MAX_SCOPE_LENGTH) {
      std::ostringstream ss;
      ss << "Scope length " << samplesPerChannel << " is too long. "
         << "Max scope length is " << SCOPE_MODE_MAX_SCOPE_LENGTH << ".";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  const FEDStripData::ChannelData& FEDStripData::channel(const uint8_t internalFEDChannelNum) const
  {
    try {
      return data_.at(internalFEDChannelNum);
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Channel index out of range. (" << uint16_t(internalFEDChannelNum) << ") "
         << "Index should be in internal numbering scheme (0-95). ";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  uint16_t FEDStripData::ChannelData::getSample(const uint16_t sampleNumber) const
  {
    try {
      return data_.at(sampleNumber);
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Sample index out of range. "
         << "Requesting sample " << sampleNumber
         << " when channel has only " << data_.size() << " samples.";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  uint8_t FEDStripData::ChannelData::get8BitSample(const uint16_t sampleNumber) const
  {
    if (dataIs8Bit_) return (0xFF & getSample(sampleNumber));
    else {
      const uint16_t sample = getSample(sampleNumber);
      if (sample < 0xFE) return sample;
      else if (sample == 0x3FF) return 0xFF;
      else return 0xFE;
    }
  }
  
  void FEDStripData::ChannelData::setSample(const uint16_t sampleNumber, const uint16_t value)
  {
    if (value > 0x3FF) {
      std::ostringstream ss;
      ss << "Sample value (" << value << ") is too large. Maximum allowed is 1023. ";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
    try {
      data_.at(sampleNumber) = value;
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Sample index out of range. "
         << "Requesting sample " << sampleNumber
         << " when channel has only " << data_.size() << " samples.";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  //FEDBufferPayload
  
  FEDBufferPayload::FEDBufferPayload(const std::vector< std::vector<uint8_t> >& channelBuffers)
  {
    //calculate size of buffer and allocate enough memory
    uint32_t totalSize = 0;
    for (uint8_t iFE = 0; iFE < FEUNITS_PER_FED; iFE++) {
      for (uint8_t iCh = 0; iCh < FEDCH_PER_FEUNIT; iCh++) {
        totalSize += channelBuffers[iFE*FEDCH_PER_FEUNIT+iCh].size();
      }
      //if it does not finish on a 64Bit word boundary then take into account padding
      if (totalSize%8) {
        totalSize = ((totalSize/8) + 1)*8;
      }
    }
    data_.resize(totalSize);
    size_t indexInBuffer = 0;
    feLengths_.reserve(FEUNITS_PER_FED);
    //copy channel data into buffer with padding and update lengths
    for (uint8_t iFE = 0; iFE < FEUNITS_PER_FED; iFE++) {
      const size_t lengthAtStartOfFEUnit = indexInBuffer;
      //insert data for FE unit
      for (uint8_t iCh = 0; iCh < FEDCH_PER_FEUNIT; iCh++) {
        appendToBuffer(&indexInBuffer,channelBuffers[iFE*FEDCH_PER_FEUNIT+iCh].begin(),channelBuffers[iFE*FEDCH_PER_FEUNIT+iCh].end());
      }
      //store length
      feLengths_.push_back(indexInBuffer-lengthAtStartOfFEUnit);
      //add padding
      while (indexInBuffer % 8) appendToBuffer(&indexInBuffer,0);
    }
  }
  
  const uint8_t* FEDBufferPayload::data() const
  {
    //vectors are guarenteed to be contiguous
    if (lengthInBytes()) return &data_[0];
    //return NULL if there is no data yet
    else return NULL;
  }
  
  uint16_t FEDBufferPayload::getFELength(const uint8_t internalFEUnitNum) const
  {
    try{
      return feLengths_.at(internalFEUnitNum);
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Invalid FE unit number " << internalFEUnitNum << ". "
         << "Number should be in internal numbering scheme (0-7). ";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  FEDBufferPayload FEDBufferPayloadCreator::createPayload(const FEDReadoutMode mode, const FEDStripData& data) const
  {
    std::vector< std::vector<uint8_t> > channelBuffers(FEDCH_PER_FED,std::vector<uint8_t>());
    for (size_t iCh = 0; iCh < FEDCH_PER_FED; iCh++) {
      if (!feUnitsEnabled_[iCh/FEDCH_PER_FEUNIT]) continue;
      fillChannelBuffer(&channelBuffers[iCh],mode,data.channel(iCh),channelsEnabled_[iCh]);
    }
    return FEDBufferPayload(channelBuffers);
  }
  
  void FEDBufferPayloadCreator::fillChannelBuffer(std::vector<uint8_t>* channelBuffer, const FEDReadoutMode mode,
                                                 const FEDStripData::ChannelData& data, const bool channelEnabled) const
  {
    switch (mode) {
    case READOUT_MODE_SCOPE:
      fillRawChannelBuffer(channelBuffer,PACKET_CODE_SCOPE,data,channelEnabled,false);
      break;
    case READOUT_MODE_VIRGIN_RAW:
      fillRawChannelBuffer(channelBuffer,PACKET_CODE_VIRGIN_RAW,data,channelEnabled,true);
      break;
    case READOUT_MODE_PROC_RAW:
      fillRawChannelBuffer(channelBuffer,PACKET_CODE_PROC_RAW,data,channelEnabled,true);
      break;
    case READOUT_MODE_ZERO_SUPPRESSED:
      fillZeroSuppressedChannelBuffer(channelBuffer,data,channelEnabled);
      break;
    case READOUT_MODE_ZERO_SUPPRESSED_LITE:
      fillZeroSuppressedLiteChannelBuffer(channelBuffer,data,channelEnabled);
      break;
    default:
      std::ostringstream ss;
      ss << "Invalid readout mode " << mode;
      throw cms::Exception("FEDBufferGenerator") << ss.str();
      break;
    }
  }
  
  void FEDBufferPayloadCreator::fillRawChannelBuffer(std::vector<uint8_t>* channelBuffer,
                                                    const uint8_t packetCode,
                                                    const FEDStripData::ChannelData& data,
                                                    const bool channelEnabled,
                                                    const bool reorderData) const
  {
    const uint16_t channelLength = data.size();
    //2 bytes per sample + packet code + 2 bytes for length
    channelBuffer->reserve(channelLength*2 + 3);
    //length (max length is 0xFFF)
    channelBuffer->push_back( channelLength & 0xFF );
    channelBuffer->push_back( (channelLength & 0xF00) >> 8 );
    //packet code
    channelBuffer->push_back(packetCode);
    //channel samples
    for (uint16_t sampleNumber = 0; sampleNumber < channelLength; sampleNumber++) {
      const uint16_t sampleIndex = ( reorderData ? FEDStripOrdering::physicalOrderForStripInChannel(sampleNumber) : sampleNumber );
      uint16_t sampleValue = (channelEnabled ? data.getSample(sampleIndex) : 0);
      channelBuffer->push_back(sampleValue & 0xFF);
      channelBuffer->push_back((sampleValue & 0x300) >> 8);
    }
  }
  
  void FEDBufferPayloadCreator::fillZeroSuppressedChannelBuffer(std::vector<uint8_t>* channelBuffer,
                                                               const FEDStripData::ChannelData& data,
                                                               const bool channelEnabled) const
  {
    channelBuffer->reserve(50);
    //if channel is disabled then create empty channel header and return
    if (!channelEnabled) {
      //min length 7
      channelBuffer->push_back(7);
      channelBuffer->push_back(0);
      //packet code
      channelBuffer->push_back(PACKET_CODE_ZERO_SUPPRESSED);
      //4 bytes of medians
      channelBuffer->insert(channelBuffer->end(),4,0);
      return;
    }
    //if channel is not empty
    //add space for channel length
    channelBuffer->push_back(0xFF); channelBuffer->push_back(0xFF);
    //packet code
    channelBuffer->push_back(PACKET_CODE_ZERO_SUPPRESSED);
    //add medians
    const std::pair<uint16_t,uint16_t> medians = data.getMedians();
    channelBuffer->push_back(medians.first & 0xFF);
    channelBuffer->push_back((medians.first & 0x300) >> 8);
    channelBuffer->push_back(medians.second & 0xFF);
    channelBuffer->push_back((medians.second & 0x300) >> 8);
    //clusters
    fillClusterData(channelBuffer,data);
    //set length
    const uint16_t length = channelBuffer->size();
    (*channelBuffer)[0] = (length & 0xFF);
    (*channelBuffer)[1] = ((length & 0x300) >> 8);
  }
  
  void FEDBufferPayloadCreator::fillZeroSuppressedLiteChannelBuffer(std::vector<uint8_t>* channelBuffer,
                                                                   const FEDStripData::ChannelData& data,
                                                                   const bool channelEnabled) const
  {
    //if channel is disabled then create empty channel header and return
    if (!channelEnabled) {
      //min length 2
      channelBuffer->push_back(2);
      channelBuffer->push_back(0);
      return;
    }
    //if channel is not empty
    //add space for channel length
    channelBuffer->push_back(0xFF); channelBuffer->push_back(0xFF);
    //clusters
    fillClusterData(channelBuffer,data);
    //set length
    const uint16_t length = channelBuffer->size();
    (*channelBuffer)[0] = (length & 0xFF);
    (*channelBuffer)[1] = ((length & 0x300) >> 8);
  }
  
  void FEDBufferPayloadCreator::fillClusterData(std::vector<uint8_t>* channelBuffer, const FEDStripData::ChannelData& data) const
  {
    //current cluster info
    uint8_t clusterAddress = 0;
    std::list<uint8_t> clusterADCCounts;
    //loop over samples creating clusters
    for (size_t strip = 0; strip < data.size(); strip++) {
      const uint8_t adc = data.get8BitSample(strip);
      //check if strip is compatible with previous cluster
      if ( adc && (strip!=STRIPS_PER_APV) ) {
        //if this is the first strip in the cluster then update the address
        if (!clusterADCCounts.size()) clusterAddress = strip;
        clusterADCCounts.push_back(adc);
      }
      //if strip is not compatible with old cluster then write the old one to the buffer and start a new one
      else {
        writeClusterToBuffer(channelBuffer,clusterAddress,clusterADCCounts);
        clusterADCCounts.clear();
        if (adc) {
          if (!clusterADCCounts.size()) clusterAddress = strip;
          clusterADCCounts.push_back(adc);
        }
      }
    }
    //write last cluster
    writeClusterToBuffer(channelBuffer,clusterAddress,clusterADCCounts);
  }
  
  void FEDBufferPayloadCreator::writeClusterToBuffer(std::vector<uint8_t>* buffer, const uint8_t address,
                                                     const std::list<uint8_t> adcCounts) const
  {
    if (adcCounts.size()) {
      buffer->push_back(address);
      buffer->push_back(adcCounts.size());
      buffer->insert(buffer->end(),adcCounts.begin(),adcCounts.end());
    }
  }
  
  //FEDBufferGenerator
  
  FEDBufferGenerator::FEDBufferGenerator(const uint32_t l1ID, const uint16_t bxID,
                                         const std::vector<bool>& feUnitsEnabled, const std::vector<bool>& channelsEnabled,
                                         const FEDReadoutMode readoutMode, const FEDHeaderType headerType, const FEDBufferFormat bufferFormat,
                                         const FEDDAQEventType evtType, const FEDDataType dataType)
    : defaultDAQHeader_(l1ID,bxID,0,evtType),
      defaultDAQTrailer_(0,0),
      defaultTrackerSpecialHeader_(bufferFormat,readoutMode,headerType,dataType),
      defaultFEHeader_(FEDFEHeader::newFEHeader(headerType)),
      feUnitsEnabled_(feUnitsEnabled),
      channelsEnabled_(channelsEnabled)
  {
    if (!defaultFEHeader_.get()) {
      std::ostringstream ss;
      ss << "Bad header format: " << headerType;
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  bool FEDBufferGenerator::getFEUnitEnabled(const uint8_t internalFEUnitNumber) const
  {
    try {
      return feUnitsEnabled_.at(internalFEUnitNumber);
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Invalid FE unit number " << internalFEUnitNumber << ". Should be in internal numbering scheme (0-7)";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  bool FEDBufferGenerator::getChannelEnabled(const uint8_t internalFEDChannelNumber) const
  {
    try {
      return channelsEnabled_.at(internalFEDChannelNumber);
    } catch (const std::out_of_range&) {
      
      std::ostringstream ss;
      ss << "Invalid channel number " << internalFEDChannelNumber << ". "
         << "Should be in internal numbering scheme (0-95)";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
  }
  
  FEDBufferGenerator& FEDBufferGenerator::setFEUnitEnable(const uint8_t internalFEUnitNumber, const bool enabled)
  {
    try {
      feUnitsEnabled_.at(internalFEUnitNumber) = enabled;
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Invalid FE unit number " << internalFEUnitNumber << ". "
         << "Should be in internal numbering scheme (0-7)";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
    return *this;
  }
  
  FEDBufferGenerator& FEDBufferGenerator::setChannelEnable(const uint8_t internalFEDChannelNumber, const bool enabled)
  {
    try {
      channelsEnabled_.at(internalFEDChannelNumber) = enabled;
    } catch (const std::out_of_range&) {
      std::ostringstream ss;
      ss << "Invalid channel number " << internalFEDChannelNumber << ". "
         <<"Should be in internal numbering scheme (0-95)";
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
    return *this;
  }
  
  FEDBufferGenerator& FEDBufferGenerator::setFEUnitEnables(const std::vector<bool>& feUnitEnables)
  {
    if (feUnitEnables.size() != FEUNITS_PER_FED) {
      std::ostringstream ss;
      ss << "Setting FE enable vector with vector which is the wrong size. Size is " << feUnitEnables.size()
         << " it must be " << FEUNITS_PER_FED << "." << std::endl;
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
    feUnitsEnabled_ = feUnitEnables;
    return *this;
  }
  
  FEDBufferGenerator& FEDBufferGenerator::setChannelEnables(const std::vector<bool>& channelEnables)
  {
    if (channelEnables.size() != FEDCH_PER_FED) {
      std::ostringstream ss;
      ss << "Setting FED channel enable vector with vector which is the wrong size. Size is " << channelEnables.size()
         << " it must be " << FEDCH_PER_FED << "." << std::endl;
      throw cms::Exception("FEDBufferGenerator") << ss.str();
    }
    channelsEnabled_ = channelEnables;
    return *this;
  }
  
  void FEDBufferGenerator::generateBuffer(FEDRawData* rawDataObject, const FEDStripData& data, const uint16_t sourceID) const
  {
    //deal with disabled FE units and channels properly (FE enables, status bits)
    TrackerSpecialHeader tkSpecialHeader(defaultTrackerSpecialHeader_);
    std::auto_ptr<FEDFEHeader> fedFeHeader(defaultFEHeader_->clone());
    for (uint8_t iFE = 0; iFE < FEUNITS_PER_FED; iFE++) {
      const bool enabled = feUnitsEnabled_[iFE];
      tkSpecialHeader.setFEEnableForFEUnit(iFE,enabled);
      if (!enabled) {
        for (uint8_t iFEUnitChannel = 0; iFEUnitChannel < FEDCH_PER_FEUNIT; iFEUnitChannel++) {
          fedFeHeader->setChannelStatus(iFE,iFEUnitChannel,FEDChannelStatus(0));
        }
      }
    }
    for (uint8_t iCh = 0; iCh < FEDCH_PER_FED; iCh++) {
      if (!channelsEnabled_[iCh]) {
        fedFeHeader->setChannelStatus(iCh,FEDChannelStatus(0));
      }
    }
    //set the source ID
    FEDDAQHeader daqHeader(defaultDAQHeader_);
    daqHeader.setSourceID(sourceID);
    //build payload
    const FEDBufferPayloadCreator payloadPacker(feUnitsEnabled_,channelsEnabled_);
    const FEDBufferPayload payload = payloadPacker(getReadoutMode(),data);
    //fill FE lengths
    for (uint8_t iFE = 0; iFE < FEUNITS_PER_FED; iFE++) {
      fedFeHeader->setFEUnitLength(iFE,payload.getFELength(iFE));
    }
    //resize buffer
    rawDataObject->resize(bufferSizeInBytes(*fedFeHeader,payload));
    //fill buffer
    fillBuffer(rawDataObject->data(),daqHeader,defaultDAQTrailer_,tkSpecialHeader,*fedFeHeader,payload);
  }
  
  void FEDBufferGenerator::fillBuffer(uint8_t* pointerToStartOfBuffer,
                                      const FEDDAQHeader& daqHeader,
                                      const FEDDAQTrailer& daqTrailer,
                                      const TrackerSpecialHeader& tkSpecialHeader,
                                      const FEDFEHeader& feHeader,
                                      const FEDBufferPayload& payload)
  {
    //set the length in the DAQ trailer
    const size_t lengthInBytes = bufferSizeInBytes(feHeader,payload);
    FEDDAQTrailer updatedDAQTrailer(daqTrailer);
    updatedDAQTrailer.setEventLengthIn64BitWords(lengthInBytes/8);
    //copy pieces into buffer in order
    uint8_t* bufferPointer = pointerToStartOfBuffer;
    memcpy(bufferPointer,daqHeader.data(),8);
    bufferPointer += 8;
    memcpy(bufferPointer,tkSpecialHeader.data(),8);
    bufferPointer += 8;
    memcpy(bufferPointer,feHeader.data(),feHeader.lengthInBytes());
    bufferPointer += feHeader.lengthInBytes();
    memcpy(bufferPointer,payload.data(),payload.lengthInBytes());
    bufferPointer += payload.lengthInBytes();
    memcpy(bufferPointer,updatedDAQTrailer.data(),8);
    //word swap if necessary
    if (tkSpecialHeader.wasSwapped()) {
      for (size_t i = 0; i < lengthInBytes; i++) {
        pointerToStartOfBuffer[i] = pointerToStartOfBuffer[i^4];
      }
    }
    //update CRC
    const uint16_t crc = calculateFEDBufferCRC(pointerToStartOfBuffer,lengthInBytes);
    updatedDAQTrailer.setCRC(crc);
    memcpy(bufferPointer,updatedDAQTrailer.data(),8);
    if (tkSpecialHeader.wasSwapped()) {
      for (size_t i = 0; i < 8; i++) {
        bufferPointer[i] = bufferPointer[i^4];
      }
    }
  }
  
}
