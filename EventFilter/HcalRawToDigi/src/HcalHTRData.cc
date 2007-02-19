//#include "Utilities/Configuration/interface/Architecture.h"
/*  
 *  $Date: 2007/02/19 23:26:44 $
 *  $Revision: 1.5 $
 *  \author J. Mans -- UMD
 */
#ifndef HTBDAQ_DATA_STANDALONE
#include "EventFilter/HcalRawToDigi/interface/HcalHTRData.h"
#else
#include "HcalHTRData.h"
#endif
#include <string.h>
#include <stdio.h>
const int HcalHTRData::CHANNELS_PER_SPIGOT         = 24;
const int HcalHTRData::MAXIMUM_SAMPLES_PER_CHANNEL = 20;



HcalHTRData::HcalHTRData() : m_formatVersion(-2), m_rawLength(0), m_rawConst(0), m_ownData(0) { }
HcalHTRData::HcalHTRData(const unsigned short* data, int length) {
  adoptData(data,length);
  m_ownData=0;
}
HcalHTRData::HcalHTRData(const HcalHTRData& hd) : m_formatVersion(hd.m_formatVersion), m_rawLength(hd.m_rawLength), m_rawConst(hd.m_rawConst), m_ownData(0) { }

HcalHTRData::HcalHTRData(int version_to_create) : m_formatVersion(version_to_create) {
  allocate(version_to_create);
}

void HcalHTRData::allocate(int version_to_create) {
  m_formatVersion=version_to_create;
  // the needed space is for the biggest possible event...
  const int needed=0x200;
  // create a buffer big enough...
  m_ownData=new unsigned short[needed];
  m_rawLength=0;
  m_rawConst=m_ownData;
}

HcalHTRData& HcalHTRData::operator=(const HcalHTRData& hd) {
  if (m_ownData==0) {
    m_formatVersion=hd.m_formatVersion;
    m_rawLength=hd.m_rawLength;
    m_rawConst=hd.m_rawConst;
  }
  return (*this);
}

void HcalHTRData::adoptData(const unsigned short* data, int length) {
  m_rawLength=length;
  m_rawConst=data;
  if (m_rawLength<5) {
    m_formatVersion=-2; // invalid!
  } else {
    // determine format version
    if ((m_rawConst[2]&0x8000)==0) m_formatVersion=-1; // original format before versions
    else m_formatVersion=(m_rawConst[4]>>12)&0xF;
  }
}

// check :: not EE, length is reasonable, length matches wordcount
//          length required for tp+daq is correct

bool HcalHTRData::check() const {
  if (m_formatVersion==-1) {
    // length checks
    //  minimum length
    if (m_rawLength<6+12) return false;
    //  matches wordcount
    if (m_rawLength!=m_rawConst[m_rawLength-3]) return false;
    // empty event check
    if (m_rawConst[2]&0x20) return false;
  } else {
    // length checks
    //  minimum length
    if (m_rawLength<8+4) return false;
    //  matches wordcount
    if (m_rawLength!=m_rawConst[m_rawLength-3]) {
      if (isHistogramEvent() && m_rawConst[m_rawLength-3]==786) {
	// known bug!
      } else
	return false;
    }
    // empty event check (redundant...)
    if (m_rawConst[2]&0x4) return false;
  }

  if (!isHistogramEvent()) {
    // daq/tp length check
    int tp, daq, header, trailer;
    determineSectionLengths(tp,daq,header,trailer);
    if (tp+daq+header+trailer>m_rawLength) return false;
  }

  return true;
}

void HcalHTRData::determineSectionLengths(int& tpWords, int& daqWords, int& headerWords, int& trailerWords) const {
  if (m_formatVersion==-1) {
    tpWords=m_rawConst[5]>>8;
    daqWords=CHANNELS_PER_SPIGOT*(m_rawConst[m_rawLength-4]>>8); // always 24 channels, no zero suppresion
    headerWords=6;
    trailerWords=12;
  } else {
    tpWords=m_rawConst[5]>>8;
    if (m_rawLength>4) 
      daqWords=m_rawConst[m_rawLength-4]&0x7FF; // zero suppression supported
    headerWords=8;
    trailerWords=4; // minimum, may be more...
  }
}

void HcalHTRData::determineStaticLengths(int& headerWords, int& trailerWords) const {
  if (m_formatVersion==-1) {
    headerWords=6;
    trailerWords=12;
  } else {
    headerWords=8;
    trailerWords=4; // minimum, may be more...
  }
}

void HcalHTRData::dataPointers(const unsigned short** daq_first, 
			       const unsigned short** daq_last, 
			       const unsigned short** tp_first, 
			       const unsigned short** tp_last) {
  int tp_words_total, daq_words_total, headerLen, trailerLen;
  determineSectionLengths(tp_words_total,daq_words_total,headerLen,trailerLen);
  
  *tp_first=m_rawConst+headerLen;
  *tp_last=*tp_first+(tp_words_total-1);
  *daq_first=*tp_last+1;
  *daq_last=*daq_first+(daq_words_total-1);
}

/* using FiberAd[2:0] ChanId[1:0] */
static const int channelDecoder[32] = { 0, 1, 2, 99, 3, 4, 5, 99, 
                                        6, 7, 8, 99, 9,10,11, 99,
                                        12,13,14,99,15,16,17, 99,
                                        18,19,20,99,21,22,23, 99};

void HcalHTRData::unpack(unsigned char* daq_lengths, unsigned short* daq_samples,
			 unsigned char* tp_lengths, unsigned short* tp_samples) const {

  if (daq_lengths!=0) memset(daq_lengths,0,CHANNELS_PER_SPIGOT);
  if (tp_lengths!=0) memset(tp_lengths,0,CHANNELS_PER_SPIGOT);

  // currently, the major differences between the versions are
  //  -1 : 6 word header, no zero suppression, trailer setup
  //   0 : 8 word header, zero suppression, 

  int tp_words_total, daq_words_total, headerLen, trailerLen;
  determineSectionLengths(tp_words_total,daq_words_total,headerLen,trailerLen);

  //  printf("%d %d %d %d\n",tp_words_total,daq_words_total,headerLen,trailerLen);
  int wordPtr;
  const unsigned short* tpBase=m_rawConst+headerLen;
  // process the trigger primitive words
  if (tp_lengths!=0) {
    for (wordPtr=0; wordPtr<tp_words_total; wordPtr++) {
      int ichan=channelDecoder[tpBase[wordPtr]>>11];
      if (ichan>=24) continue;
      tp_samples[ichan*MAXIMUM_SAMPLES_PER_CHANNEL+tp_lengths[ichan]]=tpBase[wordPtr]&0x3ff;
      tp_lengths[ichan]++;
    }
  }
 
  const unsigned short* daqBase=m_rawConst+headerLen+tp_words_total;
  // process the DAQ words [ assumes that data from one channel will always be together ]
  int lastChan=-1;
  int lastCapid=0;
  if (daq_lengths!=0) {
    for (wordPtr=0; wordPtr<daq_words_total; wordPtr++) {
      int ichan=channelDecoder[daqBase[wordPtr]>>11];
      if (ichan>=24) continue;
      int capid=(daqBase[wordPtr]&0x180)>>7;
      int erdv=(daqBase[wordPtr]&0x600)>>9;
      if (erdv!=0x1 || 
	  (lastChan==ichan && (capid!=((lastCapid+1)%4)))) {
	daq_lengths[ichan]|=0x80;
      } 
      lastChan=ichan;
      lastCapid=capid;

      int useLength=daq_lengths[ichan]&0x1F;
      //     printf("%d %d\n",ichan,useLength);
      daq_samples[ichan*MAXIMUM_SAMPLES_PER_CHANNEL+useLength]=daqBase[wordPtr]&0x3ff;
      daq_lengths[ichan]=(useLength+1)|(daq_lengths[ichan]&0xE0); // keep the error bits
    }
  }

}

void HcalHTRData::pack(unsigned char* daq_lengths, unsigned short* daq_samples,
		       unsigned char* tp_lengths, unsigned short* tp_samples, bool do_capid) {
  
  int tp_words_total=0, daq_words_total=0, headerLen, trailerLen;
  determineStaticLengths(headerLen,trailerLen);

  tp_words_total=0;
  daq_words_total=0;
  int ichan,isample;

  // trigger primitive words
  unsigned short* ptr=m_ownData+headerLen;
  if (tp_samples!=0 && tp_lengths!=0) {
    for (ichan=0; ichan<24; ichan++) {
      unsigned short chanid=((ichan%3)+((ichan/3)<<2))<<11;
      for (isample=0; isample<tp_lengths[ichan] && isample<MAXIMUM_SAMPLES_PER_CHANNEL; isample++) {
	ptr[tp_words_total]=chanid|(tp_samples[ichan*MAXIMUM_SAMPLES_PER_CHANNEL+isample]&0x3FF);
	tp_words_total++;
      }
    }
  }

  // daq words
  ptr=m_ownData+headerLen+tp_words_total;
  for (ichan=0; ichan<24; ichan++) {
    unsigned short chanid=((ichan%3)+((ichan/3)<<2))<<11;
    for (isample=0; isample<daq_lengths[ichan] && isample<MAXIMUM_SAMPLES_PER_CHANNEL; isample++) {
      unsigned short basedata=daq_samples[ichan*MAXIMUM_SAMPLES_PER_CHANNEL+isample]&0x3FF;
      if (do_capid) basedata=(basedata&0x7F)|(0x200)|((isample%4)<<7);
      ptr[daq_words_total]=chanid|basedata;
      daq_words_total++;
    }
  }

  if (m_formatVersion==-1) {
    m_ownData[5]=(tp_words_total<<8)|0x1;
    unsigned short totalLen=headerLen+tp_words_total+daq_words_total+trailerLen;
    m_rawLength=totalLen;
    m_ownData[totalLen-3]=totalLen;
    m_ownData[totalLen-4]=(tp_words_total/CHANNELS_PER_SPIGOT)|((daq_words_total/CHANNELS_PER_SPIGOT)<<8);
  } else {
    m_ownData[5]=(tp_words_total<<8)|0x1;
    unsigned short totalLen=headerLen+tp_words_total+daq_words_total+trailerLen;
    if ((totalLen%2)==1) {
      m_ownData[totalLen-4]=0xFFFF; // parity word
      totalLen++; // round to even number of 16-bit words
    }
    m_rawLength=totalLen;
    m_ownData[totalLen-2]=totalLen/2; // 32-bit words
    m_ownData[totalLen-3]=totalLen;
    m_ownData[totalLen-4]=daq_words_total;
  }

}

void HcalHTRData::packHeaderTrailer(int L1Anumber, int bcn, int submodule, int orbitn, int pipeline, int ndd, int nps, int firmwareRev) {
  m_ownData[0]=L1Anumber&0xFF;
  m_ownData[1]=(L1Anumber&0xFFFF00)>>8;
  if (m_formatVersion==-1) {
    m_ownData[2]=((pipeline&0x7F)<<8); // no error bits
    m_ownData[3]=((orbitn&0xFF)<<8)|(submodule&0xFF);
    m_ownData[4]=bcn&0xFFF;
    //    m_ownData[5]&=0xFF01;
  } else {
    m_ownData[2]=0x8000; // Version is valid, no error bits
    m_ownData[3]=((orbitn&0x3F)<<10)|(submodule&0x3FF);
    m_ownData[4]=((m_formatVersion&0xF)<<12)|(bcn&0xFFF);
    m_ownData[5]|=((nps&0x1F)<<3)|0x1;
    m_ownData[6]=((firmwareRev&0x70000)>>3)|(firmwareRev&0x1FFF);
    m_ownData[7]=pipeline&0xFF;
    m_ownData[m_rawLength-4]&=0x7FF;
    m_ownData[m_rawLength-4]|=(ndd&0x1F)<<11;
  }
  m_ownData[m_rawLength-2]=m_rawLength/2; // 32-bit words
  m_ownData[m_rawLength-1]=(L1Anumber&0xFF)<<8;
}

unsigned int HcalHTRData::getOrbitNumber() const { 
  return (m_formatVersion==-1)?(m_rawConst[3]>>8):(m_rawConst[3]>>10);
}
unsigned int HcalHTRData::getSubmodule() const {
  return (m_formatVersion==-1)?(m_rawConst[3]&0xFF):(m_rawConst[3]&0x3FF);
}
unsigned int HcalHTRData::htrSlot() const{
  const unsigned int smid = getSubmodule();
  return ((smid>>1)&0x1F);
} 
unsigned int HcalHTRData::htrTopBottom() const{
  const unsigned int smid = getSubmodule();
  return (smid&0x01);
} 
unsigned int HcalHTRData::readoutVMECrateId() const{
  const unsigned int smid = getSubmodule();
  return ((smid>>6)&0x1F);
} 
bool HcalHTRData::isCalibrationStream() const {
  return (m_formatVersion==-1)?(false):(m_rawConst[2]&0x4000);
}
bool HcalHTRData::isPatternRAMEvent() const {
  return (m_formatVersion==-1)?(false):(m_rawConst[2]&0x1000);
}
bool HcalHTRData::isHistogramEvent() const {
  return (m_formatVersion==-1)?(m_rawConst[2]&0x2):(m_rawConst[2]&0x2000);
}
int HcalHTRData::getNDD() const {
  return (m_formatVersion==-1)?(m_rawConst[m_rawLength-4]>>8):(m_rawConst[m_rawLength-4]>>11);
}
int HcalHTRData::getNTP() const {
  return (m_formatVersion==-1)?(m_rawConst[m_rawLength-4]&0xFF):(m_rawConst[m_rawLength-4]>>11);
}
int HcalHTRData::getNPS() const {
  return (m_formatVersion==-1)?(0):((m_rawConst[5]>>3)&0x1F);
}
unsigned int HcalHTRData::getPipelineLength() const {
  return (m_formatVersion==-1)?(m_rawConst[2]>>8):(m_rawConst[7]&0xFF);
}
unsigned int HcalHTRData::getFirmwareRevision() const {
  return (m_formatVersion==-1)?(0):((m_rawConst[6]&0x1FFF)+((m_rawConst[6]&0xE000)<<3));
}

void HcalHTRData::getHistogramFibers(int& a, int& b) const {
  a=-1;
  b=-1;
  if (m_formatVersion==-1) {
    a=((m_rawConst[2]&0x0F00)>>8);
    b=((m_rawConst[2]&0xF000)>>12);
  } else {
    a=((m_rawConst[5]&0x0F00)>>8);
    b=((m_rawConst[5]&0xF000)>>12);
  }
}

bool HcalHTRData::wasHistogramError(int ifiber) const {
  bool retval=!isHistogramEvent();
  if (!retval) {
    retval=((m_rawConst[7])&(1<<ifiber))!=0;
  }
  return retval;
}

bool HcalHTRData::unpackHistogram(int myfiber, int mysc, int capid, unsigned short* histogram) const {
  // check for histogram mode
  if (!isHistogramEvent()) return false;

  int fiber1, fiber2;
  getHistogramFibers(fiber1,fiber2);
  if (fiber1!=myfiber && fiber2!=myfiber) return false;

  if (m_formatVersion==-1) {
    int offset=6+mysc*4*32+capid*32;
    if (myfiber==fiber2) offset+=3*4*32; // skip to the second half...
    for (int i=0; i<32; i++)
      histogram[i]=m_rawConst[offset+i];
    return true;
  } else {
    int offset=8+mysc*4*32+capid*32;
    if (myfiber==fiber2) offset+=3*4*32; // skip to the second half...
    for (int i=0; i<32; i++)
      histogram[i]=m_rawConst[offset+i];
    return true;
  }
}
