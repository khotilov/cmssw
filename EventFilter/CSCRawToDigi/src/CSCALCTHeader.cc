#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include <iomanip>

bool CSCALCTHeader::debug=false;
short unsigned int CSCALCTHeader::firmwareVersion=2006; 

CSCALCTHeader::CSCALCTHeader(int chamberType)
: header2006(chamberType)
{ //constructor for digi->raw packing based on header2006 
  memcpy(theOriginalBuffer, &header2006, header2006.sizeForPacking());

}

CSCALCTHeader::CSCALCTHeader(const unsigned short * buf) {
  ///collision and hot channel masks are variable sized
  ///the sizes vary depending on type of the ALCT board
  ///                                        number of words for various
  ///                                        alct board types:  1  2  3     5  6
  static unsigned short int collisionMaskWordcount[7]    = { 8, 8,12,16,16,24,28};
  static unsigned short int hotChannelMaskWordcount[7]   = {18,18,24,36,36,48,60};


  ///first determine the correct format  
  if (buf[0]==0xDB0A) {
    firmwareVersion=2007;
  }
  else if ( (buf[0]&0xF800)==0x6000 ) {
    firmwareVersion=2006;
  }
  else {
    edm::LogError("CSCALCTHeader|CSCRawToDigi") << "failed to determine ALCT firmware version!!";
  }

  LogTrace("CSCALCTHeader|CSCRawToDigi") << "firmware version - " << firmwareVersion;

  ///Now fill data 
  switch (firmwareVersion) {
  case 2006:
    memcpy(&header2006, buf, header2006.sizeInWords()*2);///the header part
    buf +=header2006.sizeInWords();
    memcpy(&alcts2006, buf, alcts2006.sizeInWords()*2);///the alct0 and alct1
    buf +=alcts2006.sizeInWords();
    break;

  case 2007:
    memcpy(&header2007, buf, header2007.sizeInWords()*2); ///the fixed sized header part
    buf +=header2007.sizeInWords();
    sizeInWords2007_ = header2007.sizeInWords();
    ///now come the variable parts
    if (header2007.configPresent==1) {
      memcpy(&virtexID, buf, virtexID.sizeInWords()*2);
      buf +=virtexID.sizeInWords();
      sizeInWords2007_ = virtexID.sizeInWords();
      memcpy(&configRegister, buf, configRegister.sizeInWords()*2);
      buf +=configRegister.sizeInWords();
      sizeInWords2007_ += configRegister.sizeInWords();
      
      collisionMasks.resize(collisionMaskWordcount[header2007.boardType]);
      for (unsigned int i=0; i<collisionMaskWordcount[header2007.boardType]; ++i){
	memcpy(&collisionMasks[i], buf, collisionMasks[i].sizeInWords()*2);
	buf += collisionMasks[i].sizeInWords();
	sizeInWords2007_ += collisionMasks[i].sizeInWords();
      }

      hotChannelMasks.resize(hotChannelMaskWordcount[header2007.boardType]);
      for (unsigned int i=0; i<hotChannelMaskWordcount[header2007.boardType]; ++i) {
	memcpy(&hotChannelMasks[i], buf, hotChannelMasks[i].sizeInWords()*2);
	buf += hotChannelMasks[i].sizeInWords();
	sizeInWords2007_ += hotChannelMasks[i].sizeInWords();
      }
    }

    alcts.resize(header2007.lctBins*2); ///2007 has LCTbins * 2 alct words
    for (unsigned int i=0; i<header2007.lctBins*2; ++i) {
      memcpy(&alcts[i], buf, alcts[i].sizeInWords()*2);
      buf += alcts[i].sizeInWords(); 
      sizeInWords2007_ += alcts[i].sizeInWords();
    }
    break;

  default:
    edm::LogError("CSCALCTHeader|CSCRawToDigi")
      <<"couldn't construct: ALCT firmware version is bad/not defined!";
    break;
  }

  ///also store raw data buffer too; it is later returned by data() method
  if ((firmwareVersion==2006)||(firmwareVersion==2007))
    memcpy(theOriginalBuffer, buf-sizeInWords(), sizeInWords()*2);
  
}


CSCALCTHeader::CSCALCTHeader(const CSCALCTStatusDigi & digi){
  CSCALCTHeader(digi.header());
}

void CSCALCTHeader::setEventInformation(const CSCDMBHeader & dmb) {
 header2006.l1Acc = dmb.l1a();
 header2006.cscID = dmb.dmbID();
 header2006.nTBins = 16;
 header2006.bxnCount = dmb.bxn();
 ///have to re-copy header into original buffer to update the data
 memcpy(theOriginalBuffer, &header2006, header2006.sizeForPacking());
}

unsigned short CSCALCTHeader::nLCTChipRead() const {///header2006 method
int count = 0;
 for(int i=0; i<7; ++i) {
   if( (header2006.lctChipRead>>i) & 1) ++count;
 }
 return count;
}


std::vector<CSCALCTDigi> CSCALCTHeader::ALCTDigis() const 
{ 
  std::vector<CSCALCTDigi> result;
  result.reserve(alcts.size());

  switch (firmwareVersion) {
  case 2006:
    {
      result = alcts2006.ALCTDigis();
      break;
    }
  case 2007:
    {
      for (unsigned int i=0; i<alcts.size(); ++i) {///loop over all alct words
	CSCALCTDigi digi(alcts[i].valid, alcts[i].quality, alcts[i].accel, alcts[i].pattern,
			 alcts[i].keyWire, (int)i/2, i%2+1);
	result.push_back(digi);
      }
      break;
    }
  default:
    edm::LogError("CSCALCTHeader|CSCRawToDigi")
      <<"Empty Digis: ALCT firmware version is bad/not defined!"; 
    break;
  }
  for(unsigned i = 0; i < result.size(); ++i) {result[i].setFullBX(BXNCount());}
  return result;

}


void CSCALCTHeader::add(const std::vector<CSCALCTDigi> & digis)
{
  if(firmwareVersion != 2006) {
    throw cms::Exception("CSCDigi2Raw")
      << "The ALCTDigis do not live in the ALCT header past the 2006 firmware version";
  }
  alcts2006.add(digis);
}


void CSCALCTHeader::selfTest()
{
  // tests packing and unpacking
  for(int station = 1; station <= 4; ++station)
  {
    CSCDetId detId(1, station, 1, 1, 0);
    CSCALCTDigi alct0(true, 1, 1, 1, 10, 6, 1);
    CSCALCTDigi alct1(true, 1, 1, 0, 11, 6, 2);

    std::vector<CSCALCTDigi> oldAlcts;
    oldAlcts.push_back(alct0);
    oldAlcts.push_back(alct1);
    CSCALCTHeader alctHeader(detId.iChamberType());

    alctHeader.add(oldAlcts);

    std::vector<CSCALCTDigi> alcts = alctHeader.ALCTDigis();
    assert(alcts[0] == alct0);
    assert(alcts[1] == alct1);
  }
}

std::ostream & operator<<(std::ostream & os, const CSCALCTHeader & header) 
{
  os << "ALCT HEADER CSCID " << header.CSCID()
     << "  L1ACC " << header.L1Acc() << std::endl;
  os << " time samples " << header.NTBins() << std::endl;
  return os;
}


