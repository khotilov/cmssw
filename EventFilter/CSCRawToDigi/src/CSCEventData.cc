#include "EventFilter/CSCRawToDigi/interface/CSCEventData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCFEBData.h"
#include "DataFormats/CSCDigi/interface/CSCStripDigi.h"
#include "EventFilter/CSCRawToDigi/interface/CSCDMBTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCCFEBData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTHeader.h"
#include "EventFilter/CSCRawToDigi/interface/CSCAnodeData.h"
#include "EventFilter/CSCRawToDigi/interface/CSCALCTTrailer.h"
#include "EventFilter/CSCRawToDigi/interface/CSCTMBData.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iterator>


bool CSCEventData::debug = false;

CSCEventData::CSCEventData(int chamberType) : 
  theDMBHeader(), 
  theALCTHeader(0), 
  theAnodeData(0),
  theALCTTrailer(0),
  theTMBData(0),
  theDMBTrailer(),
  theChamberType(chamberType)
{
  for(unsigned i = 0; i < 5; ++i) {
    theCFEBData[i] = 0;
  }
}


CSCEventData::CSCEventData(unsigned short * buf)
 {
    // zero everything
    init();
    unsigned short * pos = buf;
    if(debug) {
      edm::LogInfo ("CSCEventData") << "The event data ";
      for(int i = 0; i < 16; ++i) {
	edm::LogInfo ("CSCEventData") << std::ios::hex << pos[i ] << " ";
      }
    }
   
    theDMBHeader = CSCDMBHeader(pos);
    if(!(theDMBHeader.check())) {
      edm::LogError ("CSCEventData")  << "Bad DMB Header??? " << " first four words: ";
      for(int i = 0; i < 4; ++i) {
	edm::LogError ("CSCEventData") << std::ios::hex << pos[i ] << " ";
      }
    }
  
    nalct_ = theDMBHeader.nalct();
    nclct_ = theDMBHeader.nclct();
   
      
    if (debug){
      edm::LogInfo ("CSCEventData") << "nalct = " << nalct();
      edm::LogInfo ("CSCEventData") << "nclct = " << nclct();
    }

    if (debug) {
     edm::LogInfo ("CSCEventData") << "size in words of DMBHeader" << theDMBHeader.sizeInWords();
     edm::LogInfo ("CSCEventData") << "sizeof(DMBHeader)" << sizeof(theDMBHeader); 
    }
   
    pos += theDMBHeader.sizeInWords();

    if (nalct() > 0) {
      theALCTHeader = new CSCALCTHeader( pos );
      pos += theALCTHeader->sizeInWords(); //size of the header

      if(!theALCTHeader->check()) {  
         edm::LogError ("CSCEventData") <<"+++WARNING: Corrupt ALCT data - won't attempt to decode";
      } else {
	//fill ALCT Digis
        theALCTHeader->ALCTDigis();


        theAnodeData = new CSCAnodeData(*theALCTHeader, pos);  
        pos += theAnodeData->sizeInWords(); // size of the data is determined during unpacking
        theALCTTrailer = new CSCALCTTrailer( pos );
        pos += theALCTTrailer->sizeInWords();
      }
    }

    if (nclct() > 0) {
      theTMBData = new CSCTMBData(pos);  //fill all TMB data
      pos += theTMBData->size();
    }

    for(int icfeb = 0; icfeb < 5; ++icfeb) {
      theCFEBData[icfeb] = 0;
      int cfeb_available = theDMBHeader.cfebAvailable(icfeb);
      // There seemed to be a bug in the DMB firmware used during 2003 test
      // beam, which resulted in occasional zero values of cfebAvailable
      // word in DMBHeader whereas in reality it was not zero.  So TMBHeader
      // might be a better place to get this info from.  But I am not quite
      // sure of it, so I leave it as is for now. -Sl.
      // int cfeb_available = 0;
      // if (nclct() > 0)
      //   cfeb_available = (theTMBData->tmbHeader().ActiveCFEBs() >> icfeb) & 1;
      if (cfeb_available) {
	// Fill CFEB data and convert it into cathode digis
	theCFEBData[icfeb] = new CSCCFEBData(icfeb, pos);
	theCFEBData[icfeb]->check();
	pos += theCFEBData[icfeb]->sizeInWords();
      }
    }
    
    if(debug) {
      edm::LogInfo ("CSCEventData") << "ready to unpack trailer ";
      for(unsigned i = 0; i < 8; ++i) {
	edm::LogInfo ("CSCEventData") << std::ios::hex << pos[i] <<" ";
      }
    }
    theDMBTrailer = *( (CSCDMBTrailer *) pos );
    pos += theDMBTrailer.sizeInWords();
    size_ = pos-buf;
 }


CSCEventData::CSCEventData(const CSCEventData & data) {
  copy(data);
}

CSCEventData::~CSCEventData() {
  destroy();
}


CSCEventData CSCEventData::operator=(const CSCEventData & data) {
  // check for self-assignment before destructing
  if(&data != this) destroy();
  copy(data);
  return *this;
}


void CSCEventData::init() {
  theALCTHeader = 0;
  theAnodeData = 0;
  theALCTTrailer = 0;
  theTMBData = 0;
  for(int icfeb = 0; icfeb < 5; ++icfeb) {
    theCFEBData[icfeb] = 0;
  }
}


void CSCEventData::copy(const CSCEventData & data) {
  init();
  theDMBHeader  = data.theDMBHeader;
  theDMBTrailer = data.theDMBTrailer;
  if(data.theALCTHeader != NULL)
    theALCTHeader  = new CSCALCTHeader(*(data.theALCTHeader));
  if(data.theAnodeData != NULL) 
    theAnodeData   = new CSCAnodeData(*(data.theAnodeData));
  if(data.theALCTTrailer != NULL) 
    theALCTTrailer = new CSCALCTTrailer(*(data.theALCTTrailer));
  if(data.theTMBData != NULL) 
    theTMBData     = new CSCTMBData(*(data.theTMBData));
  for(int icfeb = 0; icfeb < 5; ++icfeb) {
    theCFEBData[icfeb] = 0;
    if(data.theCFEBData[icfeb] != NULL) 
      theCFEBData[icfeb] = new CSCCFEBData(*(data.theCFEBData[icfeb]));
  }   
  nalct_ = data.nalct_;
  nclct_ = data.nclct_;
  size_  = data.size_;
  theChamberType = data.theChamberType;
}


void CSCEventData::destroy() {
  delete theALCTHeader;
  delete theAnodeData;
  delete theALCTTrailer;
  delete theTMBData;
  for(int icfeb = 0; icfeb < 5; ++icfeb) {
    delete theCFEBData[icfeb];
  }
}


std::vector<CSCStripDigi> CSCEventData::stripDigis(unsigned ilayer) const {
  assert(ilayer > 0 && ilayer <= 6);
  std::vector<CSCStripDigi> result;
  for(unsigned icfeb = 0; icfeb < 5; ++icfeb) {
    if(theCFEBData[icfeb] != NULL) {
      std::vector<CSCStripDigi> newDigis = theCFEBData[icfeb]->digis(ilayer);
      result.insert(result.end(), newDigis.begin(), newDigis.end());
    }
  }

  return result;
}

std::vector<CSCWireDigi> CSCEventData::wireDigis(unsigned ilayer) const {
  if(theAnodeData == 0) {
    return std::vector<CSCWireDigi>();
  } else {
    return theAnodeData->wireDigis(ilayer);
  }
}


std::vector < std::vector<CSCStripDigi> > CSCEventData::stripDigis() const {
  std::vector < std::vector<CSCStripDigi> > result;
  for (int layer = 1; layer <= 6; ++layer) {
    std::vector<CSCStripDigi> digis = stripDigis(layer);
    result.push_back(digis);
  }
  return result;
}

std::vector < std::vector<CSCWireDigi> > CSCEventData::wireDigis() const {
  std::vector < std::vector<CSCWireDigi> > result;
  for (int layer = 1; layer <= 6; ++layer) {
    result.push_back(wireDigis(layer));
  }
  return result;
}


CSCCFEBData* CSCEventData::cfebData(unsigned icfeb) const {
  return theCFEBData[icfeb];
}


CSCALCTHeader CSCEventData::alctHeader() const{
   if(nalct() == 0) throw("No ALCT for this chamber");
   return *theALCTHeader;
}

CSCALCTTrailer CSCEventData::alctTrailer() const{
   if(nalct() == 0) throw("No ALCT for this chamber");
   return *theALCTTrailer;
}


CSCAnodeData & CSCEventData::alctData() const {
   if(nalct() == 0) throw("No ALCT for this chamber");
   return *theAnodeData;
}

CSCTMBData & CSCEventData::tmbData() const {
  if(nclct() == 0) throw("No CLCT for this chamber");
   return *theTMBData;
}


CSCTMBHeader & CSCEventData::tmbHeader() const {
  if(nclct() == 0) throw("No CLCT for this chamber");
  return tmbData().tmbHeader();
}

CSCCLCTData & CSCEventData::clctData() const {
  if(nclct() == 0) throw("No CLCT for this chamber");
  return tmbData().clctData();
}


void CSCEventData::setEventInformation(int bxnum, int lvl1num) {
  theDMBHeader.setBXN(bxnum);
  theDMBHeader.setL1A(lvl1num);
  if(theALCTHeader) {
    theALCTHeader->setEventInformation(theDMBHeader);
  }
  if(theTMBData) {
    theTMBData->tmbHeader().setEventInformation(theDMBHeader);
  }
}
    

void CSCEventData::createALCTClasses() {
  assert(theChamberType>0);
  theALCTHeader = new CSCALCTHeader(theChamberType);
  theALCTHeader->setEventInformation(theDMBHeader);
  theAnodeData = new CSCAnodeData(*theALCTHeader);
  theALCTTrailer = new CSCALCTTrailer();
  // set data available flag
  theDMBHeader.addNALCT();

}

std::ostream & operator<<(std::ostream & os, const CSCEventData & evt) {
  for(int ilayer = 1; ilayer <= 6; ++ilayer) {

    std::vector<CSCStripDigi> stripDigis = evt.stripDigis(ilayer);
    //copy(stripDigis.begin(), stripDigis.end(), std::ostream_iterator<CSCStripDigi>(os, "\n"));
    //print your scas here
    std::vector<CSCWireDigi> wireDigis = evt.wireDigis(ilayer);
    //copy(wireDigis.begin(), wireDigis.end(), std::ostream_iterator<CSCWireDigi>(os, "\n"));
  }
  return os;
}

