/// CSCDCCEventData.cc
/// 01/20/05 
/// A.Tumanov

#include "EventFilter/CSCRawToDigi/interface/CSCDCCEventData.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>
#include <cstdio>
#include "EventFilter/CSCRawToDigi/src/bitset_append.h"

bool CSCDCCEventData::debug = false;

CSCDCCEventData::CSCDCCEventData(int sourceId, int nDDUs, int bx, int l1a) 
: theDCCHeader(bx, l1a, sourceId) 
{

  theDDUData.reserve(nDDUs);
  CSCDDUHeader dduHeader(bx, l1a, sourceId);
  for(int i = 0; i < nDDUs; ++i) 
    {
      theDDUData.push_back(CSCDDUEventData(dduHeader));
    }
} 

CSCDCCEventData::CSCDCCEventData(unsigned short *buf, CSCDCCExaminer* examiner)
{
  unpack_data(buf, examiner);
}

CSCDCCEventData::~CSCDCCEventData() 
{
}


void CSCDCCEventData::unpack_data(unsigned short *buf, CSCDCCExaminer* examiner) 
{
  //for (int i=0;i<200;i++) {
  //  printf("%04x %04x %04x %04x\n",buf[i+3],buf[i+2],buf[i+1],buf[i]); 
  //  i+=3;
  //}
  theDDUData.clear();
  if (debug) 
    LogTrace ("CSCDCCEventData|CSCRawToDigi") << "CSCDCCEventData::unpack_data() is called";

  // decode DCC header (128 bits)
  if (debug) 
    LogTrace ("CSCDCCEventData|CSCRawToDigi") << "unpacking dcc header...";
  memcpy(&theDCCHeader, buf, theDCCHeader.sizeInWords()*2);
  //theDCCHeader = CSCDCCHeader(buf); // direct unpacking instead of bitfields
  buf += theDCCHeader.sizeInWords();

  //std::cout <<"Sandrik DCC Id = " << theDCCHeader.getCDFSourceId()  << std::endl;
 
  ///loop over DDUEventDatas
  while ( (buf[7]==0x8000)&&(buf[6]==0x0001)&&(buf[5]==0x8000))
    {
       CSCDDUEventData dduEventData(buf, examiner);
//	CSCDDUEventData dduEventData(buf);

      if (debug) LogTrace ("CSCDCCEventData|CSCRawToDigi") << " checking ddu data integrity ";
      if (dduEventData.check()) 
	{
	  theDDUData.push_back(dduEventData);
	  buf += dduEventData.sizeInWords();
	} 
      else
	{
	  if (debug) LogTrace("CSCDCCEventData|CSCRawToDigi") <<"DDU Data Check failed!  ";
	  break;
	}
      
    }
  
  if (debug)
    {
      LogTrace ("CSCDCCEventData|CSCRawToDigi") << "unpacking dcc trailer ";
      LogTrace ("CSCDCCEventData|CSCRawToDigi") << std::hex << buf[3] <<" "
				       << buf[2]<<" " << buf[1]<<" " << buf[0];
    }
	    
  //decode dcc trailer (128 bits)
  if (debug) LogTrace ("CSCDCCEventData|CSCRawToDigi") <<"decoding DCC trailer";
  memcpy(&theDCCTrailer, buf, theDCCTrailer.sizeInWords()*2);
  if (debug) LogTrace("CSCDCCEventData|CSCRawToDigi") << "checking DDU Trailer" << theDCCTrailer.check(); 
  buf += theDCCTrailer.sizeInWords();
  
}
	  

bool CSCDCCEventData::check() const 
{
  // the trailer counts in 64-bit words
  if (debug) 
    {
      LogTrace ("CSCDCCEventData|CSCRawToDigi") << "size in Words () = " << std::dec << sizeInWords();
    }

  return  theDCCHeader.check() && theDCCTrailer.check();
}


void CSCDCCEventData::addChamber(CSCEventData & chamber, int dduID, int dmbID)
{
  // first, find this DDU
  std::vector<CSCDDUEventData>::iterator dduItr;
  for(dduItr = theDDUData.begin(); dduItr != theDDUData.end(); ++dduItr)
  {
    if(dduItr->header().source_id() == dduID) continue;
  }
  if(dduItr == theDDUData.end())
  {
    // make a new one
    CSCDDUHeader newDDUHeader(dccHeader().getCDFBunchCounter(), 
                              dccHeader().getCDFEventNumber(), dduID);
    theDDUData.push_back(CSCDDUEventData(newDDUHeader));
    dduItr = theDDUData.end()-1;
    //FIXME should be ddu slot ID
    //dccHeader().setDAV(dduSlot);
  }
  dduItr->add(chamber, dmbID);
}
 

boost::dynamic_bitset<> CSCDCCEventData::pack() 
{
  boost::dynamic_bitset<> result( theDCCHeader.sizeInWords()*16);
  result = bitset_utilities::ushortToBitset(theDCCHeader.sizeInWords()*16, theDCCHeader.data());
  //std::cout <<"SANDRIK DCC size of header  in words"<< theDCCHeader.sizeInWords()*16 <<std::endl;  
  //std::cout <<"SANDRIK DCC size of header in bits"<< result.size()<<std::endl;
  //for(size_t i = 0; i < result.size(); ++i) {
  //  std::cout<<result[i];
  //  if (((i+1)%32)==0) std::cout<<std::endl;
  //}
  
  for(size_t i = 0; i < theDDUData.size(); ++i) 
    {
      result = bitset_utilities::append(result,theDDUData[i].pack());
      //std::cout <<"SANDRIK here is ddu data check ";
      //theDDUData[i].header().check();
      //std::cout <<std::endl;
      //bitset_utilities::printWords(result);
    }
  
  //std::cout <<"SANDRIK packed dcc size is "<<result.size()<<std::endl;
  //for(size_t i = 0; i < result.size(); ++i) {
  //  std::cout<<result[i];
  //  if (((i+1)%32)==0) std::cout<<std::endl;
  //}

  boost::dynamic_bitset<> dccTrailer = bitset_utilities::ushortToBitset(theDCCTrailer.sizeInWords()*16,
									theDCCTrailer.data());
  result = bitset_utilities::append(result,dccTrailer);
  //  bitset_utilities::printWords(result);
  return result;
}

