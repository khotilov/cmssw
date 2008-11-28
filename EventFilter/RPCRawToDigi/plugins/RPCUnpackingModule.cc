#include "RPCUnpackingModule.h"
#include "EventFilter/RPCRawToDigi/interface/RPCRecordFormatter.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "EventFilter/RPCRawToDigi/interface/RPCRawDataCounts.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/RPCObjects/interface/RPCEMap.h"
#include "CondFormats/DataRecord/interface/RPCEMapRcd.h"
#include "RPCReadOutMappingWithFastSearch.h"
#include "EventFilter/RPCRawToDigi/interface/DataRecord.h"
#include "EventFilter/RPCRawToDigi/interface/ReadoutError.h"
#include "EventFilter/RPCRawToDigi/interface/EventRecords.h"
#include "EventFilter/RPCRawToDigi/interface/DebugDigisPrintout.h"
#include "EventFilter/RPCRawToDigi/interface/RPCRawSynchro.h"

#include <iostream>
#include <bitset>

using namespace edm;
using namespace std;
using namespace rpcrawtodigi;

typedef uint64_t Word64;


RPCUnpackingModule::RPCUnpackingModule(const edm::ParameterSet& pset) 
  : dataLabel_(pset.getUntrackedParameter<edm::InputTag>("InputLabel",edm::InputTag("source"))),
    doSynchro_(pset.getParameter<bool>("doSynchro")),
    eventCounter_(0),
    theCabling(0)
{
  produces<RPCDigiCollection>();
  produces<RPCRawDataCounts>();
  if (doSynchro_) produces<RPCRawSynchro::ProdItem>();
}

RPCUnpackingModule::~RPCUnpackingModule()
{ 
  delete theCabling;
}


void RPCUnpackingModule::produce(Event & ev, const EventSetup& es)
{
  static bool debug = edm::MessageDrop::instance()->debugEnabled;
  eventCounter_++; 
  if (debug) LogDebug ("RPCUnpacker::produce") <<"Beginning To Unpack Event: "<<eventCounter_;
 
  Handle<FEDRawDataCollection> allFEDRawData; 
  ev.getByLabel(dataLabel_,allFEDRawData); 

  static edm::ESWatcher<RPCEMapRcd> recordWatcher;
  static RPCReadOutMappingWithFastSearch readoutMappingSearch;

  if (recordWatcher.check(es)) {  
    delete theCabling; 
    ESHandle<RPCEMap> readoutMapping;
    if (debug) LogTrace("") << "record has CHANGED!!, initialise readout map!";
    es.get<RPCEMapRcd>().get(readoutMapping);
    theCabling = readoutMapping->convert();
    if (debug) LogTrace("") <<" READOUT MAP VERSION: " << theCabling->version() << endl;
    readoutMappingSearch.init(theCabling);
  }

  std::auto_ptr<RPCDigiCollection> producedRPCDigis(new RPCDigiCollection);
  std::auto_ptr<RPCRawDataCounts> producedRawDataCounts(new RPCRawDataCounts);
  std::auto_ptr<RPCRawSynchro::ProdItem> producedRawSynchoCounts;
  if (doSynchro_) producedRawSynchoCounts.reset(new RPCRawSynchro::ProdItem);

  std::pair<int,int> rpcFEDS=FEDNumbering::getRPCFEDIds();
 

  int status = 0;
  for (int fedId= rpcFEDS.first; fedId<=rpcFEDS.second; ++fedId){  

    const FEDRawData & rawData = allFEDRawData->FEDData(fedId);
    RPCRecordFormatter interpreter = 
        theCabling ? RPCRecordFormatter(fedId,&readoutMappingSearch) : RPCRecordFormatter(fedId,0);
    int triggerBX =0;
    int nWords = rawData.size()/sizeof(Word64);
    if (nWords==0) continue;

    //
    // check headers
    //
    const Word64* header = reinterpret_cast<const Word64* >(rawData.data()); header--;
    bool moreHeaders = true;
    while (moreHeaders) {
      header++;
      FEDHeader fedHeader( reinterpret_cast<const unsigned char*>(header));
      if (!fedHeader.check()) {
        producedRawDataCounts->addReadoutError(fedId, ReadoutError(ReadoutError::HeaderCheckFail)); 
        if (debug) LogTrace("") <<" ** PROBLEM **, header.check() failed, break"; 
        break; 
      }
      if ( fedHeader.sourceID() != fedId) {
        producedRawDataCounts->addReadoutError(fedId, ReadoutError(ReadoutError::InconsitentFedId)); 
        if (debug) LogTrace ("") <<" ** PROBLEM **, fedHeader.sourceID() != fedId"
            << "fedId = " << fedId<<" sourceID="<<fedHeader.sourceID(); 
      }
      triggerBX = fedHeader.bxID();
      moreHeaders = fedHeader.moreHeaders();
      if (debug) {
        ostringstream str;
        str <<"  header: "<< *reinterpret_cast<const bitset<64>*> (header) << endl;
        str <<"  header triggerType: " << fedHeader.triggerType()<<endl;
        str <<"  header lvl1ID:      " << fedHeader.lvl1ID() << endl;
        str <<"  header bxID:        " << fedHeader.bxID() << endl;
        str <<"  header sourceID:    " << fedHeader.sourceID() << endl;
        str <<"  header version:     " << fedHeader.version() << endl;
        LogTrace("") << str.str();
      }
    }

    //
    // check trailers
    //
    const Word64* trailer=reinterpret_cast<const Word64* >(rawData.data())+(nWords-1); trailer++;
    bool moreTrailers = true;
    while (moreTrailers) {
      trailer--;
      FEDTrailer fedTrailer(reinterpret_cast<const unsigned char*>(trailer));
      if ( !fedTrailer.check()) {
        producedRawDataCounts->addReadoutError(fedId, ReadoutError(ReadoutError::TrailerCheckFail));
        if (debug) LogTrace("") <<" ** PROBLEM **, trailer.check() failed, break";
        break;
      }
      if ( fedTrailer.lenght()!= nWords) {
        producedRawDataCounts->addReadoutError(fedId, ReadoutError(ReadoutError::InconsistentDataSize)); 
        if (debug) LogTrace("")<<" ** PROBLEM **, fedTrailer.lenght()!= nWords, break";
        break;
      }
      moreTrailers = fedTrailer.moreTrailers();
      if (debug) {
        ostringstream str;
        str <<" trailer: "<<  *reinterpret_cast<const bitset<64>*> (trailer) << endl; 
        str <<"  trailer lenght:    "<<fedTrailer.lenght()<<endl;
        str <<"  trailer crc:       "<<fedTrailer.crc()<<endl;
        str <<"  trailer evtStatus: "<<fedTrailer.evtStatus()<<endl;
        str <<"  trailer ttsBits:   "<<fedTrailer.ttsBits()<<endl;
        LogTrace("") << str.str();
      }
    }

    //
    // data records
    //
    if (debug) {
      ostringstream str;
      for (const Word64* word = header+1; word != trailer; word++) {
        str<<"    data: "<<*reinterpret_cast<const bitset<64>*>(word) << endl; 
      }
      LogTrace("") << str.str();
    }
    EventRecords event(triggerBX);
    for (const Word64* word = header+1; word != trailer; word++) {
      for( int iRecord=1; iRecord<=4; iRecord++){
        const DataRecord::Data* pRecord = reinterpret_cast<const DataRecord::Data* >(word+1)-iRecord;
        DataRecord record(*pRecord);
        event.add(record);
        if (debug) {
          std::ostringstream str;
          str <<"record: "<<record.print()<<" hex: "<<hex<<*pRecord<<dec;
          str <<" type:"<<record.type()<<DataRecord::print(record);
          if (event.complete()) {
            str<< " --> dccId: "<<fedId
               << " dccInputChannelNum: " <<event.recordSLD().rmb()
               << " tbLinkInputNum: "<<event.recordSLD().tbLinkInputNumber()
               << " lbNumInLink: "<<event.recordCD().lbInLink()
               << " partition "<<event.recordCD().partitionNumber()
               << " cdData "<<event.recordCD().partitionData();
          }
          LogTrace("") << str.str();
        }
        producedRawDataCounts->addDccRecord(fedId, record);
        int statusTMP = 0;
        if (event.complete()) statusTMP= 
            interpreter.recordUnpack( event, 
            producedRPCDigis.get(), producedRawDataCounts.get(), producedRawSynchoCounts.get()); 
        if (statusTMP != 0) status = statusTMP;
      }
    }
  }
  if (status) LogWarning(" RPCUnpackingModule - There was unpacking PROBLEM in this event");
  if (debug) LogTrace("") << DebugDigisPrintout()(producedRPCDigis.get()) << endl;
  ev.put(producedRPCDigis);  
  ev.put(producedRawDataCounts);
  if (doSynchro_) ev.put(producedRawSynchoCounts);
}
