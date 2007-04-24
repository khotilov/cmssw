/** \file
 *
 *  $Date: 2007/03/12 00:44:18 $
 *  $Revision: 1.22 $
 *  \author S. Argiro - N. Amapane - M. Zanetti 
 * FRC 060906
 */


#include <FWCore/Framework/interface/Event.h>
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <EventFilter/DTRawToDigi/plugins/DTUnpackingModule.h>
#include <DataFormats/FEDRawData/interface/FEDRawData.h>
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/DTDigi/interface/DTDigiCollection.h>
#include <DataFormats/DTDigi/interface/DTLocalTriggerCollection.h>

#include <CondFormats/DTObjects/interface/DTReadOutMapping.h>
#include <CondFormats/DataRecord/interface/DTReadOutMappingRcd.h>

#include <EventFilter/DTRawToDigi/interface/DTDDUWords.h>
#include <EventFilter/DTRawToDigi/plugins/DTDDUUnpacker.h>
#include <EventFilter/DTRawToDigi/plugins/DTROS25Unpacker.h>
#include <EventFilter/DTRawToDigi/plugins/DTROS8Unpacker.h>


using namespace edm;
using namespace std;

#include <iostream>


#define SLINK_WORD_SIZE 8 


DTUnpackingModule::DTUnpackingModule(const edm::ParameterSet& ps) :
  unpacker(0), numOfEvents(0)
{

  eventScanning = ps.getUntrackedParameter<int>("eventScanning",1000);

  const string &  dataType = ps.getParameter<string>("dataType");

  if (dataType == "DDU") {
    unpacker = new DTDDUUnpacker(ps);
  } else if (dataType == "ROS8") {
    unpacker = new DTROS8Unpacker(ps);
  } else if (dataType == "ROS25") {
    unpacker = new DTROS25Unpacker(ps);
  } 
  else {
    throw cms::Exception("InvalidParameter") << "DTUnpackingModule: dataType "
					     << dataType << " is unknown";
  }

  produces<DTDigiCollection>();
  produces<DTLocalTriggerCollection>();
}

DTUnpackingModule::~DTUnpackingModule(){
  cout<<"[DTUnpackingModule]: Destructor"<<endl;
  delete unpacker;
}


void DTUnpackingModule::produce(Event & e, const EventSetup& context){

  // Get the data from the event 
//   Handle<FEDRawDataCollection> rawdata;
//   e.getByLabel("DaqSource", rawdata);

  Handle<FEDRawDataCollection> rawdata;
  e.getByType(rawdata);

  // Get the mapping from the setup
  ESHandle<DTReadOutMapping> mapping;
  context.get<DTReadOutMappingRcd>().get(mapping);
  
  // Create the result i.e. the collections of MB Digis and SC local triggers
  auto_ptr<DTDigiCollection> product(new DTDigiCollection);
  auto_ptr<DTLocalTriggerCollection> product2(new DTLocalTriggerCollection);


  // Loop over the DT FEDs
  for (int id=FEDNumbering::getDTFEDIds().first; id<=FEDNumbering::getDTFEDIds().second; ++id){ 
    
    const FEDRawData& feddata = rawdata->FEDData(id);
    
    if (feddata.size()){
      
      // Unpack the DDU data

      unpacker->interpretRawData(reinterpret_cast<const unsigned int*>(feddata.data()), 
 				 feddata.size(), id, mapping, product, product2);
      
      numOfEvents++;      
      if (numOfEvents%eventScanning == 0) 
	cout<<"[DTUnpackingModule]: "<<numOfEvents<<" events analyzed"<<endl;
      
    }
  }

  // commit to the event  
  e.put(product);
  e.put(product2);
}

