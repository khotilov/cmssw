#include <EventFilter/DTRawToDigi/plugins/DTDigiToRawModule.h>
#include <EventFilter/DTRawToDigi/plugins/DTDigiToRaw.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>
#include <DataFormats/DTDigi/interface/DTDigiCollection.h>
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>


#include <CondFormats/DataRecord/interface/DTReadOutMappingRcd.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <DataFormats/DTDigi/interface/DTLocalTriggerCollection.h>

#include <iostream>

using namespace edm;
using namespace std;

DTDigiToRawModule::DTDigiToRawModule(const edm::ParameterSet& ps) {
  produces<FEDRawDataCollection>();
  
  dduID = ps.getUntrackedParameter<int>("dduID", 770);
  debug = ps.getUntrackedParameter<bool>("debugMode", false);
  digicoll = ps.getUntrackedParameter<string>("digiColl", "dtunpacker");
  digibyType = ps.getUntrackedParameter<bool>("digibytype", true);
  
  useStandardFEDid_ = ps.getUntrackedParameter<bool>("useStandardFEDid", true);
  minFEDid_ = ps.getUntrackedParameter<int>("minFEDid", 731);
  maxFEDid_ = ps.getUntrackedParameter<int>("maxFEDid", 735);
  
  packer = new DTDigiToRaw(ps);
  if (debug) cout << "[DTDigiToRawModule]: constructor" << endl;
}

DTDigiToRawModule::~DTDigiToRawModule(){
  delete packer;
  if (debug) cout << "[DTDigiToRawModule]: destructor" << endl;
}


void DTDigiToRawModule::produce(Event & e, const EventSetup& iSetup) {

  auto_ptr<FEDRawDataCollection> fed_buffers(new FEDRawDataCollection);
  
  // Take digis from the event
  Handle<DTDigiCollection> digis;
  if (digibyType) {
    e.getByType(digis);
  }
  else {
    e.getByLabel(digicoll, digis);
  }
  
  // Load DTMap
  edm::ESHandle<DTReadOutMapping> map;
  iSetup.get<DTReadOutMappingRcd>().get( map );
  
  // Create the packed data
  int FEDIDmin = 0, FEDIDMax = 0;
  if (useStandardFEDid_){
    FEDIDmin = FEDNumbering::getDTFEDIds().first;
    FEDIDMax = FEDNumbering::getDTFEDIds().second;
  }
  else {
    FEDIDmin = minFEDid_;
    FEDIDMax = maxFEDid_;
  }
  
  for (int id=FEDIDmin; id<=FEDIDMax; ++id){
    
    packer->SetdduID(id);
    FEDRawData* rawData = packer->createFedBuffers(*digis, map);
    
    FEDRawData& fedRawData = fed_buffers->FEDData(id);
    fedRawData = *rawData;
  }
  // Put the raw data to the event
  e.put(fed_buffers);
  
}

