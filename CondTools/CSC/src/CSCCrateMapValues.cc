#include <memory>
#include "boost/shared_ptr.hpp"
#include <fstream>

#include "CondFormats/CSCObjects/interface/CSCCrateMap.h"
#include "CondFormats/DataRecord/interface/CSCCrateMapRcd.h"
#include "CondTools/CSC/interface/CSCCrateMapValues.h"
#include "CondTools/CSC/interface/CSCMap1.h"

CSCCrateMapValues::CSCCrateMapValues(const edm::ParameterSet& iConfig)
{
  //the following line is needed to tell the framework what
  // data is being produced
  mapObj = fillCrateMap();
  setWhatProduced(this,&CSCCrateMapValues::produceCrateMap);
  findingRecord<CSCCrateMapRcd>();
  //now do what ever other initialization is needed
}


CSCCrateMapValues::~CSCCrateMapValues()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete mapObj;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
CSCCrateMapValues::ReturnType
CSCCrateMapValues::produceCrateMap(const CSCCrateMapRcd& iRecord)
{
  //need a new object so to not be deleted at exit
  CSCCrateMap* mydata=new CSCCrateMap( *mapObj );
  return mydata;
  
}

 void CSCCrateMapValues::setIntervalFor(const edm::eventsetup::EventSetupRecordKey &, const edm::IOVSyncValue&,
 edm::ValidityInterval & oValidity)
 {
 oValidity = edm::ValidityInterval(edm::IOVSyncValue::beginOfTime(),edm::IOVSyncValue::endOfTime());
 
 }
