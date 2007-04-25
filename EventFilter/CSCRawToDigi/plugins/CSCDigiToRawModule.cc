/** \file
 *
 *  $Date: 2007/03/07 15:18:08 $
 *  $Revision: 1.5 $
 *  \author A. Tumanov - Rice
 */

#include <EventFilter/CSCRawToDigi/src/CSCDigiToRawModule.h>
#include <EventFilter/CSCRawToDigi/src/CSCDigiToRaw.h>
#include <DataFormats/FEDRawData/interface/FEDRawDataCollection.h>
#include <DataFormats/CSCDigi/interface/CSCStripDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
#include "DataFormats/Common/interface/Handle.h"
#include <FWCore/Framework/interface/Event.h>
#include "CondFormats/CSCObjects/interface/CSCReadoutMappingFromFile.h"

using namespace edm;
using namespace std;

CSCDigiToRawModule::CSCDigiToRawModule(const edm::ParameterSet & pset): 
  packer(new CSCDigiToRaw) {
  theMapping  = CSCReadoutMappingFromFile(pset);
  produces<FEDRawDataCollection>("CSCRawData"); 
}


CSCDigiToRawModule::~CSCDigiToRawModule(){
  delete packer;
}


void CSCDigiToRawModule::produce(Event & e, const EventSetup& c){


  auto_ptr<FEDRawDataCollection> fed_buffers(new FEDRawDataCollection);
  // Take digis from the event
  Handle<CSCStripDigiCollection> stripDigis;
  e.getByLabel("cscunpacker","MuonCSCStripDigi", stripDigis);
  Handle<CSCWireDigiCollection> wireDigis;
  e.getByLabel("cscunpacker","MuonCSCWireDigi", wireDigis);


  // Create the packed data
  packer->createFedBuffers(*stripDigis, *wireDigis, *(fed_buffers.get()), theMapping);
  
  // put the raw data to the event
  e.put(fed_buffers, "CSCRawData");
}


