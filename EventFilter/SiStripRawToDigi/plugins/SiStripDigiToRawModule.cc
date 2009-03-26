// Last commit: $Id: OldSiStripDigiToRawModule.cc,v 1.6 2008/07/04 14:07:47 bainbrid Exp $

#include "EventFilter/SiStripRawToDigi/plugins/SiStripDigiToRawModule.h"
#include "EventFilter/SiStripRawToDigi/interface/SiStripDigiToRaw.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "CondFormats/DataRecord/interface/SiStripFedCablingRcd.h"
#include <cstdlib>

using namespace std;

// -----------------------------------------------------------------------------
/** 
    Creates instance of DigiToRaw converter, defines EDProduct type.
*/
OldSiStripDigiToRawModule::OldSiStripDigiToRawModule( const edm::ParameterSet& pset ) :
  inputModuleLabel_( pset.getParameter<std::string>( "InputModuleLabel" ) ),
  inputDigiLabel_( pset.getParameter<std::string>( "InputDigiLabel" ) ),
  digiToRaw_(0),
  eventCounter_(0)
{
  if ( edm::isDebugEnabled() ) {
    LogDebug("DigiToRaw") 
      << "[OldSiStripDigiToRawModule::OldSiStripDigiToRawModule]"
      << " Constructing object...";
  }

  // Create instance of DigiToRaw formatter
  std::string mode = pset.getUntrackedParameter<std::string>("FedReadoutMode","VIRGIN_RAW");
  int16_t nbytes = pset.getUntrackedParameter<int>("AppendedBytes",0);
  bool use_fed_key = pset.getUntrackedParameter<bool>("UseFedKey",false);
  digiToRaw_ = new OldSiStripDigiToRaw( mode, nbytes, use_fed_key );
  
  produces<FEDRawDataCollection>();

}

// -----------------------------------------------------------------------------
/** */
OldSiStripDigiToRawModule::~OldSiStripDigiToRawModule() {
  if ( edm::isDebugEnabled() ) {
    LogDebug("DigiToRaw")
      << "[OldSiStripDigiToRawModule::~OldSiStripDigiToRawModule]"
      << " Destructing object...";
  }
  if ( digiToRaw_ ) delete digiToRaw_;
}

// -----------------------------------------------------------------------------
/** 
    Retrieves cabling map from EventSetup, retrieves a DetSetVector of
    SiStripDigis from Event, creates a FEDRawDataCollection
    (EDProduct), uses DigiToRaw converter to fill
    FEDRawDataCollection, attaches FEDRawDataCollection to Event.
*/
void OldSiStripDigiToRawModule::produce( edm::Event& iEvent, 
				      const edm::EventSetup& iSetup ) {

  eventCounter_++; 
  
  edm::ESHandle<SiStripFedCabling> cabling;
  iSetup.get<SiStripFedCablingRcd>().get( cabling );

  edm::Handle< edm::DetSetVector<SiStripDigi> > digis;
  iEvent.getByLabel( inputModuleLabel_, inputDigiLabel_, digis );

  std::auto_ptr<FEDRawDataCollection> buffers( new FEDRawDataCollection );
  
  digiToRaw_->createFedBuffers( iEvent, cabling, digis, buffers );

  iEvent.put( buffers );
  
}
