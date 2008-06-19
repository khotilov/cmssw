#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CalibFormats/CastorObjects/interface/CastorDbService.h"
#include "CalibFormats/CastorObjects/interface/CastorDbRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>


#include "EventFilter/CastorRawToDigi/plugins/CastorDigiToRaw.h"

using namespace std;


CastorDigiToRaw::CastorDigiToRaw(edm::ParameterSet const& conf) :
  castorTag_(conf.getUntrackedParameter("CASTOR",edm::InputTag())),
  calibTag_(conf.getUntrackedParameter("CALIB",edm::InputTag())),
  trigTag_(conf.getUntrackedParameter("TRIG",edm::InputTag()))
{
  produces<FEDRawDataCollection>();
}

// Virtual destructor needed.
CastorDigiToRaw::~CastorDigiToRaw() { }  

// Functions that gets called by framework every event
void CastorDigiToRaw::produce(edm::Event& e, const edm::EventSetup& es)
{
  CastorPacker::Collections colls;

  
  // Step A: Get Inputs 
  edm::Handle<CastorDigiCollection> castor;
  if (!castorTag_.label().empty()) {
    e.getByLabel(castorTag_,castor);
    colls.castorCont=castor.product();
  }
//  edm::Handle<CastorCalibDigiCollection> Calib;
//  if (!calibTag_.label().empty()) {
//    e.getByLabel(calibTag_,Calib);
//    colls.calibCont=Calib.product();
//  }
//  edm::Handle<CastorTrigPrimDigiCollection> htp;
//  if (!trigTag_.label().empty()) {
//    e.getByLabel(trigTag_,htp);
//    colls.tpCont=htp.product();
//  }
  // get the mapping
  edm::ESHandle<CastorDbService> pSetup;
  es.get<CastorDbRecord>().get( pSetup );
  const CastorElectronicsMap* readoutMap=pSetup->getCastorMapping();
  // Step B: Create empty output
  std::auto_ptr<FEDRawDataCollection> raw=std::auto_ptr<FEDRawDataCollection>(new FEDRawDataCollection());

// change to this when getCastorFEDIds is added to FEDNumbering
//  const int ifed_first=FEDNumbering::getCastorFEDIds().first;
//  const int ifed_last=FEDNumbering::getCastorFEDIds().second;
  const int ifed_first=FEDNumbering::getHcalFEDIds().first;
  const int ifed_last=FEDNumbering::getHcalFEDIds().second;

  int orbitN=e.id().event();
  int bcnN=2000;

  // Step C: pack all requested FEDs
  for (int ifed=ifed_first; ifed<=ifed_last; ++ifed) {
    FEDRawData& fed = raw->FEDData(ifed);
    try {
      packer_.pack(ifed,ifed-ifed_first, e.id().event(),
		   orbitN, bcnN, colls, *readoutMap, fed);
    } catch (cms::Exception& e) {
      edm::LogWarning("Unpacking error") << e.what();
    } catch (...) {
      edm::LogWarning("Unpacking exception");
    }
  }


  e.put(raw);
}


