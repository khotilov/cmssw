// system include files
#include <memory>

#include "Validation/Performance/test/SimDigiDumper.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/DetSetVector.h"

// ecal calorimeter info
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// hcal calorimeter info
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

// silicon strip info
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"

// silicon pixel info
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

// muon geometry info
//#include "Geometry/DTGeometry/interface/DTGeometry.h"
//#include "Geometry/Records/interface/MuonGeometryRecord.h"
//#include "Geometry/Records/interface/MuonNumberingRecord.h"

// muon DT info
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "CondFormats/DTObjects/interface/DTT0.h"

// muon CSC Strip info
#include "DataFormats/CSCDigi/interface/CSCStripDigiCollection.h"

// muon CSC Wire info
#include "DataFormats/CSCDigi/interface/CSCWireDigiCollection.h"

// muon RPC info
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"

SimDigiDumper::SimDigiDumper( const edm::ParameterSet& iPSet )
{
  //get Labels to use to extract information
  ECalEBSrc_ = iPSet.getParameter<edm::InputTag>("ECalEBSrc");
  ECalEESrc_ = iPSet.getParameter<edm::InputTag>("ECalEESrc");
  ECalESSrc_ = iPSet.getParameter<edm::InputTag>("ECalESSrc");
  HCalDigi_ = iPSet.getParameter<edm::InputTag>("HCalDigi");
  SiStripSrc_ = iPSet.getParameter<edm::InputTag>("SiStripSrc"); 
  SiPxlSrc_ = iPSet.getParameter<edm::InputTag>("SiPxlSrc");
  MuDTSrc_ = iPSet.getParameter<edm::InputTag>("MuDTSrc");
  MuCSCStripSrc_ = iPSet.getParameter<edm::InputTag>("MuCSCStripSrc");
  MuCSCWireSrc_ = iPSet.getParameter<edm::InputTag>("MuCSCWireSrc");
  MuRPCSrc_ = iPSet.getParameter<edm::InputTag>("MuRPCSrc");

  // print out Parameter Set information being used
  std::cout
    << "\n===============================\n"
    << "Dumping event digis for the collections:\n"
    << "    ECalEBSrc     = " << ECalEBSrc_.label()
    << ":" << ECalEBSrc_.instance() << "\n"
    << "    ECalEESrc     = " << ECalEESrc_.label()
    << ":" << ECalEESrc_.instance() << "\n"
    << "    ECalESSrc     = " << ECalESSrc_.label()
    << ":" << ECalESSrc_.instance() << "\n"
    << "    HCalDigi       = " << HCalDigi_.label()
    << ":" << HCalDigi_.instance() << "\n"
    << "    SiStripSrc    = " << SiStripSrc_.label()
    << ":" << SiStripSrc_.instance() << "\n"
    << "    SiPixelSrc    = " << SiPxlSrc_.label()
    << ":" << SiPxlSrc_.instance() << "\n"
    << "    MuDTSrc       = " << MuDTSrc_.label()
    << ":" << MuDTSrc_.instance() << "\n"
    << "    MuCSCStripSrc = " << MuCSCStripSrc_.label()
    << ":" << MuCSCStripSrc_.instance() << "\n"
    << "    MuCSCWireSrc  = " << MuCSCWireSrc_.label()
    << ":" << MuCSCWireSrc_.instance() << "\n"
    << "    MuRPCSrc      = " << MuRPCSrc_.label()
    << ":" << MuRPCSrc_.instance() << "\n"
    << "===============================\n" 
    << std::endl;

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SimDigiDumper::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  using namespace edm;
  
  // ECAL Barrel
  
  bool isBarrel = true;
  edm::Handle<EBDigiCollection> EcalDigiEB;  
  const EBDigiCollection *EBdigis = 0;
  iEvent.getByLabel(ECalEBSrc_, EcalDigiEB);
  if (!EcalDigiEB.isValid()) {
    std::cout << "Unable to find EcalDigiEB in event!" << std::endl;
  }  
  else {
    EBdigis = EcalDigiEB.product();
    if ( EcalDigiEB->size() == 0) isBarrel = false;
    std::cout << "Ecal Barrel, digi multiplicity = " << EcalDigiEB->size() << std::endl;
    
    if (isBarrel) {
      // loop over digis
      for (unsigned int digis=0; digis<EcalDigiEB->size(); ++digis) {
        
        EBDataFrame ebdf = (*EBdigis)[digis];
        std::cout << ebdf << std::endl;
      }    
    }
  }

  // ECAL Endcap
  bool isEndcap = true;
  edm::Handle<EEDigiCollection> EcalDigiEE;  
  const EEDigiCollection *EEdigis = 0;
  iEvent.getByLabel(ECalEESrc_, EcalDigiEE);
  if (!EcalDigiEE.isValid()) {
    std::cout << "Unable to find EcalDigiEE in event!" << std::endl;
  }  
  else {
    EEdigis = EcalDigiEE.product();
    if ( EcalDigiEE->size() == 0) isEndcap = false;
    std::cout << "Ecal Endcap, digi multiplicity = " << EcalDigiEE->size() << std::endl;
    
    if (isEndcap) {
      // loop over digis
      for (unsigned int digis=0; digis<EcalDigiEE->size(); ++digis) {
        
        EEDataFrame eedf = (*EEdigis)[digis];
        std::cout << eedf << std::endl;
      }    
    }
  }
  
  // ECAL Preshower
  bool isPreshower = true;
  edm::Handle<ESDigiCollection> EcalDigiES;  
  const ESDigiCollection *ESdigis = 0;
  iEvent.getByLabel(ECalESSrc_, EcalDigiES);
  if (!EcalDigiES.isValid()) {
    std::cout << "Unable to find EcalDigiES in event!" << std::endl;
  }  
  else {
    ESdigis = EcalDigiES.product();
    if ( EcalDigiES->size() == 0) isPreshower = false;
    std::cout << "Ecal Preshower, digi multiplicity = " << EcalDigiES->size() << std::endl;
     
    if (isPreshower) {
      // loop over digis
      for (unsigned int digis=0; digis<EcalDigiES->size(); ++digis) {
         
        ESDataFrame esdf = (*ESdigis)[digis];
        std::cout << esdf << std::endl;
      }    
    }
  }
  
  // HBHE
  bool isHBHE = true;
  edm::Handle<HBHEDigiCollection> hbhe;
  const HBHEDigiCollection *HBHEdigis = 0;
  iEvent.getByLabel(HCalDigi_,hbhe);
  if (!hbhe.isValid()) {
    std::cout << "Unable to find HBHEDataFrame in event!" << std::endl;
  }
  else {
    HBHEdigis = hbhe.product();
    if ( hbhe->size() == 0 ) isHBHE = false;
    std::cout << "HBHE, digi multiplicity = " << hbhe->size() << std::endl;
    
    if (isHBHE) {
      //loop over digis
      for (unsigned int digis=0; digis<hbhe->size(); ++digis) {	
        HBHEDataFrame hehbdf = (*HBHEdigis)[digis];
        std::cout << hehbdf << std::endl;
	//edm::SortedCollection<HBHEDataFrame>::const_iterator ihbhe;
	//for  (ihbhe == hbhe->begin(); ihbhe != hbhe->end(); ihbhe++) {
	//std::cout << "Nothing" << std::endl;
	//std::cout << (*ihbhe) << std::endl;
      }
    }
  }

  // HO
  bool isHO = true;
  edm::Handle<HODigiCollection> ho;
  const HODigiCollection *HOdigis = 0;
  iEvent.getByLabel(HCalDigi_,ho);
  if (!ho.isValid()) {
    std::cout << "Unable to find HODataFrame in event!" << std::endl;
  }
  else {
    HOdigis = ho.product();
    if ( ho->size() == 0 ) isHO = false;
    std::cout << "HO, digi multiplicity = " << ho->size() << std::endl;
    
    if (isHO) {
      //loop over digis
      for (unsigned int digis=0; digis<ho->size(); ++digis) {	
        HODataFrame hodf = (*HOdigis)[digis];
        std::cout << hodf << std::endl;
      }
    }
  }

  // HF
  bool isHF = true;
  edm::Handle<HFDigiCollection> hf;
  const HFDigiCollection *HFdigis = 0;
  iEvent.getByLabel(HCalDigi_,hf);
  if (!hf.isValid()) {
    std::cout << "Unable to find HFDataFrame in event!" << std::endl;
  }
  else {
    HFdigis = hf.product();
    if ( hf->size() == 0 ) isHF = false;
    std::cout << "HF, digi multiplicity = " << hf->size() << std::endl;
    
    if (isHF) {
      //loop over digis
      for (unsigned int digis=0; digis<hf->size(); ++digis) {	
        HFDataFrame hodf = (*HFdigis)[digis];
        std::cout << hodf << std::endl;
      }
    }
  }



  // Strip Tracker
  bool isStrip = true;
  edm::Handle<edm::DetSetVector<SiStripDigi> > stripDigis;
  iEvent.getByLabel(SiStripSrc_, stripDigis);
  if (!stripDigis.isValid()) {
    std::cout << "Unable to find stripDigis in event!" << std::endl; 
  }
  else {
    if ( stripDigis->size() == 0 ) isStrip = false;
    std::cout << "Strip Tracker, digi multiplicity = " <<  stripDigis->size() << std::endl;
    if (isStrip) {
      edm::DetSetVector<SiStripDigi>::const_iterator DSViter;
      for (DSViter = stripDigis->begin(); DSViter != stripDigis->end();
           ++DSViter) {
        edm::DetSet<SiStripDigi>::const_iterator begin = DSViter->data.begin();
        edm::DetSet<SiStripDigi>::const_iterator end = DSViter->data.end();
        edm::DetSet<SiStripDigi>::const_iterator iter;
        unsigned int id = DSViter->id;
        DetId detId(id);
         
        if (detId.subdetId() == sdSiTIB) {
          std::cout << "TIB " << DSViter->data.size() << std::endl;
        }
        else if (detId.subdetId() == sdSiTOB) {
          std::cout << "TOB " << DSViter->data.size() << std::endl;
        }
        else if (detId.subdetId() == sdSiTID) {
          std::cout << "TID " << DSViter->data.size() << std::endl;
        }
        if (detId.subdetId() == sdSiTEC) {
          std::cout << "TEC " << DSViter->data.size() << std::endl;
        }
        for (iter = begin; iter != end; ++iter) {
          std::cout << (*iter) << std::endl;
        }
      }
    }
  }

  // Pixel Tracker
  bool isPixel = true;
  edm::Handle<edm::DetSetVector<PixelDigi> > pixelDigis;
  iEvent.getByLabel(SiPxlSrc_, pixelDigis);
  if (!pixelDigis.isValid()) {
    std::cout << "Unable to find pixelDigis in event!" << std::endl; 
  }
  else {
    if ( pixelDigis->size() == 0 ) isPixel = false;
    std::cout << "Pixel Tracker, digi multiplicity = " <<  pixelDigis->size() << std::endl;

    if (isPixel) {
      edm::DetSetVector<PixelDigi>::const_iterator DSViter;
      for (DSViter = pixelDigis->begin(); DSViter != pixelDigis->end();
           ++DSViter) {
        edm::DetSet<PixelDigi>::const_iterator begin = DSViter->data.begin();
        edm::DetSet<PixelDigi>::const_iterator end = DSViter->data.end();
        edm::DetSet<PixelDigi>::const_iterator iter;
        unsigned int id = DSViter->id;
        DetId detId(id);
         
        if (detId.subdetId() == sdPxlBrl) {
          std::cout << "Pixel barrel " << DSViter->data.size() << std::endl;
        }
        else if (detId.subdetId() == sdPxlFwd) {
          std::cout << "Pixel forward " << DSViter->data.size() << std::endl;
        }
        for (iter = begin; iter != end; ++iter) {
          std::cout << (*iter) << std::endl;
        }
         
      }
    }
  }

  // DT 
  bool isDT = true;
  edm::Handle<DTDigiCollection> dtDigis;
  //edm::Handle<DTLayerIdDTDigiMuonDigiCollection> dtDigis;
  iEvent.getByLabel(MuDTSrc_, dtDigis);
  if (!dtDigis.isValid()) {
    std::cout << "Unable to find dtDigis in event!" << std::endl;
  }
  unsigned int nDT = 0;
  if ( dtDigis->begin() == dtDigis->end() ) {
    isDT = false;
    std::cout << "dtDigis seem empty!" << std::endl;
  } 
  if (isDT) {
    DTDigiCollection::DigiRangeIterator dtLayerIt;
    for (dtLayerIt = dtDigis->begin(); 
	 dtLayerIt != dtDigis->end();
         ++dtLayerIt) {
      const DTDigiCollection::Range& range = (*dtLayerIt).second;
      std::cout << "DT layer = " << (*dtLayerIt).first << " digi " << std::endl;
      for (DTDigiCollection::const_iterator digiIt = range.first;
           digiIt != range.second; ++digiIt) {
        std::cout << (*digiIt) << std::endl;
        nDT++;
      }
    }
   }
  std::cout << "DT, digi multiplicity = " << nDT << std::endl;

  
  // CSC strip
  bool isCSCStrip = true;
  edm::Handle<CSCStripDigiCollection> cscStripDigis;
  iEvent.getByLabel(MuCSCStripSrc_, cscStripDigis);
  if (!cscStripDigis.isValid()) {
    std::cout << "Unable to find cscStripDigis in event!" << std::endl;
  }
  if ( cscStripDigis->begin() == cscStripDigis->end() ) isCSCStrip = false;
  unsigned int nCSCStrip = 0;
  
  if (isCSCStrip) {
    CSCStripDigiCollection::DigiRangeIterator detUnitIt;
    for (detUnitIt = cscStripDigis->begin(); detUnitIt != cscStripDigis->end();
         ++detUnitIt) {
      const CSCStripDigiCollection::Range& range = (*detUnitIt).second;
      std::cout << "CSC detid = " << (*detUnitIt).first << " digi " << std::endl;
      for (CSCStripDigiCollection::const_iterator digiIt = range.first;
           digiIt != range.second; ++digiIt) {
        std::cout << (*digiIt) << std::endl;
        nCSCStrip++;
      }
    }
  }
  std::cout << "CSC strip, digi multiplicity = " << nCSCStrip << std::endl;

  // CSC wire
  bool isCSCWire = true;
  edm::Handle<CSCWireDigiCollection> cscWireDigis;
  iEvent.getByLabel(MuCSCWireSrc_, cscWireDigis);
  if (!cscWireDigis.isValid()) {
    std::cout << "Unable to find cscWireDigis in event!" << std::endl;
  }
  if ( cscWireDigis->begin() == cscWireDigis->end() ) isCSCWire = false;
  unsigned int nCSCWire = 0;
  
  if (isCSCWire) {
    CSCWireDigiCollection::DigiRangeIterator detUnitIt;
    for (detUnitIt = cscWireDigis->begin(); detUnitIt != cscWireDigis->end();
         ++detUnitIt) {
      const CSCWireDigiCollection::Range& range = (*detUnitIt).second;
      std::cout << "CSC detid = " << (*detUnitIt).first << " digi " << std::endl;
      for (CSCWireDigiCollection::const_iterator digiIt = range.first;
           digiIt != range.second; ++digiIt) {
        std::cout << (*digiIt) << std::endl;
        nCSCWire++;
      }
    }
  }
  std::cout << "CSC wire, digi multiplicity = " << nCSCWire << std::endl;
  
  // RPC 
  bool isRPC = true;
  edm::Handle<RPCDigiCollection> rpcDigis;
  iEvent.getByLabel(MuRPCSrc_, rpcDigis);
  if (!rpcDigis.isValid()) {
    std::cout << "Unable to find rpcDigis in event!" << std::endl;
  }
  if ( rpcDigis->begin() == rpcDigis->end() ) isRPC = false;
  unsigned int nRPC = 0;
  
  if (isRPC) {
    RPCDigiCollection::DigiRangeIterator detUnitIt;
    for (detUnitIt = rpcDigis->begin(); detUnitIt != rpcDigis->end();
         ++detUnitIt) {
      const RPCDigiCollection::Range& range = (*detUnitIt).second;
      std::cout << "RPC detid = " << (*detUnitIt).first << " digi " << std::endl;
      for (RPCDigiCollection::const_iterator digiIt = range.first;
           digiIt != range.second; ++digiIt) {
        std::cout << (*digiIt) << std::endl;
        nRPC++;
      }
    }
  }
  std::cout << "RPC, digi multiplicity = " << nRPC << std::endl;
  
  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimDigiDumper);

