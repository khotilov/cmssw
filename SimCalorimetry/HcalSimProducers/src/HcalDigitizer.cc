#include "SimCalorimetry/HcalSimProducers/interface/HcalDigitizer.h"
#include "SimCalorimetry/HcalSimProducers/src/HcalTestHitGenerator.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalShape.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HFShape.h"
#include "SimCalorimetry/HcalSimAlgos/interface/ZDCShape.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalElectronicsSim.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalAmplifier.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalCoderFactory.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalHitCorrection.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSiPMShape.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSiPMHitResponse.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloTDigitizer.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloShapeIntegrator.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HPDNoiseGenerator.h"
#include <boost/foreach.hpp>
using namespace std;

namespace HcalDigitizerImpl {

  template<typename SIPMDIGITIZER>
  void fillSiPMCells(const vector<int> & siPMCells, SIPMDIGITIZER * siPMDigitizer)
  {
    std::vector<DetId> siPMDetIds;
    siPMDetIds.reserve(siPMCells.size());
    for(std::vector<int>::const_iterator idItr = siPMCells.begin();
        idItr != siPMCells.end(); ++idItr)
    {
      siPMDetIds.push_back(DetId(*idItr));
    }
    siPMDigitizer->setDetIds(siPMDetIds);
  }

  // if both exist, assume the SiPM one has cells filled, and
  // assign the rest to the HPD
  template<typename HPDDIGITIZER, typename SIPMDIGITIZER>
  void fillCells(const vector<DetId>& allCells,
                 HPDDIGITIZER * hpdDigitizer,
                 SIPMDIGITIZER * siPMDigitizer)
  {
    // if both digitizers exist, split up the cells
    if(siPMDigitizer && hpdDigitizer)
    {
      std::vector<DetId> siPMDetIds = siPMDigitizer->detIds();
      std::sort(siPMDetIds.begin(), siPMDetIds.end());
      std::vector<DetId> sortedCells = allCells;
      std::sort(sortedCells.begin(), sortedCells.end());
      std::vector<DetId> hpdCells;
      std::set_difference(sortedCells.begin(), sortedCells.end(),
                          siPMDetIds.begin(), siPMDetIds.end(),
                          std::back_inserter(hpdCells) );
      hpdDigitizer->setDetIds(hpdCells);
    }
    else
    {
      if(siPMDigitizer) siPMDigitizer->setDetIds(allCells);
      if(hpdDigitizer) hpdDigitizer->setDetIds(allCells);
    }
  }
} // namespace HcaiDigitizerImpl


HcalDigitizer::HcalDigitizer(const edm::ParameterSet& ps) 
: theParameterMap(new HcalSimParameterMap(ps)),
  theHcalShape(0),
  theSiPMShape(0),
  theHFShape(new HFShape()),
  theZDCShape(new ZDCShape()),
  theHcalIntegratedShape(0),
  theSiPMIntegratedShape(0),
  theHFIntegratedShape(new CaloShapeIntegrator(theHFShape)),
  theZDCIntegratedShape(new CaloShapeIntegrator(theZDCShape)),
  theHBHEResponse(0),
  theHBHESiPMResponse(0),
  theHOResponse(0),   
  theHOSiPMResponse(0),
  theHFResponse(new CaloHitResponse(theParameterMap, theHFIntegratedShape)),
  theZDCResponse(new CaloHitResponse(theParameterMap, theZDCIntegratedShape)),
  theHBHEAmplifier(0),
  theHFAmplifier(0),
  theHOAmplifier(0),
  theZDCAmplifier(0),
  theCoderFactory(0),
  theHBHEElectronicsSim(0),
  theHFElectronicsSim(0),
  theHOElectronicsSim(0),
  theZDCElectronicsSim(0),
  theHBHEHitFilter(),
  theHFHitFilter(ps.getParameter<bool>("doHFWindow")),
  theHOHitFilter(),
  theHOSiPMHitFilter(HcalOuter),
  theZDCHitFilter(),
  theHitCorrection(0),
  theNoiseGenerator(0),
  theNoiseHitGenerator(0),
  theHBHEDigitizer(0),
  theHBHESiPMDigitizer(0),
  theHODigitizer(0),
  theHOSiPMDigitizer(0),
  theHFDigitizer(0),
  theZDCDigitizer(0),
  isZDC(true),
  isHCAL(true),
  zdcgeo(true),
  hbhegeo(true),
  hogeo(true),
  hfgeo(true),
  theHOSiPMCode(ps.getParameter<edm::ParameterSet>("ho").getParameter<int>("siPMCode"))
{
  bool doNoise = ps.getParameter<bool>("doNoise");
  bool doEmpty = ps.getParameter<bool>("doEmpty");
  // need to make copies, because they might get different noise generators
  theHBHEAmplifier = new HcalAmplifier(theParameterMap, doNoise);
  theHFAmplifier = new HcalAmplifier(theParameterMap, doNoise);
  theHOAmplifier = new HcalAmplifier(theParameterMap, doNoise);
  theZDCAmplifier = new HcalAmplifier(theParameterMap, doNoise);
  theCoderFactory = new HcalCoderFactory(HcalCoderFactory::DB);
  theHBHEElectronicsSim = new HcalElectronicsSim(theHBHEAmplifier, theCoderFactory);
  theHFElectronicsSim = new HcalElectronicsSim(theHFAmplifier, theCoderFactory);
  theHOElectronicsSim = new HcalElectronicsSim(theHOAmplifier, theCoderFactory);
  theZDCElectronicsSim = new HcalElectronicsSim(theZDCAmplifier, theCoderFactory);

  // a code of 1 means make all cells SiPM
  std::vector<int> hbSiPMCells(ps.getParameter<edm::ParameterSet>("hb").getParameter<std::vector<int> >("siPMCells"));
  //std::vector<int> hoSiPMCells(ps.getParameter<edm::ParameterSet>("ho").getParameter<std::vector<int> >("siPMCells"));
  // 0 means none, 1 means all, and 2 means use hardcoded

  bool doHBHEHPD = hbSiPMCells.empty() || (hbSiPMCells[0] != 1);
  bool doHOHPD = (theHOSiPMCode != 1);
  bool doHBHESiPM = !hbSiPMCells.empty();
  bool doHOSiPM = (theHOSiPMCode != 0);

  if(doHBHEHPD || doHOHPD )
  {
    theHcalShape = new HcalShape();
    theHcalIntegratedShape = new CaloShapeIntegrator(theHcalShape);
  }
  if(doHBHESiPM || doHOSiPM )
  {
    theSiPMShape = new HcalSiPMShape();
    theSiPMIntegratedShape = new CaloShapeIntegrator(theSiPMShape);
  }

  if(doHBHEHPD)
  {
    theHBHEResponse = new CaloHitResponse(theParameterMap, theHcalIntegratedShape);
    theHBHEResponse->setHitFilter(&theHBHEHitFilter);
    theHBHEDigitizer = new HBHEDigitizer(theHBHEResponse, theHBHEElectronicsSim, doEmpty);
  }
  if(doHOHPD) 
  {
    theHOResponse = new CaloHitResponse(theParameterMap, theHcalIntegratedShape);
    theHOResponse->setHitFilter(&theHOHitFilter);
    theHODigitizer = new HODigitizer(theHOResponse, theHOElectronicsSim, doEmpty);
  }

  if(doHBHESiPM)
  {
    theHBHESiPMResponse = new HcalSiPMHitResponse(theParameterMap, theSiPMIntegratedShape);
    theHBHESiPMResponse->setHitFilter(&theHBHEHitFilter);
    theHBHESiPMDigitizer = new HBHEDigitizer(theHBHESiPMResponse, theHBHEElectronicsSim, doEmpty);
  }
  if(doHOSiPM)
  {
    theHOSiPMResponse = new HcalSiPMHitResponse(theParameterMap, theSiPMIntegratedShape);
    theHOSiPMResponse->setHitFilter(&theHOSiPMHitFilter);
    theHOSiPMDigitizer = new HODigitizer(theHOSiPMResponse, theHOElectronicsSim, doEmpty);
  }

  // if both are present, fill the SiPM cells now
  if(doHBHEHPD && doHBHESiPM)
  {
    HcalDigitizerImpl::fillSiPMCells(hbSiPMCells, theHBHESiPMDigitizer);
  }


  theHFResponse->setHitFilter(&theHFHitFilter);
  theZDCResponse->setHitFilter(&theZDCHitFilter);

  bool doTimeSlew = ps.getParameter<bool>("doTimeSlew");
  if(doTimeSlew) {
    // no time slewing for HF
    theHitCorrection = new HcalHitCorrection(theParameterMap);
    if(theHBHEResponse) theHBHEResponse->setHitCorrection(theHitCorrection);
    if(theHBHESiPMResponse) theHBHESiPMResponse->setHitCorrection(theHitCorrection);
    if(theHOResponse) theHOResponse->setHitCorrection(theHitCorrection);
    if(theHOSiPMResponse) theHOSiPMResponse->setHitCorrection(theHitCorrection);
    theZDCResponse->setHitCorrection(theHitCorrection);
  }

  theHFDigitizer = new HFDigitizer(theHFResponse, theHFElectronicsSim, doEmpty);
  theZDCDigitizer = new ZDCDigitizer(theZDCResponse, theZDCElectronicsSim, doEmpty);

  bool doHPDNoise = ps.getParameter<bool>("doHPDNoise");
  if(doHPDNoise) {
    //edm::ParameterSet hpdNoisePset = ps.getParameter<edm::ParameterSet>("HPDNoiseLibrary");
    theNoiseGenerator = new HPDNoiseGenerator(ps); 
    if(theHBHEDigitizer) theHBHEDigitizer->setNoiseSignalGenerator(theNoiseGenerator);
    if(theHBHESiPMDigitizer) theHBHESiPMDigitizer->setNoiseSignalGenerator(theNoiseGenerator);
  }

  if(ps.getParameter<bool>("injectTestHits") ){
    theNoiseHitGenerator = new HcalTestHitGenerator(ps);
    if(theHBHEDigitizer) theHBHEDigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
    if(theHBHESiPMDigitizer) theHBHESiPMDigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
    if(theHODigitizer) theHODigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
    if(theHOSiPMDigitizer) theHOSiPMDigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
    theHFDigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
    theZDCDigitizer->setNoiseHitGenerator(theNoiseHitGenerator);
  }

  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << "HcalDigitizer requires the RandomNumberGeneratorService\n"
         "which is not present in the configuration file.  You must add the service\n"
         "in the configuration file or remove the modules that require it.";
  }

  CLHEP::HepRandomEngine& engine = rng->getEngine();
  if(theHBHEDigitizer) theHBHEDigitizer->setRandomEngine(engine);
  if(theHBHESiPMDigitizer) theHBHESiPMDigitizer->setRandomEngine(engine);
  if(theHODigitizer) theHODigitizer->setRandomEngine(engine);
  if(theHOSiPMDigitizer) theHOSiPMDigitizer->setRandomEngine(engine);
  theHFDigitizer->setRandomEngine(engine);
  theZDCDigitizer->setRandomEngine(engine);

  if (theHitCorrection!=0) theHitCorrection->setRandomEngine(engine);

  hitsProducer_ = ps.getParameter<std::string>("hitsProducer");

}


HcalDigitizer::~HcalDigitizer() {
  delete theHBHEDigitizer;
  delete theHBHESiPMDigitizer;
  delete theHODigitizer;
  delete theHOSiPMDigitizer;
  delete theHFDigitizer;
  delete theZDCDigitizer;
  delete theParameterMap;
  delete theHcalShape;
  delete theSiPMShape;
  delete theHFShape;
  delete theZDCShape;
  delete theHcalIntegratedShape;
  delete theSiPMIntegratedShape;
  delete theHFIntegratedShape;
  delete theZDCIntegratedShape;
  delete theHBHEResponse;
  delete theHBHESiPMResponse;
  delete theHOResponse;
  delete theHOSiPMResponse;
  delete theHFResponse;
  delete theZDCResponse;
  delete theHBHEElectronicsSim;
  delete theHFElectronicsSim;
  delete theHOElectronicsSim;
  delete theZDCElectronicsSim;
  delete theHBHEAmplifier;
  delete theHFAmplifier;
  delete theHOAmplifier;
  delete theZDCAmplifier;
  delete theCoderFactory;
  delete theHitCorrection;
  delete theNoiseGenerator;
}


void HcalDigitizer::setHBHENoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator)
{
  noiseGenerator->setParameterMap(theParameterMap);
  noiseGenerator->setElectronicsSim(theHBHEElectronicsSim);
  theHBHEDigitizer->setNoiseSignalGenerator(noiseGenerator);
  theHBHEAmplifier->setNoiseSignalGenerator(noiseGenerator);
}

void HcalDigitizer::setHFNoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator)
{
  noiseGenerator->setParameterMap(theParameterMap);
  noiseGenerator->setElectronicsSim(theHFElectronicsSim);
  theHFDigitizer->setNoiseSignalGenerator(noiseGenerator);
  theHFAmplifier->setNoiseSignalGenerator(noiseGenerator);
}

void HcalDigitizer::setHONoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator)
{
  noiseGenerator->setParameterMap(theParameterMap);
  noiseGenerator->setElectronicsSim(theHOElectronicsSim);
  theHODigitizer->setNoiseSignalGenerator(noiseGenerator);
  theHOAmplifier->setNoiseSignalGenerator(noiseGenerator);
}

void HcalDigitizer::setZDCNoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator)
{
  noiseGenerator->setParameterMap(theParameterMap);
  noiseGenerator->setElectronicsSim(theZDCElectronicsSim);
  theZDCDigitizer->setNoiseSignalGenerator(noiseGenerator);
  theZDCAmplifier->setNoiseSignalGenerator(noiseGenerator);
}


void HcalDigitizer::produce(edm::Event& e, const edm::EventSetup& eventSetup) {
  // get the appropriate gains, noises, & widths for this event
  edm::ESHandle<HcalDbService> conditions;
  eventSetup.get<HcalDbRecord>().get(conditions);
  theHBHEAmplifier->setDbService(conditions.product());
  theHFAmplifier->setDbService(conditions.product());
  theHOAmplifier->setDbService(conditions.product());
  theZDCAmplifier->setDbService(conditions.product());

  theCoderFactory->setDbService(conditions.product());
  theParameterMap->setDbService(conditions.product());

  // get the correct geometry
  checkGeometry(eventSetup);
  
  // Step A: Get Inputs
  edm::Handle<CrossingFrame<PCaloHit> > cf, zdccf;

  // test access to SimHits for HcalHits and ZDC hits
  const std::string zdcHitsName(hitsProducer_+"ZDCHITS");
  e.getByLabel("mix", zdcHitsName , zdccf);
  MixCollection<PCaloHit> * colzdc = 0 ;
  if(zdccf.isValid()){
    colzdc = new MixCollection<PCaloHit>(zdccf.product());
  }else{
    edm::LogInfo("HcalDigitizer") << "We don't have ZDC hit collection available ";
    isZDC = false;
  }

  const std::string hcalHitsName(hitsProducer_+"HcalHits");
  e.getByLabel("mix", hcalHitsName ,cf);
  MixCollection<PCaloHit> * col = 0 ;
  if(cf.isValid()){
    col = new MixCollection<PCaloHit>(cf.product());
  }else{
    edm::LogInfo("HcalDigitizer") << "We don't have HCAL hit collection available ";
    isHCAL = false;
  }

  if(theHitCorrection != 0)
  {
    theHitCorrection->clear();
    if(isHCAL)
      theHitCorrection->fillChargeSums(*col);
    if(isZDC)
      theHitCorrection->fillChargeSums(*colzdc);
  }
  // Step B: Create empty output

  std::auto_ptr<HBHEDigiCollection> hbheResult(new HBHEDigiCollection());
  std::auto_ptr<HODigiCollection> hoResult(new HODigiCollection());
  std::auto_ptr<HFDigiCollection> hfResult(new HFDigiCollection());
  std::auto_ptr<ZDCDigiCollection> zdcResult(new ZDCDigiCollection());

  // Step C: Invoke the algorithm, passing in inputs and getting back outputs.
  if(isHCAL&&hbhegeo)
  {
    if(theHBHEDigitizer) theHBHEDigitizer->run(*col, *hbheResult);
    if(theHBHESiPMDigitizer) theHBHESiPMDigitizer->run(*col, *hbheResult);
  }
  if(isHCAL&&hogeo)
  {
    if(theHODigitizer) theHODigitizer->run(*col, *hoResult);
    if(theHOSiPMDigitizer) theHOSiPMDigitizer->run(*col, *hoResult);
  }
  if(isHCAL&&hfgeo)
    theHFDigitizer->run(*col, *hfResult);  
  if(isZDC&&zdcgeo) 
    theZDCDigitizer->run(*colzdc, *zdcResult);
  
  edm::LogInfo("HcalDigitizer") << "HCAL HBHE digis : " << hbheResult->size();
  edm::LogInfo("HcalDigitizer") << "HCAL HO digis   : " << hoResult->size();
  edm::LogInfo("HcalDigitizer") << "HCAL HF digis   : " << hfResult->size();
  edm::LogInfo("HcalDigitizer") << "HCAL ZDC digis   : " << zdcResult->size();

  // Step D: Put outputs into event
  e.put(hbheResult);
  e.put(hoResult);
  e.put(hfResult);
  e.put(zdcResult);
}


void HcalDigitizer::checkGeometry(const edm::EventSetup & eventSetup) {
  // TODO find a way to avoid doing this every event
  edm::ESHandle<CaloGeometry> geometry;
  eventSetup.get<CaloGeometryRecord>().get(geometry);
  if(theHBHEResponse) theHBHEResponse->setGeometry(&*geometry);
  if(theHBHESiPMResponse) theHBHESiPMResponse->setGeometry(&*geometry);
  if(theHOResponse) theHOResponse->setGeometry(&*geometry);
  if(theHOSiPMResponse) theHOSiPMResponse->setGeometry(&*geometry);
  theHFResponse->setGeometry(&*geometry);
  theZDCResponse->setGeometry(&*geometry);

  const vector<DetId>& hbCells =  geometry->getValidDetIds(DetId::Hcal, HcalBarrel);
  const vector<DetId>& heCells =  geometry->getValidDetIds(DetId::Hcal, HcalEndcap);
  const vector<DetId>& hoCells =  geometry->getValidDetIds(DetId::Hcal, HcalOuter);
  const vector<DetId>& hfCells =  geometry->getValidDetIds(DetId::Hcal, HcalForward);
  const vector<DetId>& zdcCells = geometry->getValidDetIds(DetId::Calo, HcalZDCDetId::SubdetectorId);
  //const vector<DetId>& hcalTrigCells = geometry->getValidDetIds(DetId::Hcal, HcalTriggerTower);
  //const vector<DetId>& hcalCalib = geometry->getValidDetIds(DetId::Calo, HcalCastorDetId::SubdetectorId);
  //std::cout<<"HcalDigitizer::CheckGeometry number of cells: "<<zdcCells.size()<<std::endl;
  if(zdcCells.empty()) zdcgeo = false;
  if(hbCells.empty() && heCells.empty()) hbhegeo = false;
  if(hoCells.empty()) hogeo = false;
  if(hfCells.empty()) hfgeo = false;
  // combine HB & HE



  vector<DetId> hbheCells = hbCells;
  hbheCells.insert(hbheCells.end(), heCells.begin(), heCells.end());

  HcalDigitizerImpl::fillCells(hbheCells, theHBHEDigitizer, theHBHESiPMDigitizer);
  //HcalDigitizerImpl::fillCells(hoCells, theHODigitizer, theHOSiPMDigitizer);
  buildHOSiPMCells(hoCells);
  theHFDigitizer->setDetIds(hfCells);
  theZDCDigitizer->setDetIds(zdcCells); 
}


void HcalDigitizer::buildHOSiPMCells(const std::vector<DetId>& allCells)
{
  // all HPD
  if(theHOSiPMCode == 0)
  {
    theHODigitizer->setDetIds(allCells);
  }
  else if(theHOSiPMCode == 1)
  {
    theHOSiPMDigitizer->setDetIds(allCells);
    // FIXME pick Zecotek or hamamatsu?
  } 
  else if(theHOSiPMCode == 2)// hardcode which are SiPM 
  {
    std::vector<HcalDetId> zecotekDetIds;
    std::vector<HcalDetId> hamamatsuDetIds;
    std::vector<DetId> siPMDetIds;
    std::vector<DetId> hpdDetIds;
    for(std::vector<DetId>::const_iterator detItr = allCells.begin();
        detItr != allCells.end(); ++detItr)
    {
      HcalDetId hcalId(*detItr);
      int ieta = hcalId.ieta();
      int iphi = hcalId.iphi(); 
      if ((ieta>=5 && ieta <= 10 )  && (iphi >=47 && iphi <=52))
      {
        zecotekDetIds.push_back(hcalId);
        siPMDetIds.push_back(*detItr);
      } 
      else if(((ieta>=5 && ieta <= 10 ) && (iphi >=53 && iphi <=58))
           || ((ieta>=11 && ieta <= 15 )  && (iphi >=59 && iphi <=70))){
        hamamatsuDetIds.push_back(hcalId);
        siPMDetIds.push_back(*detItr);
      }
      else {
        hpdDetIds.push_back(*detItr);
      }
    }
    assert(theHODigitizer);
    assert(theHOSiPMDigitizer);
    theHODigitizer->setDetIds(hpdDetIds);
    theHOSiPMDigitizer->setDetIds(siPMDetIds);
    theHOSiPMHitFilter.setDetIds(siPMDetIds);
    // FIXME not applying a HitFilter to the HPDs, for now
    theParameterMap->setHOZecotekDetIds(zecotekDetIds);
    theParameterMap->setHOHamamatsuDetIds(hamamatsuDetIds);

    // make sure we don't got through this exercise again
    theHOSiPMCode = -2;
  }
}

      
    

