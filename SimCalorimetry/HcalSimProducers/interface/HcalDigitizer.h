#ifndef HcalSimProducers_HcalDigitizer_h
#define HcalSimProducers_HcalDigitizer_h

#include "SimCalorimetry/HcalSimAlgos/interface/HcalDigitizerTraits.h"
#include "SimCalorimetry/CaloSimAlgos/interface/CaloTDigitizer.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HBHEHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HFHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HOHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalHitFilter.h"
#include "SimCalorimetry/HcalSimAlgos/interface/ZDCHitFilter.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <vector>

class CaloHitResponse;
class HcalSimParameterMap;
class HcalAmplifier;
class HPDIonFeedbackSim;
class HcalCoderFactory;
class HcalElectronicsSim;
class HcalHitCorrection;
class HcalTimeSlewSim;
class HcalBaseSignalGenerator;
class HcalShapes;
class PCaloHit;
class PileUpEventPrincipal;

class HcalDigitizer
{
public:

  explicit HcalDigitizer(const edm::ParameterSet& ps);
  virtual ~HcalDigitizer();

  /**Produces the EDM products,*/
  void initializeEvent(edm::Event const& e, edm::EventSetup const& c);
  void accumulate(edm::Event const& e, edm::EventSetup const& c);
  void accumulate(PileUpEventPrincipal const& e, edm::EventSetup const& c);
  void finalizeEvent(edm::Event& e, edm::EventSetup const& c);
  void beginRun(const edm::EventSetup & es);
  void endRun();
  
  void setHBHENoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator);
  void setHFNoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator);
  void setHONoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator);
  void setZDCNoiseSignalGenerator(HcalBaseSignalGenerator * noiseGenerator);

private:
  void accumulateCaloHits(edm::Handle<std::vector<PCaloHit> > const& hcalHits, edm::Handle<std::vector<PCaloHit> > const& zdcHits, int bunchCrossing);

  /// some hits in each subdetector, just for testing purposes
  void fillFakeHits();
  /// make sure the digitizer has the correct list of all cells that
  /// exist in the geometry
  void checkGeometry(const edm::EventSetup& eventSetup);
  const CaloGeometry * theGeometry;
  void updateGeometry(const edm::EventSetup& eventSetup);

  void buildHOSiPMCells(const std::vector<DetId>& allCells, const edm::EventSetup& eventSetup);

  /** Reconstruction algorithm*/
  typedef CaloTDigitizer<HBHEDigitizerTraits> HBHEDigitizer;
  typedef CaloTDigitizer<HODigitizerTraits> HODigitizer;
  typedef CaloTDigitizer<HFDigitizerTraits> HFDigitizer;
  typedef CaloTDigitizer<ZDCDigitizerTraits> ZDCDigitizer;
 
  HcalSimParameterMap * theParameterMap;
  HcalShapes * theShapes;

  CaloHitResponse * theHBHEResponse;
  CaloHitResponse * theHBHESiPMResponse;
  CaloHitResponse * theHOResponse;
  CaloHitResponse * theHOSiPMResponse;
  CaloHitResponse * theHFResponse;
  CaloHitResponse * theZDCResponse;

  // we need separate amplifiers (and electronicssims)
  // because they might have separate noise generators
  HcalAmplifier * theHBHEAmplifier;
  HcalAmplifier * theHFAmplifier;
  HcalAmplifier * theHOAmplifier;
  HcalAmplifier * theZDCAmplifier;

  HPDIonFeedbackSim * theIonFeedback;
  HcalCoderFactory * theCoderFactory;

  HcalElectronicsSim * theHBHEElectronicsSim;
  HcalElectronicsSim * theHFElectronicsSim;
  HcalElectronicsSim * theHOElectronicsSim;
  HcalElectronicsSim * theZDCElectronicsSim;


  HBHEHitFilter theHBHEHitFilter;
  HFHitFilter   theHFHitFilter;
  HOHitFilter   theHOHitFilter;
  HcalHitFilter theHOSiPMHitFilter;
  ZDCHitFilter  theZDCHitFilter;

  HcalHitCorrection * theHitCorrection;
  HcalTimeSlewSim * theTimeSlewSim;
  CaloVNoiseSignalGenerator * theNoiseGenerator;
  CaloVNoiseHitGenerator * theNoiseHitGenerator;

  HBHEDigitizer * theHBHEDigitizer;
  HBHEDigitizer * theHBHESiPMDigitizer;
  HODigitizer* theHODigitizer;
  HODigitizer* theHOSiPMDigitizer;
  HFDigitizer* theHFDigitizer;
  ZDCDigitizer* theZDCDigitizer;

  // need to cache some DetIds for the digitizers,
  // if they don't come straight from the geometry
  std::vector<DetId> theHBHEDetIds;
  std::vector<DetId> theHOHPDDetIds;
  std::vector<DetId> theHOSiPMDetIds;

  bool isZDC,isHCAL,zdcgeo,hbhegeo,hogeo,hfgeo;

  std::string hitsProducer_;

  int theHOSiPMCode;
};

#endif


