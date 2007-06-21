#include "Utilities/Timing/interface/TimingReport.h" 
#include "SimMuon/CSCDigitizer/src/CSCDigitizer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimMuon/CSCDigitizer/src/CSCDetectorHit.h"
#include "SimMuon/CSCDigitizer/src/CSCWireHitSim.h"
#include "SimMuon/CSCDigitizer/src/CSCStripHitSim.h"
#include "SimMuon/CSCDigitizer/src/CSCDriftSim.h"
#include "SimMuon/CSCDigitizer/src/CSCWireElectronicsSim.h"
#include "SimMuon/CSCDigitizer/src/CSCStripElectronicsSim.h"
#include "SimMuon/CSCDigitizer/src/CSCNeutronReader.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <iostream>


CSCDigitizer::CSCDigitizer(const edm::ParameterSet & p)
: theDriftSim(new CSCDriftSim()),
  theWireHitSim(new CSCWireHitSim(theDriftSim)),
  theStripHitSim(new CSCStripHitSim()),
  theWireElectronicsSim(new CSCWireElectronicsSim(p.getParameter<edm::ParameterSet>("wires"))),
  theStripElectronicsSim(new CSCStripElectronicsSim(p.getParameter<edm::ParameterSet>("strips"))),
  theNeutronReader(0),
  theCSCGeometry(0)
{
  if(p.getParameter<bool>("doNeutrons"))
  {
    theNeutronReader = new CSCNeutronReader(p.getParameter<edm::ParameterSet>("neutrons"));
  }
}


CSCDigitizer::~CSCDigitizer() {
  delete theNeutronReader;
  delete theStripElectronicsSim;
  delete theWireElectronicsSim;
  delete theStripHitSim;
  delete theWireHitSim;
  delete theDriftSim;
}



void CSCDigitizer::doAction(MixCollection<PSimHit> & simHits, 
                            CSCWireDigiCollection & wireDigis, 
                            CSCStripDigiCollection & stripDigis, 
                            CSCComparatorDigiCollection & comparators,
                            DigiSimLinks & wireDigiSimLinks,
                            DigiSimLinks & stripDigiSimLinks) 
{
  // arrange the hits by layer
  std::map<int, edm::PSimHitContainer> hitMap;
  for(MixCollection<PSimHit>::MixItr hitItr = simHits.begin();
      hitItr != simHits.end(); ++hitItr) 
  {
    hitMap[hitItr->detUnitId()].push_back(*hitItr);
  }

  // add neutron background, if needed
  if(theNeutronReader != 0)
  {
    theNeutronReader->addHits(hitMap);
  }

  // now loop over layers and run the simulation for each one
  for(std::map<int, edm::PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();
      hitMapItr != hitMap.end(); ++hitMapItr)
  {
    const CSCLayer * layer = findLayer(hitMapItr->first);
    const edm::PSimHitContainer & layerSimHits = hitMapItr->second;

    std::vector<CSCDetectorHit> newWireHits, newStripHits;
  
    LogTrace("CSCDigitizer") << "CSCDigitizer: found " << layerSimHits.size() <<" hit(s) in layer"
       << " E" << layer->id().endcap() << " S" << layer->id().station() << " R" << layer->id().ring()
       << " C" << layer->id().chamber() << " L" << layer->id().layer();

    // turn the edm::PSimHits into WireHits, using the WireHitSim
    {
      TimeMe t("CSCWireHitSim");
      newWireHits.swap(theWireHitSim->simulate(layer, layerSimHits));
    }
    if(!newWireHits.empty()) {
      TimeMe t("CSCStripHitSim");
      newStripHits.swap(theStripHitSim->simulate(layer, newWireHits));
    }

    // turn the hits into wire digis, using the electronicsSim
    {
      TimeMe t("CSCWireElectronicsSim");
      theWireElectronicsSim->simulate(layer, newWireHits);
      theWireElectronicsSim->fillDigis(wireDigis);
      wireDigiSimLinks.insert( theWireElectronicsSim->digiSimLinks() );
    }  
    {
      TimeMe t("CSCStripElectronicsSim");
      theStripElectronicsSim->simulate(layer, newStripHits);
      theStripElectronicsSim->fillDigis(stripDigis, comparators);
      stripDigiSimLinks.insert( theStripElectronicsSim->digiSimLinks() );
    }
  }
}


void CSCDigitizer::setMagneticField(const MagneticField * field) {
  theDriftSim->setMagneticField(field);
}


void CSCDigitizer::setStripConditions(CSCStripConditions * cond)
{
  theStripElectronicsSim->setStripConditions(cond);
}


void CSCDigitizer::setParticleDataTable(const ParticleDataTable * pdt)
{
  theWireHitSim->setParticleDataTable(pdt);
}


void CSCDigitizer::setRandomEngine(CLHEP::HepRandomEngine& engine)
{
  theWireHitSim->setRandomEngine(engine);
  theWireElectronicsSim->setRandomEngine(engine);
  theStripElectronicsSim->setRandomEngine(engine);
  if(theNeutronReader) theNeutronReader->setRandomEngine(engine);
}


const CSCLayer * CSCDigitizer::findLayer(int detId) const {
  assert(theCSCGeometry != 0);
  const GeomDetUnit* detUnit = theCSCGeometry->idToDetUnit(CSCDetId(detId));
  if(detUnit == 0)
  {
    throw cms::Exception("CSCDigiProducer") << "Invalid DetUnit: " << CSCDetId(detId)
      << "\nPerhaps your signal or pileup dataset are not compatible with the current release?";
  }  
  return dynamic_cast<const CSCLayer *>(detUnit);
}

