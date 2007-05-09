#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimMuon/CSCDigitizer/src/CSCDigiProducer.h"
#include "SimMuon/CSCDigitizer/src/CSCDigitizer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimMuon/CSCDigitizer/src/CSCConfigurableStripConditions.h"
#include "SimMuon/CSCDigitizer/src/CSCDbStripConditions.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"


CSCDigiProducer::CSCDigiProducer(const edm::ParameterSet& ps) 
:  theDigitizer(ps),
   theStripConditions(0)
{
  produces<CSCWireDigiCollection>("MuonCSCWireDigi");
  produces<CSCStripDigiCollection>("MuonCSCStripDigi");
  produces<CSCComparatorDigiCollection>("MuonCSCComparatorDigi");
  produces<DigiSimLinks>("MuonCSCWireDigiSimLinks");
  produces<DigiSimLinks>("MuonCSCStripDigiSimLinks");
  std::string stripConditions( ps.getParameter<std::string>("stripConditions") );
  if( stripConditions == "Configurable" )
  {
    edm::ParameterSet stripPSet = ps.getParameter<edm::ParameterSet>("strips");
    theStripConditions = new CSCConfigurableStripConditions(stripPSet);
  }
  else if ( stripConditions == "Database" )
  {
    theStripConditions = new CSCDbStripConditions();
  }
  else
  {
    throw cms::Exception("CSCDigiProducer") 
      << "Bad option for strip conditions: "
      << stripConditions;
  }
  theDigitizer.setStripConditions(theStripConditions);
}


void CSCDigiProducer::produce(edm::Event& e, const edm::EventSetup& eventSetup) {

  edm::Handle<CrossingFrame> cf;
  e.getByType(cf);

  // test access to SimHits
  const std::string hitsName("MuonCSCHits");

  std::auto_ptr<MixCollection<PSimHit> > 
    hits( new MixCollection<PSimHit>(cf.product(), hitsName) );


  // Create empty output

  std::auto_ptr<CSCWireDigiCollection> pWireDigis(new CSCWireDigiCollection());
  std::auto_ptr<CSCStripDigiCollection> pStripDigis(new CSCStripDigiCollection());
  std::auto_ptr<CSCComparatorDigiCollection> pComparatorDigis(new CSCComparatorDigiCollection());
  std::auto_ptr<DigiSimLinks> pWireDigiSimLinks(new DigiSimLinks() );
  std::auto_ptr<DigiSimLinks> pStripDigiSimLinks(new DigiSimLinks() );

  // find the geometry & conditions for this event
  edm::ESHandle<CSCGeometry> hGeom;
  eventSetup.get<MuonGeometryRecord>().get( hGeom );
  const CSCGeometry *pGeom = &*hGeom;

  theDigitizer.setGeometry( pGeom );


  // find the magnetic field
  edm::ESHandle<MagneticField> magfield;
  eventSetup.get<IdealMagneticFieldRecord>().get(magfield);

  theDigitizer.setMagneticField(&*magfield);


  // set the particle table
  edm::ESHandle < ParticleDataTable > pdt;
  eventSetup.getData( pdt );
  theDigitizer.setParticleDataTable(&*pdt);

  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
   throw cms::Exception("Configuration")
     << "CSCDigitizer requires the RandomNumberGeneratorService\n"
        "which is not present in the configuration file.  You must add the service\n"
        "in the configuration file or remove the modules that require it.";
  }

  CLHEP::HepRandomEngine& engine = rng->getEngine();

  theDigitizer.setRandomEngine(engine);
  theStripConditions->setRandomEngine(engine);

  theStripConditions->initializeEvent(eventSetup);

  // run the digitizer
  theDigitizer.doAction(*hits, *pWireDigis, *pStripDigis, *pComparatorDigis,
                        *pWireDigiSimLinks, *pStripDigiSimLinks);


  // store them in the event
  e.put(pWireDigis, "MuonCSCWireDigi");
  e.put(pStripDigis, "MuonCSCStripDigi");
  e.put(pComparatorDigis, "MuonCSCComparatorDigi");
  e.put(pWireDigiSimLinks, "MuonCSCWireDigiSimLinks");
  e.put(pStripDigiSimLinks, "MuonCSCStripDigiSimLinks");
}

