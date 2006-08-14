#include "PluginManager/PluginManager.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimG4Core/Application81/interface/OscarProducer.h"
#include "SimG4Core/Application81/interface/G4SimEvent.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "SimG4Core/Watcher/interface/SimProducer.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"

#include <iostream>

OscarProducer::OscarProducer(edm::ParameterSet const & p) 
{    
    produces<edm::SimTrackContainer>();
    produces<edm::SimVertexContainer>();
    produces<edm::PSimHitContainer>("TrackerHitsPixelBarrelLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsPixelBarrelHighTof");
    produces<edm::PSimHitContainer>("TrackerHitsTIBLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsTIBHighTof");
    produces<edm::PSimHitContainer>("TrackerHitsTIDLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsTIDHighTof");
    produces<edm::PSimHitContainer>("TrackerHitsPixelEndcapLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsPixelEndcapHighTof");
    produces<edm::PSimHitContainer>("TrackerHitsTOBLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsTOBHighTof");
    produces<edm::PSimHitContainer>("TrackerHitsTECLowTof");
    produces<edm::PSimHitContainer>("TrackerHitsTECHighTof");
    
    produces<edm::PSimHitContainer>("TotemHitsT1");
    produces<edm::PSimHitContainer>("TotemHitsT2Gem");
    produces<edm::PSimHitContainer>("TotemHitsRP");
    
    produces<edm::PCaloHitContainer>("EcalHitsEB");
    produces<edm::PCaloHitContainer>("EcalHitsEE");
    produces<edm::PCaloHitContainer>("EcalHitsES");
    produces<edm::PCaloHitContainer>("HcalHits");
    produces<edm::PCaloHitContainer>("CaloHitsTk");
    produces<edm::PSimHitContainer>("MuonDTHits");
    produces<edm::PSimHitContainer>("MuonCSCHits");
    produces<edm::PSimHitContainer>("MuonRPCHits");
    produces<edm::PCaloHitContainer>("CastorPL");
    produces<edm::PCaloHitContainer>("CastorFI");
    produces<edm::PCaloHitContainer>("CastorBU");
    produces<edm::PCaloHitContainer>("CastorTU");
    
    produces<edm::PCaloHitContainer>("ZDCHITS"); 
    
    m_runManager = RunManager::init(p);

    //register any products 
    m_producers= m_runManager->producers();

    for(Producers::iterator itProd = m_producers.begin();
	itProd != m_producers.end();
	++itProd) {
       (*itProd)->registerProducts(*this);
    }
}

OscarProducer::~OscarProducer() 
{ 
  //this is causing a seg fault when an exception occurs while constructing
  // an HcalSD.  Need to check for memory problems. 
  if (m_runManager!=0) delete m_runManager; 
}

void OscarProducer::beginJob(const edm::EventSetup & es)
{
    std::cout << " OscarProducer initializing " << std::endl;
    m_runManager->initG4(es);
}
 
void OscarProducer::endJob()
{ std::cout << " OscarProducer terminating " << std::endl; }
 
void OscarProducer::produce(edm::Event & e, const edm::EventSetup & es)
{
    std::vector<SensitiveTkDetector*>& sTk = m_runManager->sensTkDetectors();
    std::vector<SensitiveCaloDetector*>& sCalo = m_runManager->sensCaloDetectors();

    try
    {
    m_runManager->produce(e,es);

    std::auto_ptr<edm::SimTrackContainer> p1(new edm::SimTrackContainer);
    std::auto_ptr<edm::SimVertexContainer> p2(new edm::SimVertexContainer);
    G4SimEvent * evt = m_runManager->simEvent();
    evt->load(*p1);
    evt->load(*p2);
    e.put(p1);
    e.put(p2);

    for (std::vector<SensitiveTkDetector*>::iterator it = sTk.begin(); it != sTk.end(); it++)
    {
	std::vector<std::string> v = (*it)->getNames();
	for (std::vector<std::string>::iterator in = v.begin(); in!= v.end(); in++)
	{
	    std::auto_ptr<edm::PSimHitContainer> product(new edm::PSimHitContainer);
 	    (*it)->fillHits(*product,*in);
	    e.put(product,*in);
	}
    }
    for (std::vector<SensitiveCaloDetector*>::iterator it = sCalo.begin(); it != sCalo.end(); it++)
    {
	std::vector<std::string>  v = (*it)->getNames();
	for (std::vector<std::string>::iterator in = v.begin(); in!= v.end(); in++)
	{
	    std::auto_ptr<edm::PCaloHitContainer> product(new edm::PCaloHitContainer);
	    (*it)->fillHits(*product,*in);
	    e.put(product,*in);
	}
    }

    for(Producers::iterator itProd = m_producers.begin();
	itProd != m_producers.end();
	++itProd) {
       (*itProd)->produce(e,es);
    }
    }
    catch ( const SimG4Exception& simg4ex )
    {
       std::cout << " SimG4Exception caght !" << std::endl ;
       std::cout << " " << simg4ex.what() << std::endl ;
       m_runManager->abortEvent() ;
       throw edm::Exception( edm::errors::EventCorruption ) ;
    }
}
 
DEFINE_FWK_MODULE(OscarProducer)
 
