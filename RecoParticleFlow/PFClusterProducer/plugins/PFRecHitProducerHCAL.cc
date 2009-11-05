#include "RecoParticleFlow/PFClusterProducer/plugins/PFRecHitProducerHCAL.h"

#include <memory>

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/DetId/interface/DetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "Geometry/CaloTopology/interface/CaloTowerTopology.h"
#include "RecoCaloTools/Navigation/interface/CaloTowerNavigator.h"


using namespace std;
using namespace edm;

PFRecHitProducerHCAL::PFRecHitProducerHCAL(const edm::ParameterSet& iConfig)
  : PFRecHitProducer( iConfig ) 
{

 

  // access to the collections of rechits 

  
  inputTagHcalRecHitsHBHE_ =
    iConfig.getParameter<InputTag>("hcalRecHitsHBHE");
    
 
  inputTagCaloTowers_ = 
    iConfig.getParameter<InputTag>("caloTowers");
   
  thresh_HF_ = 
    iConfig.getParameter<double>("thresh_HF");
  navigation_HF_ = 
    iConfig.getParameter<bool>("navigation_HF");
  weight_HFem_ =
    iConfig.getParameter<double>("weight_HFem");
  weight_HFhad_ =
    iConfig.getParameter<double>("weight_HFhad");

  HCAL_Calib_ =
    iConfig.getParameter<bool>("HCAL_Calib");
  HF_Calib_ =
    iConfig.getParameter<bool>("HF_Calib");


  //--ab
  produces<reco::PFRecHitCollection>("HFHAD").setBranchAlias("HFHADRecHits");
  produces<reco::PFRecHitCollection>("HFEM").setBranchAlias("HFEMRecHits");
  //--ab
}



PFRecHitProducerHCAL::~PFRecHitProducerHCAL() {}



void PFRecHitProducerHCAL::createRecHits(vector<reco::PFRecHit>& rechits,
					 edm::Event& iEvent, 
					 const edm::EventSetup& iSetup ) {

  
  // this map is necessary to find the rechit neighbours efficiently
  //C but I should think about using Florian's hashed index to do this.
  //C in which case the map might not be necessary anymore
  //C however the hashed index does not seem to be implemented for HCAL
  // 
  // the key of this map is detId. 
  // the value is the index in the rechits vector
  map<unsigned,  unsigned > idSortedRecHits;
  map<unsigned,  unsigned > idSortedRecHitsHFEM;
  map<unsigned,  unsigned > idSortedRecHitsHFHAD;
  typedef map<unsigned, unsigned >::iterator IDH;  


  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
  // get the hcalBarrel geometry
  const CaloSubdetectorGeometry *hcalBarrelGeometry = 
    geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);

  // get the endcap geometry
  const CaloSubdetectorGeometry *hcalEndcapGeometry = 
    geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalEndcap);

  //--ab
  auto_ptr< vector<reco::PFRecHit> > HFHADRecHits( new vector<reco::PFRecHit> ); 
  auto_ptr< vector<reco::PFRecHit> > HFEMRecHits( new vector<reco::PFRecHit> ); 
  //--ab

  // 2 possibilities to make HCAL clustering :
  // - from the HCAL rechits
  // - from the CaloTowers. 
  // ultimately, clustering will be done taking CaloTowers as an 
  // input. This possibility is currently under investigation, and 
  // was thus made optional.

  // in the first step, we will fill the map of PFRecHits hcalrechits
  // either from CaloTowers or from HCAL rechits. 

  // in the second step, we will perform clustering on this map.

  if( !(inputTagCaloTowers_ == InputTag()) ) {
      
    edm::Handle<CaloTowerCollection> caloTowers; 
    CaloTowerTopology caloTowerTopology; 
    const CaloSubdetectorGeometry *caloTowerGeometry = 0; 
    // = geometry_->getSubdetectorGeometry(id)

    // get calotowers
    bool found = iEvent.getByLabel(inputTagCaloTowers_,
				   caloTowers);

    if(!found) {
      ostringstream err;
      err<<"could not find rechits "<<inputTagCaloTowers_;
      LogError("PFRecHitProducerHCAL")<<err.str()<<endl;
    
      throw cms::Exception( "MissingProduct", err.str());
    }
    else {
      assert( caloTowers.isValid() );
      
      // create rechits
      typedef CaloTowerCollection::const_iterator ICT;
      
      for(ICT ict=caloTowers->begin(); ict!=caloTowers->end();ict++) {
	  
	const CaloTower& ct = (*ict);
	  
	//C	
	if(!caloTowerGeometry) 
	  caloTowerGeometry = geoHandle->getSubdetectorGeometry(ct.id());

	  
	// get the hadronic energy.
	
	// Mike: Just ask for the Hadronic part only now!
	// Patrick : ARGH ! While this is ok for the HCAL, this is 
	// just wrong for the HF (in which em/had are artificially 
	// separated. 
	double energy = ct.hadEnergy();
	//Auguste: Photons in HF have no hadEnergy in fastsim: -> all RecHit collections are empty with photons.
	double energyEM = ct.emEnergy(); // For HF !
	//so test the total energy to deal with the photons in  HF:
	if( (energy+energyEM) < 1e-9 ) continue;
	  
	assert( ct.constituentsSize() );	  
	//Mike: The DetId will be taken by the first Hadronic constituent
	//         of the tower. That is only what we need

	
	//get the constituents of the tower
	const std::vector<DetId>& hits = ct.constituents();


	//Reserve the DetId we are looking for:

	HcalDetId detid;
	bool foundHCALConstituent = false;

 
	//Loop on them and search for HCAL
	for(unsigned int i=0;i< hits.size();++i)
	  {
	    if(hits[i].det()==DetId::Hcal)
	      {
		foundHCALConstituent = true;
		detid = hits[i];
		break;
	      }

	  }
	  
	reco::PFRecHit* pfrh = 0;
	//---ab: need 2 rechits for the HF:
	reco::PFRecHit* pfrhHFEM = 0;
	reco::PFRecHit* pfrhHFHAD = 0;

	if(foundHCALConstituent)
	  {
	    // std::cout << ", new Energy = " << energy << std::endl;
	    switch( detid.subdet() ) {
	    case HcalBarrel: 
	      {
		if(energy < thresh_Barrel_ ) continue;
		if ( HCAL_Calib_ ) energy   *= myPFCorr->getValues(detid)->getValue();
		pfrh = createHcalRecHit( detid, 
					 energy, 
					 PFLayer::HCAL_BARREL1, 
					 hcalBarrelGeometry,
					 ct.id().rawId() );
	      }
	      break;
	    case HcalEndcap:
	      {
		if(energy < thresh_Endcap_ ) continue;
		if ( HCAL_Calib_ ) energy   *= myPFCorr->getValues(detid)->getValue();
		pfrh = createHcalRecHit( detid, 
					 energy, 
					 PFLayer::HCAL_ENDCAP, 
					 hcalEndcapGeometry,
					 ct.id().rawId() );
	      }
	      break;
	    case HcalOuter:
	      {
	      }
	      break;
	    case HcalForward:
	      {
		//---ab: 2 rechits for HF:
		//double energyemHF = weight_HFem_*ct.emEnergy();
		//double energyhadHF = weight_HFhad_*ct.hadEnergy();
		double energyemHF = weight_HFem_ * energyEM;
		double energyhadHF = weight_HFhad_ * energy;
		// Some energy in the tower !
		if((energyemHF+energyhadHF) < thresh_HF_ ) continue;
		// The EM energy might be negative, as it amounts to Long - Short
		// In that case, put the EM "energy" in the HAD energy
		// Just to avoid systematic positive bias due to "Short" high fluctuations
		if ( energyemHF < thresh_HF_ ) { 
		  energyhadHF += energyemHF;
		  energyemHF = 0.;
		// Otherwise, create an EM rechit.
		} else {
		  if ( HF_Calib_ ) { 
		    energy   *= myPFCorr->getValues(detid)->getValue();
		    energyEM *= myPFCorr->getValues(detid)->getValue();
		  }
		  pfrhHFEM = createHcalRecHit( detid, 
					   energyemHF, 
					   PFLayer::HF_EM, 
					   hcalEndcapGeometry,
					   ct.id().rawId() );
		}
		// And create a HAD rechit with either the HAD or the HAD+EM energy. 
		if(energyhadHF > thresh_HF_ ){
		  pfrhHFHAD = createHcalRecHit( detid, 
					    energyhadHF, 
					    PFLayer::HF_HAD, 
					    hcalEndcapGeometry,
					    ct.id().rawId() );
		}

	      }
	      break;
	    default:
	      LogError("PFRecHitProducerHCAL")
		<<"CaloTower constituent: unknown layer : "
		<<detid.subdet()<<endl;
	    } 
	    if(pfrh) { 
	      rechits.push_back( *pfrh );
	      delete pfrh;
	      idSortedRecHits.insert( make_pair(ct.id().rawId(), 
						rechits.size()-1 ) ); 
	    }
	    //---ab: 2 rechits for HF:	   
	    if(pfrhHFEM) { 
	      HFEMRecHits->push_back( *pfrhHFEM );
	      delete pfrhHFEM;
	      idSortedRecHitsHFEM.insert( make_pair(ct.id().rawId(), 
						HFEMRecHits->size()-1 ) ); 
	    }
	    if(pfrhHFHAD) { 
	      HFHADRecHits->push_back( *pfrhHFHAD );
	      delete pfrhHFHAD;
	      idSortedRecHitsHFHAD.insert( make_pair(ct.id().rawId(), 
						HFHADRecHits->size()-1 ) ); 
	    }
	    //---ab	   
	  }
      }
      // do navigation 
      for(unsigned i=0; i<rechits.size(); i++ ) {
	findRecHitNeighboursCT( rechits[i], 
				idSortedRecHits, 
				caloTowerTopology);
      }
      for(unsigned i=0; i<HFEMRecHits->size(); i++ ) {
	findRecHitNeighboursCT( (*HFEMRecHits)[i], 
				idSortedRecHitsHFEM, 
				caloTowerTopology);
      }
      for(unsigned i=0; i<HFHADRecHits->size(); i++ ) {
	findRecHitNeighboursCT( (*HFHADRecHits)[i], 
				idSortedRecHitsHFHAD, 
				caloTowerTopology);
      }
      iEvent.put( HFHADRecHits,"HFHAD" );	
      iEvent.put( HFEMRecHits,"HFEM" );	
    }   
  }
  else if( !(inputTagHcalRecHitsHBHE_ == InputTag()) ) { 
    // clustering is not done on CaloTowers but on HCAL rechits.
       

    // get the hcal topology
    HcalTopology hcalTopology;
    
    // HCAL rechits 
    //    vector<edm::Handle<HBHERecHitCollection> > hcalHandles;  
    edm::Handle<HBHERecHitCollection>  hcalHandle;  

    
    bool found = iEvent.getByLabel(inputTagHcalRecHitsHBHE_, 
				   hcalHandle );

    if(!found) {
      ostringstream err;
      err<<"could not find rechits "<<inputTagHcalRecHitsHBHE_;
      LogError("PFRecHitProducerHCAL")<<err.str()<<endl;
    
      throw cms::Exception( "MissingProduct", err.str());
    }
    else {
      assert( hcalHandle.isValid() );
      
      const edm::Handle<HBHERecHitCollection>& handle = hcalHandle;
      for(unsigned irechit=0; irechit<handle->size(); irechit++) {
	const HBHERecHit& hit = (*handle)[irechit];
	
	double energy = hit.energy();
	
	reco::PFRecHit* pfrh = 0;
	

	const HcalDetId& detid = hit.detid();
	switch( detid.subdet() ) {
	case HcalBarrel:
	  {
	    if(energy < thresh_Barrel_ ) continue;
	    pfrh = createHcalRecHit( detid, 
				     energy, 
				     PFLayer::HCAL_BARREL1, 
				     hcalBarrelGeometry );
 	  }
	  break;
	case HcalEndcap:
	  {
	    if(energy < thresh_Endcap_ ) continue;
	    pfrh = createHcalRecHit( detid, 
				     energy, 
				     PFLayer::HCAL_ENDCAP, 
				     hcalEndcapGeometry );	  
 	  }
	  break;
	case HcalForward:
	  {
	    if(energy < thresh_HF_ ) continue;
	    pfrh = createHcalRecHit( detid, 
				     energy, 
				     PFLayer::HF_HAD, 
				     hcalEndcapGeometry );
 	  }
	  break;
	default:
	  LogError("PFRecHitProducerHCAL")
	    <<"HCAL rechit: unknown layer : "<<detid.subdet()<<endl;
	  continue;
	} 

	if(pfrh) { 
	  rechits.push_back( *pfrh );
	  delete pfrh;
	  idSortedRecHits.insert( make_pair(detid.rawId(), 
					    rechits.size()-1 ) ); 
	}
      }
      
      
      // do navigation:
      for(unsigned i=0; i<rechits.size(); i++ ) {
	
	findRecHitNeighbours( rechits[i], idSortedRecHits, 
			      hcalTopology, 
			      *hcalBarrelGeometry, 
			      hcalTopology,
			      *hcalEndcapGeometry);
      } // loop for navigation
    }  // endif hcal rechits were found
  } // endif clustering on rechits in hcal
}






reco::PFRecHit* 
PFRecHitProducerHCAL::createHcalRecHit( const DetId& detid,
					double energy,
					PFLayer::Layer layer,
					const CaloSubdetectorGeometry* geom,
					unsigned newDetId ) {
  
  const CaloCellGeometry *thisCell = geom->getGeometry(detid);
  if(!thisCell) {
    edm::LogError("PFRecHitProducerHCAL")
      <<"warning detid "<<detid.rawId()<<" not found in layer "
      <<layer<<endl;
    return 0;
  }
  
  const GlobalPoint& position = thisCell->getPosition();
  
  
  unsigned id = detid;
  if(newDetId) id = newDetId;
  reco::PFRecHit *rh = 
    new reco::PFRecHit( id,  layer, energy, 
			position.x(), position.y(), position.z(), 
			0,0,0 );
 
  
  
  
  // set the corners
  const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();

  assert( corners.size() == 8 );

  rh->setNECorner( corners[0].x(), corners[0].y(),  corners[0].z() );
  rh->setSECorner( corners[1].x(), corners[1].y(),  corners[1].z() );
  rh->setSWCorner( corners[2].x(), corners[2].y(),  corners[2].z() );
  rh->setNWCorner( corners[3].x(), corners[3].y(),  corners[3].z() );
 
  return rh;
}




void 
PFRecHitProducerHCAL::findRecHitNeighbours
( reco::PFRecHit& rh, 
  const map<unsigned,unsigned >& sortedHits, 
  const CaloSubdetectorTopology& barrelTopology, 
  const CaloSubdetectorGeometry& barrelGeometry, 
  const CaloSubdetectorTopology& endcapTopology, 
  const CaloSubdetectorGeometry& endcapGeometry ) {
  
  //cout<<"------PFRecHitProducerHcaL:findRecHitNeighbours navigation value "<<navigation_HF_<<endl;
 if(navigation_HF_ == false){
    if( rh.layer() == PFLayer::HF_HAD )
      return;
    if( rh.layer() == PFLayer::HF_EM )
      return;
  } 
  DetId detid( rh.detId() );

  const CaloSubdetectorTopology* topology = 0;
  const CaloSubdetectorGeometry* geometry = 0;
  const CaloSubdetectorGeometry* othergeometry = 0;
  
  switch( rh.layer() ) {
  case PFLayer::ECAL_ENDCAP: 
    topology = &endcapTopology;
    geometry = &endcapGeometry;
    break;
  case PFLayer::ECAL_BARREL: 
    topology = &barrelTopology;
    geometry = &barrelGeometry;
    break;
  case PFLayer::HCAL_ENDCAP:
    topology = &endcapTopology;
    geometry = &endcapGeometry;
    othergeometry = &barrelGeometry;
    break;
  case PFLayer::HCAL_BARREL1:
    topology = &barrelTopology;
    geometry = &barrelGeometry;
    othergeometry = &endcapGeometry;
    break;
  case PFLayer::PS1:
  case PFLayer::PS2:
    topology = &barrelTopology;
    geometry = &barrelGeometry;
    othergeometry = &endcapGeometry;
    break;
  default:
    assert(0);
  }
  
  assert( topology && geometry );

  CaloNavigator<DetId> navigator(detid, topology);

  DetId north = navigator.north();  
  
  DetId northeast(0);
  if( north != DetId(0) ) {
    northeast = navigator.east();  
  }
  navigator.home();


  DetId south = navigator.south();

  

  DetId southwest(0); 
  if( south != DetId(0) ) {
    southwest = navigator.west();
  }
  navigator.home();


  DetId east = navigator.east();
  DetId southeast;
  if( east != DetId(0) ) {
    southeast = navigator.south(); 
  }
  navigator.home();
  DetId west = navigator.west();
  DetId northwest;
  if( west != DetId(0) ) {   
    northwest = navigator.north();  
  }
  navigator.home();
    
  IDH i = sortedHits.find( north.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
  
  i = sortedHits.find( northeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
  
  i = sortedHits.find( south.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( east.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( west.rawId() );
  if(i != sortedHits.end() ) 
     rh.add4Neighbour( i->second );
   
  i = sortedHits.find( northwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    

}


void 
PFRecHitProducerHCAL::findRecHitNeighboursCT
( reco::PFRecHit& rh, 
  const map<unsigned, unsigned >& sortedHits, 
  const CaloSubdetectorTopology& topology ) {
  //cout<<"------PFRecHitProducerHcaL:findRecHitNeighboursCT navigation value "<<navigation_HF_<<endl;
  //  cout<<"----------- rechit print out"<<endl;
  // if(( rh.layer() == PFLayer::HF_HAD )||(rh.layer() == PFLayer::HF_EM)) {  
    
  //    cout<<rh<<endl;
    //  }
  if(navigation_HF_ == false){
    if( rh.layer() == PFLayer::HF_HAD )
      return;
    if( rh.layer() == PFLayer::HF_EM )
      return;
  }
  CaloTowerDetId ctDetId( rh.detId() );
    

  vector<DetId> northids = topology.north(ctDetId);
  vector<DetId> westids = topology.west(ctDetId);
  vector<DetId> southids = topology.south(ctDetId);
  vector<DetId> eastids = topology.east(ctDetId);


  CaloTowerDetId badId;

  // all the following detids will be CaloTowerDetId
  CaloTowerDetId north;
  CaloTowerDetId northwest;
  CaloTowerDetId west;
  CaloTowerDetId southwest;
  CaloTowerDetId south;
  CaloTowerDetId southeast;
  CaloTowerDetId east;
  CaloTowerDetId northeast;
  
  // for north and south, there is no ambiguity : 1 or 0 neighbours
  string err("PFRecHitProducerHCAL::findRecHitNeighboursCT : incorrect number of neighbours "); 
  char n[20];
  
  switch( northids.size() ) {
  case 0: 
    break;
  case 1: 
    north = northids[0];
    break;
  default:
    sprintf(n, "north: %d", northids.size() );
    err += n;
    throw( err ); 
  }

  switch( southids.size() ) {
  case 0: 
    break;
  case 1: 
    south = southids[0];
    break;
  default:
    sprintf(n, "south %d", southids.size() );
    err += n;
    throw( err ); 
  }
  
  // for east and west, one must take care 
  // of the pitch change in HCAL endcap.

  switch( eastids.size() ) {
  case 0: 
    break;
  case 1: 
    east = eastids[0];
    northeast = getNorth(east, topology);
    southeast = getSouth(east, topology);
    break;
  case 2:  
    // in this case, 0 is more on the north than 1
    east = eastids[0];
    northeast = getNorth(east, topology );
    southeast = getSouth(eastids[1], topology);    
    break;
  default:
    sprintf(n, "%d", eastids.size() );
    err += n;
    throw( err ); 
  }
  
  
  switch( westids.size() ) {
  case 0: 
    break;
  case 1: 
    west = westids[0];
    northwest = getNorth(west, topology);
    southwest = getSouth(west, topology);
    break;
  case 2:  
    // in this case, 0 is more on the north than 1
    west = westids[0];
    northwest = getNorth(west, topology );
    southwest = getSouth(westids[1], topology );    
    break;
  default:
    sprintf(n, "%d", westids.size() );
    err += n;
    throw( err ); 
  }




  // find and set neighbours
    
  IDH i = sortedHits.find( north.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
  
  i = sortedHits.find( northeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
  
  i = sortedHits.find( south.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( east.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( west.rawId() );
  if(i != sortedHits.end() ) 
     rh.add4Neighbour( i->second );
   
  i = sortedHits.find( northwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );

  //  cout<<"----------- rechit print out"<<endl;
  // if(( rh.layer() == PFLayer::HF_HAD )||(rh.layer() == PFLayer::HF_EM)) {  
    
  //   cout<<rh<<endl;
    //  }
}



DetId 
PFRecHitProducerHCAL::getSouth(const DetId& id, 
			       const CaloSubdetectorTopology& topology) {

  DetId south;
  vector<DetId> sids = topology.south(id);
  if(sids.size() == 1)
    south = sids[0];
  
  return south;
} 



DetId 
PFRecHitProducerHCAL::getNorth(const DetId& id, 
			       const CaloSubdetectorTopology& topology) {

  DetId north;
  vector<DetId> nids = topology.north(id);
  if(nids.size() == 1)
    north = nids[0];
  
  return north;
} 


