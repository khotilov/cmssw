#include "RecoParticleFlow/PFClusterProducer/plugins/PFRecHitProducerHCAL.h"

#include <memory>

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

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

  shortFibre_Cut = iConfig.getParameter<double>("ShortFibre_Cut");
  longFibre_Fraction = iConfig.getParameter<double>("LongFibre_Fraction");

  ECAL_Compensate_ = iConfig.getParameter<bool>("ECAL_Compensate");
  ECAL_Threshold_ = iConfig.getParameter<double>("ECAL_Threshold");
  ECAL_Compensation_ = iConfig.getParameter<double>("ECAL_Compensation");
  ECAL_Dead_Code_ = iConfig.getParameter<unsigned int>("ECAL_Dead_Code");

  EM_Depth_ = iConfig.getParameter<double>("EM_Depth");
  HAD_Depth_ = iConfig.getParameter<double>("HAD_Depth");

  //--ab
  produces<reco::PFRecHitCollection>("HFHAD").setBranchAlias("HFHADRecHits");
  produces<reco::PFRecHitCollection>("HFEM").setBranchAlias("HFEMRecHits");
  //--ab
}



PFRecHitProducerHCAL::~PFRecHitProducerHCAL() {}



void PFRecHitProducerHCAL::createRecHits(vector<reco::PFRecHit>& rechits,
					 vector<reco::PFRecHit>& rechitsCleaned,
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
	const std::vector<DetId>& allConstituents = theTowerConstituentsMap->constituentsOf(ct.id());

	/*
	for(unsigned int i=0;i< hits.size();++i) {
	  if(hits[i].det()==DetId::Hcal) {
	    HcalDetId did = hits[i];
	    if ( did.subdet()==HcalEndcap || did.subdet()==HcalForward ) { 
	      //double en = hits[i].energy();
	      int ieta = did.ieta();
	      const CaloCellGeometry *thisCell = hcalEndcapGeometry->getGeometry(did);
	      const GlobalPoint& position = thisCell->getPosition();
	      if ( abs(ieta) > 27 && abs(ieta) < 33 && energy > 10. ) { 
		std::cout << "HE/HF hit " << i << " at eta = " << ieta 
			  << " with CT energy = " << energy 
			  << " at eta, z (hit) = " << position.eta() << " " << position.z()
			  << " at eta, z (cte) = " << ct.emPosition().eta() << " " << ct.emPosition().z()
			  << " at eta, z (cth) = " << ct.hadPosition().eta() << " " << ct.hadPosition().z()
			  << " at eta, z (cto) = " << ct.eta() << " " << ct.vz() 
			  << std::endl;
	      }
	    }
	  }
	}
	*/
	
	//Reserve the DetId we are looking for:

	HcalDetId detid;
	// EcalDetId edetid;
	bool foundHCALConstituent = false;

 
	//Loop on the calotower constituents and search for HCAL
	double dead = 0.;
	double alive = 0.;
	for(unsigned int i=0;i< hits.size();++i) {
	  if(hits[i].det()==DetId::Hcal) { 
	    foundHCALConstituent = true;
	    detid = hits[i];
	    // An HCAL tower was found: Look for dead ECAL channels in the same CaloTower.
	    if ( ECAL_Compensate_ && energy > ECAL_Threshold_ ) {
	      for(unsigned int j=0;j<allConstituents.size();++j) { 
		if ( allConstituents[j].det()==DetId::Ecal ) { 
		  alive += 1.;
		  EcalChannelStatus::const_iterator chIt = theEcalChStatus->find(allConstituents[j]);
		  unsigned int dbStatus = chIt != theEcalChStatus->end() ? chIt->getStatusCode() : 0;
		  if ( dbStatus > ECAL_Dead_Code_ ) dead += 1.;
		}
	      }
	    } 
	    // Protection: tower 29 in HF is merged with tower 30. 
	    // Just take the position of tower 30 in that case. 
	    if ( detid.subdet() == HcalForward && abs(detid.ieta()) == 29 ) continue; 
	    break;
	  }
	}

	// In case of dead ECAL channel, rescale the HCAL energy...
	double rescaleFactor = alive > 0. ? 1. + ECAL_Compensation_*dead/alive : 1.;
	  
	reco::PFRecHit* pfrh = 0;
	reco::PFRecHit* pfrhCleaned = 0;
	//---ab: need 2 rechits for the HF:
	reco::PFRecHit* pfrhHFEM = 0;
	reco::PFRecHit* pfrhHFHAD = 0;
	reco::PFRecHit* pfrhHFEMCleaned = 0;
	reco::PFRecHit* pfrhHFHADCleaned = 0;

	if(foundHCALConstituent)
	  {
	    // std::cout << ", new Energy = " << energy << std::endl;
	    switch( detid.subdet() ) {
	    case HcalBarrel: 
	      {
		if(energy < thresh_Barrel_ ) continue;
		if ( HCAL_Calib_ ) energy   *= myPFCorr->getValues(detid)->getValue();
		//if ( rescaleFactor > 1. ) 
		// std::cout << "Barrel HCAL energy rescaled from = " << energy << " to " << energy*rescaleFactor << std::endl;
		if ( rescaleFactor > 1. ) { 
		  pfrhCleaned = createHcalRecHit( detid, 
						  energy, 
						  PFLayer::HCAL_BARREL1, 
						  hcalBarrelGeometry,
						  ct.id().rawId() );
		  pfrhCleaned->setRescale(rescaleFactor);
		  energy *= rescaleFactor;
		}
		pfrh = createHcalRecHit( detid, 
					 energy, 
					 PFLayer::HCAL_BARREL1, 
					 hcalBarrelGeometry,
					 ct.id().rawId() );
		pfrh->setRescale(rescaleFactor);
	      }
	      break;
	    case HcalEndcap:
	      {
		if(energy < thresh_Endcap_ ) continue;
		if ( HCAL_Calib_ ) energy   *= myPFCorr->getValues(detid)->getValue();
		//if ( rescaleFactor > 1. ) 
		// std::cout << "End-cap HCAL energy rescaled from = " << energy << " to " << energy*rescaleFactor << std::endl;
		if ( rescaleFactor > 1. ) { 
		  pfrhCleaned = createHcalRecHit( detid, 
						  energy, 
						  PFLayer::HCAL_BARREL1, 
						  hcalBarrelGeometry,
						  ct.id().rawId() );
		  pfrhCleaned->setRescale(rescaleFactor);
		  energy *= rescaleFactor;
		}
		pfrh = createHcalRecHit( detid, 
					 energy, 
					 PFLayer::HCAL_ENDCAP, 
					 hcalEndcapGeometry,
					 ct.id().rawId() );
		pfrh->setRescale(rescaleFactor);
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

		// Some energy must be in the long fibres is there is some energy in the short fibres ! 
		double longFibre = energyemHF + energyhadHF/2.;
		double shortFibre = energyhadHF/2.;
		int ieta = detid.ieta();
		int iphi = detid.iphi();
		if ( shortFibre > shortFibre_Cut && longFibre/shortFibre < longFibre_Fraction ) {
		  // Check if the long-fibre hit was not cleaned already (because hot)
		  // In this case don't apply the cleaning
		  HcalDetId theLongDetId (HcalForward, ieta, iphi, 1);
		  const HcalChannelStatus* theStatus = theHcalChStatus->getValues(theLongDetId);
		  unsigned theStatusValue = theStatus->getValue();
		  // The channel is killed
		  if ( !theStatusValue ) { 
		    rescaleFactor = 0. ;
		    pfrhHFHADCleaned = createHcalRecHit( detid, 
							 shortFibre, 
							 PFLayer::HF_HAD, 
							 hcalEndcapGeometry,
							 ct.id().rawId() );
		    pfrhHFHADCleaned->setRescale(rescaleFactor);
		    /*
		    std::cout << "ieta/iphi = " << ieta << " " << iphi 
			      << ", Energy em/had/long/short = " 
			      << energyemHF << " " << energyhadHF << " "
			      << longFibre << " " << shortFibre << " " 
			      << ". The status value is " << theStatusValue
			      << ". Short fibres were cleaned." << std::endl;
		    */
		    shortFibre *= rescaleFactor;
		  }
		}

		// Determine EM and HAD after cleaning of short and long fibres
		energyhadHF = 2.*shortFibre;
		energyemHF = longFibre - shortFibre;

		// The EM energy might be negative, as it amounts to Long - Short
		// In that case, put the EM "energy" in the HAD energy
		// Just to avoid systematic positive bias due to "Short" high fluctuations
		if ( energyemHF < thresh_HF_ ) { 
		  energyhadHF += energyemHF;
		  energyemHF = 0.;
		}

		// Apply HCAL DPG calibration factors, if requested
		if ( HF_Calib_ ) { 
		  energyhadHF   *= myPFCorr->getValues(detid)->getValue();
		  energyemHF *= myPFCorr->getValues(detid)->getValue();
		}
		
		
		// Create an EM and a HAD rechit if above threshold.
		if ( energyemHF > thresh_HF_ || energyhadHF > thresh_HF_ ) { 
		  pfrhHFEM = createHcalRecHit( detid, 
					       energyemHF, 
					       PFLayer::HF_EM, 
					       hcalEndcapGeometry,
					       ct.id().rawId() );
		  pfrhHFHAD = createHcalRecHit( detid, 
						energyhadHF, 
						PFLayer::HF_HAD, 
						hcalEndcapGeometry,
						ct.id().rawId() );
		  pfrhHFEM->setEnergyUp(energyhadHF);
		  pfrhHFHAD->setEnergyUp(energyemHF);
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
	    if(pfrhCleaned) { 
	      rechitsCleaned.push_back( *pfrhCleaned );
	      delete pfrhCleaned;
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
	    if(pfrhHFEMCleaned) { 
	      rechitsCleaned.push_back( *pfrhHFEMCleaned );
	      delete pfrhHFEMCleaned;
	    }
	    if(pfrhHFHADCleaned) { 
	      rechitsCleaned.push_back( *pfrhHFHADCleaned );
	      delete pfrhHFHADCleaned;
	    }
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
  
  double depth_correction = 0.;
  switch ( layer ) { 
  case PFLayer::HF_EM:
    depth_correction = position.z() > 0. ? EM_Depth_ : -EM_Depth_;
    break;
  case PFLayer::HF_HAD:
    depth_correction = position.z() > 0. ? HAD_Depth_ : -HAD_Depth_;
    break;
  default:
    break;
  }

  unsigned id = detid;
  if(newDetId) id = newDetId;
  reco::PFRecHit *rh = 
    new reco::PFRecHit( id,  layer, energy, 
			position.x(), position.y(), position.z()+depth_correction, 
			0,0,0 );
 
  
  
  
  // set the corners
  const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();

  assert( corners.size() == 8 );

  rh->setNECorner( corners[0].x(), corners[0].y(),  corners[0].z()+depth_correction );
  rh->setSECorner( corners[1].x(), corners[1].y(),  corners[1].z()+depth_correction );
  rh->setSWCorner( corners[2].x(), corners[2].y(),  corners[2].z()+depth_correction );
  rh->setNWCorner( corners[3].x(), corners[3].y(),  corners[3].z()+depth_correction );
 
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
  CaloTowerDetId northwest2;
  CaloTowerDetId west;
  CaloTowerDetId west2;
  CaloTowerDetId southwest;
  CaloTowerDetId southwest2;
  CaloTowerDetId south;
  CaloTowerDetId southeast;
  CaloTowerDetId southeast2;
  CaloTowerDetId east;
  CaloTowerDetId east2;
  CaloTowerDetId northeast;
  CaloTowerDetId northeast2;
  
  // for north and south, there is no ambiguity : 1 or 0 neighbours
  stringstream err("PFRecHitProducerHCAL::findRecHitNeighboursCT : incorrect number of neighbours "); 
  
  switch( northids.size() ) {
  case 0: 
    break;
  case 1: 
    north = northids[0];
    break;
  default:
    err<<"north: "<<northids.size();
    throw( err.str() ); 
  }

  switch( southids.size() ) {
  case 0: 
    break;
  case 1: 
    south = southids[0];
    break;
  default:
    err<<"south: "<<southids.size();
    throw( err.str() ); 
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
    east2 = eastids[1];
    northeast = getNorth(east, topology );
    southeast = getSouth(east2, topology);    
    northeast2 = getNorth(northeast, topology );
    southeast2 = getSouth(southeast, topology);    
    break;
  default:
    err<<"eastids: "<<eastids.size();
    throw( err.str() ); 
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
    west2 = westids[1];
    northwest = getNorth(west, topology );
    southwest = getSouth(west2, topology );    
    northwest2 = getNorth(northwest, topology );
    southwest2 = getSouth(southwest, topology );    
    break;
  default:
    err<<"westids: "<< westids.size();
    throw( err.str() ); 
  }




  // find and set neighbours
    
  IDH i = sortedHits.find( north.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
  
  i = sortedHits.find( northeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
  
  i = sortedHits.find( northeast2.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
  
  i = sortedHits.find( south.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( southwest2.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( east.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( east2.rawId() );
  if(i != sortedHits.end() ) 
    rh.add4Neighbour( i->second );
    
  i = sortedHits.find( southeast.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( southeast2.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );
    
  i = sortedHits.find( west.rawId() );
  if(i != sortedHits.end() ) 
     rh.add4Neighbour( i->second );
   
  i = sortedHits.find( west2.rawId() );
  if(i != sortedHits.end() ) 
     rh.add4Neighbour( i->second );
   
  i = sortedHits.find( northwest.rawId() );
  if(i != sortedHits.end() ) 
    rh.add8Neighbour( i->second );

  i = sortedHits.find( northwest2.rawId() );
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


