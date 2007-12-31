#include "Validation/EcalClusters/interface/EgammaSuperClusters.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"

#include "DataFormats/GeometryVector/interface/Pi.h"

EgammaSuperClusters::EgammaSuperClusters( const edm::ParameterSet& ps )
{
	outputFile_ = ps.getUntrackedParameter<std::string>("outputFile", "");
	CMSSW_Version_ = ps.getUntrackedParameter<std::string>("CMSSW_Version", "");

	verboseDBE_ = ps.getUntrackedParameter<bool>("verboseDBE", false);

	hist_min_Size_ = ps.getParameter<double>("hist_min_Size");
	hist_max_Size_ = ps.getParameter<double>("hist_max_Size");
	hist_bins_Size_ = ps.getParameter<int>   ("hist_bins_Size");

	hist_min_NumBC_ = ps.getParameter<double>("hist_min_NumBC");
	hist_max_NumBC_ = ps.getParameter<double>("hist_max_NumBC");
	hist_bins_NumBC_ = ps.getParameter<int>   ("hist_bins_NumBC");

	hist_min_ET_ = ps.getParameter<double>("hist_min_ET");
	hist_max_ET_ = ps.getParameter<double>("hist_max_ET");
	hist_bins_ET_ = ps.getParameter<int>   ("hist_bins_ET");

	hist_min_Eta_ = ps.getParameter<double>("hist_min_Eta");
	hist_max_Eta_ = ps.getParameter<double>("hist_max_Eta");
	hist_bins_Eta_ = ps.getParameter<int>   ("hist_bins_Eta");

	hist_min_Phi_ = ps.getParameter<double>("hist_min_Phi");
	hist_max_Phi_ = ps.getParameter<double>("hist_max_Phi");
	hist_bins_Phi_ = ps.getParameter<int>   ("hist_bins_Phi");

	hist_min_S1toS9_ = ps.getParameter<double>("hist_min_S1toS9");
	hist_max_S1toS9_ = ps.getParameter<double>("hist_max_S1toS9");
	hist_bins_S1toS9_ = ps.getParameter<int>   ("hist_bins_S1toS9");

	hist_min_S25toE_ = ps.getParameter<double>("hist_min_S25toE");
	hist_max_S25toE_ = ps.getParameter<double>("hist_max_S25toE");
	hist_bins_S25toE_ = ps.getParameter<int>   ("hist_bins_S25toE");

	hist_min_EToverTruth_ = ps.getParameter<double>("hist_min_EToverTruth");
	hist_max_EToverTruth_ = ps.getParameter<double>("hist_max_EToverTruth");
	hist_bins_EToverTruth_ = ps.getParameter<int>   ("hist_bins_EToverTruth");
	
	hist_min_deltaEta_ = ps.getParameter<double>("hist_min_deltaEta");
	hist_max_deltaEta_ = ps.getParameter<double>("hist_max_deltaEta");
	hist_bins_deltaEta_ = ps.getParameter<int>   ("hist_bins_deltaEta");

	MCTruthCollection_ = ps.getParameter<edm::InputTag>("MCTruthCollection");
	hybridBarrelSuperClusterCollection_ = ps.getParameter<edm::InputTag>("hybridBarrelSuperClusterCollection");
  	islandBarrelSuperClusterCollection_ = ps.getParameter<edm::InputTag>("islandBarrelSuperClusterCollection");
  	islandEndcapSuperClusterCollection_ = ps.getParameter<edm::InputTag>("islandEndcapSuperClusterCollection");
	hybridBarrelClusterShapeAssociation_ = ps.getParameter<edm::InputTag>("hybridBarrelClusterShapeAssociation");
	islandBarrelClusterShapeAssociation_ = ps.getParameter<edm::InputTag>("islandBarrelClusterShapeAssociation");
	islandEndcapClusterShapeAssociation_ = ps.getParameter<edm::InputTag>("islandEndcapClusterShapeAssociation");
}

EgammaSuperClusters::~EgammaSuperClusters() {}

void EgammaSuperClusters::beginJob(edm::EventSetup const&) 
{
  	dbe_ = edm::Service<DaqMonitorBEInterface>().operator->();                   

  	if ( verboseDBE_ )
	{
  	dbe_->setVerbose(1);
		dbe_->showDirStructure();
	}
	else 
		dbe_->setVerbose(0);

	dbe_->setCurrentFolder("CMSSW_"+CMSSW_Version_+"/EcalClusters/SuperClusters/");

	hist_HybridEB_SC_Size_ 
		= dbe_->book1D("hist_HybridEB_SC_Size_","# Super Clusters from Hybrid in Barrel",
			hist_bins_Size_,hist_min_Size_,hist_max_Size_);
  	hist_IslandEB_SC_Size_ 
		= dbe_->book1D("hist_IslandEB_SC_Size_","# Super Clusters from Island in Barrel",
			hist_bins_Size_,hist_min_Size_,hist_max_Size_);
  	hist_IslandEE_SC_Size_ 
		= dbe_->book1D("hist_IslandEE_SC_Size_","# Super Clusters from Island in Endcap",
			hist_bins_Size_,hist_min_Size_,hist_max_Size_);

	hist_HybridEB_SC_NumBC_ 
		= dbe_->book1D("hist_HybridEB_SC_NumBC_","# of Basic Clusters in Super Clusters from Hybrid in Barrel",
			hist_bins_NumBC_,hist_min_NumBC_,hist_max_NumBC_);
  	hist_IslandEB_SC_NumBC_ 
		= dbe_->book1D("hist_IslandEB_SC_NumBC_","# of Basic Clusters in Super Clusters from Island in Barrel",
			hist_bins_NumBC_,hist_min_NumBC_,hist_max_NumBC_);
  	hist_IslandEE_SC_NumBC_ 
		= dbe_->book1D("hist_IslandEE_SC_NumBC_","# of Basic Clusters in Super Clusters from Island in Endcap",
		hist_bins_NumBC_,hist_min_NumBC_,hist_max_NumBC_);

  	hist_HybridEB_SC_ET_ 
		= dbe_->book1D("hist_HybridEB_SC_ET_","ET of Super Clusters with Hybrid in Barrel",
			hist_bins_ET_,hist_min_ET_,hist_max_ET_);
  	hist_IslandEB_SC_ET_ 
		= dbe_->book1D("hist_IslandEB_SC_ET_","ET of Super Clusters with Island in Barrel",
			hist_bins_ET_,hist_min_ET_,hist_max_ET_);
  	hist_IslandEE_SC_ET_ 
		= dbe_->book1D("hist_IslandEE_SC_ET_","ET of Super Clusters with Island in Endcap",
			hist_bins_ET_,hist_min_ET_,hist_max_ET_);

  	hist_HybridEB_SC_Eta_ 
		= dbe_->book1D("hist_HybridEB_SC_Eta_","Eta of Super Clusters with Hybrid in Barrel",
			hist_bins_Eta_,hist_min_Eta_,hist_max_Eta_);
  	hist_IslandEB_SC_Eta_ 
		= dbe_->book1D("hist_IslandEB_SC_Eta_","Eta of Super Clusters with Island in Barrel",
			hist_bins_Eta_,hist_min_Eta_,hist_max_Eta_);
 	hist_IslandEE_SC_Eta_ 
		= dbe_->book1D("hist_IslandEE_SC_Eta_","Eta of Super Clusters with Island in Endcap",
			hist_bins_Eta_,hist_min_Eta_,hist_max_Eta_);

  	hist_HybridEB_SC_Phi_ 
		= dbe_->book1D("hist_HybridEB_SC_Phi_","Phi of Super Clusters with Hybrid in Barrel",
			hist_bins_Phi_,hist_min_Phi_,hist_max_Phi_);
  	hist_IslandEB_SC_Phi_ 
		= dbe_->book1D("hist_IslandEB_SC_Phi_","Phi of Super Clusters with Island in Barrel",
			hist_bins_Phi_,hist_min_Phi_,hist_max_Phi_);
  	hist_IslandEE_SC_Phi_ 
		= dbe_->book1D("hist_IslandEE_SC_Phi_","Phi of Super Clusters with Island in Endcap",
			hist_bins_Phi_,hist_min_Phi_,hist_max_Phi_);

  	hist_HybridEB_SC_S1toS9_ 
		= dbe_->book1D("hist_HybridEB_SC_S1toS9_","S1/S9 of Super Clusters with Hybrid in Barrel",
			hist_bins_S1toS9_,hist_min_S1toS9_,hist_max_S1toS9_);
  	hist_IslandEB_SC_S1toS9_ 
		= dbe_->book1D("hist_IslandEB_SC_S1toS9_","S1/S9 of Super Clusters with Island in Barrel",
			hist_bins_S1toS9_,hist_min_S1toS9_,hist_max_S1toS9_);
  	hist_IslandEE_SC_S1toS9_ 
		= dbe_->book1D("hist_IslandEE_SC_S1toS9_","S1/S9 of Super Clusters with Island in Endcap",
			hist_bins_S1toS9_,hist_min_S1toS9_,hist_max_S1toS9_);

  	hist_HybridEB_SC_S25toE_ 
		= dbe_->book1D("hist_HybridEB_SC_S25toE_","S25/E of Super Clusters with Hybrid in Barrel",
			hist_bins_S25toE_,hist_min_S25toE_,hist_max_S25toE_);
  	hist_IslandEB_SC_S25toE_ 
		= dbe_->book1D("hist_IslandEB_SC_S25toE_","S25/E of Super Clusters with Island in Barrel",
			hist_bins_S25toE_,hist_min_S25toE_,hist_max_S25toE_);
  	hist_IslandEE_SC_S25toE_ 
		= dbe_->book1D("hist_IslandEE_SC_S25toE_","S25/E of Super Clusters with Island in Endcap",
			hist_bins_S25toE_,hist_min_S25toE_,hist_max_S25toE_);

  	hist_HybridEB_SC_EToverTruth_ 
		= dbe_->book1D("hist_HybridEB_SC_EToverTruth_","ET/True ET of Super Clusters with Hybrid in Barrel",	
			hist_bins_EToverTruth_,hist_min_EToverTruth_,hist_max_EToverTruth_);
  	hist_IslandEB_SC_EToverTruth_ 
		= dbe_->book1D("hist_IslandEB_SC_EToverTruth_","ET/True ET of Super Clusters with Island in Barrel",
			hist_bins_EToverTruth_,hist_min_EToverTruth_,hist_max_EToverTruth_);
  	hist_IslandEE_SC_EToverTruth_ 
		= dbe_->book1D("hist_IslandEE_SC_EToverTruth_","ET/True ET of Super Clusters with Island in Endcap",
			hist_bins_EToverTruth_,hist_min_EToverTruth_,hist_max_EToverTruth_);

  	hist_HybridEB_SC_deltaEta_ 
		= dbe_->book1D("hist_HybridEB_SC_deltaEta_","Eta-True Eta of Super Clusters with Hybrid in Barrel",
			hist_bins_deltaEta_,hist_min_deltaEta_,hist_max_deltaEta_);
  	hist_IslandEB_SC_deltaEta_ 
		= dbe_->book1D("hist_IslandEB_SC_deltaEta_","Eta-True Eta of Super Clusters with Island in Barrel",
			hist_bins_deltaEta_,hist_min_deltaEta_,hist_max_deltaEta_);
  	hist_IslandEE_SC_deltaEta_ 
		= dbe_->book1D("hist_IslandEE_SC_deltaEta_","Eta-True Eta of Super Clusters with Island in Endcap",
			hist_bins_deltaEta_,hist_min_deltaEta_,hist_max_deltaEta_);
}


void EgammaSuperClusters::analyze( const edm::Event& evt, const edm::EventSetup& es )
{
  	edm::Handle<reco::SuperClusterCollection> pHybridBarrelSuperClusters;
	evt.getByLabel(hybridBarrelSuperClusterCollection_, pHybridBarrelSuperClusters);
	if (!pHybridBarrelSuperClusters.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << hybridBarrelSuperClusterCollection_.label();
  	}

  	edm::Handle<reco::BasicClusterShapeAssociationCollection> pHybridBarrelClusterShapeAssociation;
	evt.getByLabel(hybridBarrelClusterShapeAssociation_, pHybridBarrelClusterShapeAssociation);
	if (!pHybridBarrelClusterShapeAssociation.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << hybridBarrelClusterShapeAssociation_.label();
  	}

  	const reco::SuperClusterCollection* hybridBarrelSuperClusters = pHybridBarrelSuperClusters.product();
  	hist_HybridEB_SC_Size_->Fill(hybridBarrelSuperClusters->size());

  	for(reco::SuperClusterCollection::const_iterator aClus = hybridBarrelSuperClusters->begin(); 
		aClus != hybridBarrelSuperClusters->end(); aClus++)
	{
		hist_HybridEB_SC_NumBC_->Fill(aClus->clustersSize());
    		hist_HybridEB_SC_ET_->Fill(aClus->energy()/std::cosh(aClus->position().eta()));
		hist_HybridEB_SC_Eta_->Fill(aClus->position().eta());
		hist_HybridEB_SC_Phi_->Fill(aClus->position().phi());

		const reco::ClusterShapeRef& tempClusterShape = pHybridBarrelClusterShapeAssociation->find(aClus->seed())->val;
		hist_HybridEB_SC_S1toS9_->Fill(tempClusterShape->eMax()/tempClusterShape->e3x3());
		hist_HybridEB_SC_S25toE_->Fill(tempClusterShape->e5x5()/aClus->energy());
  	}

  	edm::Handle<reco::SuperClusterCollection> pIslandBarrelSuperClusters;
	evt.getByLabel(islandBarrelSuperClusterCollection_, pIslandBarrelSuperClusters);
	if (!pIslandBarrelSuperClusters.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << islandBarrelSuperClusterCollection_.label();
  	}

  	edm::Handle<reco::BasicClusterShapeAssociationCollection> pIslandBarrelClusterShapeAssociation;
	evt.getByLabel(islandBarrelClusterShapeAssociation_, pIslandBarrelClusterShapeAssociation);
	if (!pIslandBarrelClusterShapeAssociation.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << islandBarrelClusterShapeAssociation_.label();
  	}

  	const reco::SuperClusterCollection* islandBarrelSuperClusters = pIslandBarrelSuperClusters.product();
 	hist_IslandEB_SC_Size_->Fill(islandBarrelSuperClusters->size());

  	for(reco::SuperClusterCollection::const_iterator aClus = islandBarrelSuperClusters->begin(); 
		aClus != islandBarrelSuperClusters->end(); aClus++)
	{
		hist_IslandEB_SC_NumBC_->Fill(aClus->clustersSize());
    		hist_IslandEB_SC_ET_->Fill(aClus->energy()/std::cosh(aClus->position().eta()));
		hist_IslandEB_SC_Eta_->Fill(aClus->position().eta());
		hist_IslandEB_SC_Phi_->Fill(aClus->position().phi());

		const reco::ClusterShapeRef& tempClusterShape = pIslandBarrelClusterShapeAssociation->find(aClus->seed())->val;
		hist_IslandEB_SC_S1toS9_->Fill(tempClusterShape->eMax()/tempClusterShape->e3x3());
		hist_IslandEB_SC_S25toE_->Fill(tempClusterShape->e5x5()/aClus->energy());
  	}

  	edm::Handle<reco::SuperClusterCollection> pIslandEndcapSuperClusters;
	evt.getByLabel(islandEndcapSuperClusterCollection_, pIslandEndcapSuperClusters);
	if (!pIslandEndcapSuperClusters.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << islandEndcapSuperClusterCollection_.label();
  	}

  	edm::Handle<reco::BasicClusterShapeAssociationCollection> pIslandEndcapClusterShapeAssociation;
	evt.getByLabel(islandEndcapClusterShapeAssociation_, pIslandEndcapClusterShapeAssociation);
	if (!pIslandEndcapClusterShapeAssociation.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << islandEndcapClusterShapeAssociation_.label();
  	}

  	const reco::SuperClusterCollection* islandEndcapSuperClusters = pIslandEndcapSuperClusters.product();
  	hist_IslandEE_SC_Size_->Fill(islandEndcapSuperClusters->size());

  	for(reco::SuperClusterCollection::const_iterator aClus = islandEndcapSuperClusters->begin(); 
		aClus != islandEndcapSuperClusters->end(); aClus++)
	{
		hist_IslandEE_SC_NumBC_->Fill(aClus->clustersSize());
    		hist_IslandEE_SC_ET_->Fill(aClus->energy()/std::cosh(aClus->position().eta()));
		hist_IslandEE_SC_Eta_->Fill(aClus->position().eta());
		hist_IslandEE_SC_Phi_->Fill(aClus->position().phi());

		const reco::ClusterShapeRef& tempClusterShape = pIslandEndcapClusterShapeAssociation->find(aClus->seed())->val;
		hist_IslandEE_SC_S1toS9_->Fill(tempClusterShape->eMax()/tempClusterShape->e3x3());
		hist_IslandEE_SC_S25toE_->Fill(tempClusterShape->e5x5()/aClus->energy());
  	}

 	edm::Handle<edm::HepMCProduct> pMCTruth ;
	evt.getByLabel(MCTruthCollection_, pMCTruth);
	if (!pMCTruth.isValid()) {
	  edm::LogError("EgammaSuperClusters") << "Error! can't get collection with label " 
					       << MCTruthCollection_.label();
  	}

	const HepMC::GenEvent* genEvent = pMCTruth->GetEvent();
  	for(HepMC::GenEvent::particle_const_iterator currentParticle = genEvent->particles_begin(); 
		currentParticle != genEvent->particles_end(); currentParticle++ )
  	{
	  	if((*currentParticle)->status()==1) 
		{
			HepMC::FourVector vtx = (*currentParticle)->production_vertex()->position();
			double phiTrue = (*currentParticle)->momentum().phi();
			double etaTrue = ecalEta((*currentParticle)->momentum().eta(), vtx.z()/10., vtx.perp()/10.);
			double etTrue  = (*currentParticle)->momentum().e()/cosh(etaTrue);

			if(std::fabs(etaTrue) < 1.479)
			{
				{
					double etaCurrent, etaFound = 0;
					double phiCurrent;
					double etCurrent,  etFound  = 0;

					double closestParticleDistance = 999; 
				
					for(reco::SuperClusterCollection::const_iterator aClus = hybridBarrelSuperClusters->begin(); 
						aClus != hybridBarrelSuperClusters->end(); aClus++)
					{
						etaCurrent = 	aClus->position().eta();
						phiCurrent = 	aClus->position().phi();
						etCurrent  =  aClus->energy()/std::cosh(etaCurrent);
											
						double deltaR = std::sqrt(std::pow(etaCurrent-etaTrue,2)+std::pow(phiCurrent-phiTrue,2));

						if(deltaR < closestParticleDistance)
						{
							etFound  = etCurrent;
							etaFound = etaCurrent;
							closestParticleDistance = deltaR;
						}
					}
					
					if(closestParticleDistance < 0.3)
					{
						hist_HybridEB_SC_EToverTruth_->Fill(etFound/etTrue);
						hist_HybridEB_SC_deltaEta_->Fill(etaFound-etaTrue);
					}
				}
				{
					double etaCurrent, etaFound = 0;
					double phiCurrent;
					double etCurrent,  etFound  = 0;

					double closestParticleDistance = 999; 
				
				  	for(reco::SuperClusterCollection::const_iterator aClus = islandBarrelSuperClusters->begin(); 
						aClus != islandBarrelSuperClusters->end(); aClus++)
					{
						etaCurrent = 	aClus->position().eta();
						phiCurrent = 	aClus->position().phi();
						etCurrent  =  aClus->energy()/std::cosh(etaCurrent);

						double deltaR = std::sqrt(std::pow(etaCurrent-etaTrue,2)+std::pow(phiCurrent-phiTrue,2)); 

						if(deltaR < closestParticleDistance)
						{
							etFound  = etCurrent;
							etaFound = etaCurrent;
							closestParticleDistance = deltaR;
						}
					}
					
					if(closestParticleDistance < 0.3)
					{
						hist_IslandEB_SC_EToverTruth_->Fill(etFound/etTrue);
						hist_IslandEB_SC_deltaEta_->Fill(etaFound-etaTrue);
					}
				}
			}
			else
			{
				double etaCurrent, etaFound = 0;
				double phiCurrent;
				double etCurrent,  etFound  = 0;

				double closestParticleDistance = 999; 

			  	for(reco::SuperClusterCollection::const_iterator aClus = islandEndcapSuperClusters->begin(); 
					aClus != islandEndcapSuperClusters->end(); aClus++)
				{
					etaCurrent = 	aClus->position().eta();
					phiCurrent = 	aClus->position().phi();
					etCurrent  =  aClus->energy()/std::cosh(etaCurrent);

					double deltaR = std::sqrt(std::pow(etaCurrent-etaTrue,2)+std::pow(phiCurrent-phiTrue,2)); 

					if(deltaR < closestParticleDistance)
					{
						etFound  = etCurrent;
						etaFound = etaCurrent;
						closestParticleDistance = deltaR;
					}
				}
				
				if(closestParticleDistance < 0.3)
				{
					hist_IslandEE_SC_EToverTruth_->Fill(etFound/etTrue);
					hist_IslandEE_SC_deltaEta_->Fill(etaFound-etaTrue);
				}
			}
		}
	}
}

void EgammaSuperClusters::endJob()
{
	if (outputFile_.size() != 0 && dbe_) dbe_->save(outputFile_);
}

float EgammaSuperClusters::ecalEta(float EtaParticle , float Zvertex, float plane_Radius)
{  
	const float R_ECAL           = 136.5;
	const float Z_Endcap         = 328.0;
	const float etaBarrelEndcap  = 1.479;

	if(EtaParticle != 0.)
	{
		float Theta = 0.0  ;
		float ZEcal = (R_ECAL-plane_Radius)*sinh(EtaParticle)+Zvertex;

		if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
		if(Theta<0.0) Theta = Theta+Geom::pi() ;

		float ETA = - log(tan(0.5*Theta));

		if( fabs(ETA) > etaBarrelEndcap )
		{
			float Zend = Z_Endcap ;
			if(EtaParticle<0.0 )  Zend = -Zend ;
			float Zlen = Zend - Zvertex ;
			float RR = Zlen/sinh(EtaParticle);
			Theta = atan((RR+plane_Radius)/Zend);
			if(Theta<0.0) Theta = Theta+Geom::pi() ;
			ETA = - log(tan(0.5*Theta));
		}
		
		return ETA;
	}
	else
	{
		edm::LogWarning("")  << "[EgammaSuperClusters::ecalEta] Warning: Eta equals to zero, not correcting" ;
		return EtaParticle;
	}
}
