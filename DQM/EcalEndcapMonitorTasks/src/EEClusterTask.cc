/*
 * \file EEClusterTask.cc
 *
 * $Date: 2007/04/05 13:56:49 $
 * $Revision: 1.2 $
 * \author G. Della Ricca
 * \author E. Di Marco
 *
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"

#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <DQM/EcalEndcapMonitorTasks/interface/EEClusterTask.h>

#include <TLorentzVector.h>

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

EEClusterTask::EEClusterTask(const ParameterSet& ps){

  init_ = false;

  // get hold of back-end interface
  dbe_ = Service<DaqMonitorBEInterface>().operator->();

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", true);

  // parameters...
  islandBarrelBasicClusterCollection_ = ps.getParameter<edm::InputTag>("islandBarrelBasicClusterCollection");

  islandBarrelSuperClusterCollection_ = ps.getParameter<edm::InputTag>("islandBarrelSuperClusterCollection");

  hybridSuperClusterCollection_ = ps.getParameter<edm::InputTag>("hybridSuperClusterCollection");

  hybridBarrelClusterShapeAssociation_ = ps.getParameter<edm::InputTag>("hybridBarrelClusterShapeAssociation");

  // histograms...
  meIslBEne_ = 0;
  meIslBNum_ = 0;
  meIslBCry_ = 0;

  meIslBEneMap_ = 0;
  meIslBNumMap_ = 0;
  meIslBETMap_  = 0;
  meIslBCryMap_ = 0;

  meIslSEne_ = 0;
  meIslSNum_ = 0;
  meIslSSiz_ = 0;

  meIslSEneMap_ = 0;
  meIslSNumMap_ = 0;
  meIslSETMap_ = 0;
  meIslSSizMap_ = 0;

  meHybS1toE_  = 0;
  meInvMass_ = 0;

}

EEClusterTask::~EEClusterTask(){

}

void EEClusterTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEClusterTask");
    dbe_->rmdir("EcalEndcap/EEClusterTask");
  }

}

void EEClusterTask::setup(void){

  init_ = true;

  Char_t histo[200];

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEClusterTask");

    sprintf(histo, "EECLT island basic cluster energy");
    meIslBEne_ = dbe_->book1D(histo, histo, 100, 0., 150.);

    sprintf(histo, "EECLT island basic cluster number");
    meIslBNum_ = dbe_->book1D(histo, histo, 100, 0., 100.);

    sprintf(histo, "EECLT island basic cluster crystals");
    meIslBCry_ = dbe_->book1D(histo, histo, 100, 0., 100.);

    sprintf(histo, "EECLT island basic cluster energy map");
    meIslBEneMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 500., "s");

    sprintf(histo, "EECLT island basic cluster number map");
    meIslBNumMap_ = dbe_->book2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI);

    sprintf(histo, "EECLT island basic cluster ET map");
    meIslBETMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 500., "s");

    sprintf(histo, "EECLT island basic cluster size map");
    meIslBCryMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 100., "s");

    sprintf(histo, "EECLT island super cluster energy");
    meIslSEne_ = dbe_->book1D(histo, histo, 100, 0., 150.);

    sprintf(histo, "EECLT island super cluster number");
    meIslSNum_ = dbe_->book1D(histo, histo, 50, 0., 50.);

    sprintf(histo, "EECLT island super cluster size");
    meIslSSiz_ = dbe_->book1D(histo, histo, 10, 0., 10.);

    sprintf(histo, "EECLT island super cluster energy map");
    meIslSEneMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 500., "s");

    sprintf(histo, "EECLT island super cluster number map");
    meIslSNumMap_ = dbe_->book2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI);

    sprintf(histo, "EECLT island super cluster ET map");
    meIslSETMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 500., "s");

    sprintf(histo, "EECLT island super cluster size map");
    meIslSSizMap_ = dbe_->bookProfile2D(histo, histo, 34, -1.479, 1.479, 72, -M_PI, M_PI, 100, 0., 500., "s");

    sprintf(histo, "EECLT hybrid S1toE");
    meHybS1toE_ = dbe_->book1D(histo, histo, 50, 0., 1.);

    sprintf(histo, "EECLT dicluster invariant mass");
    meInvMass_ = dbe_->book1D(histo, histo, 50, 60., 120.);

  }

}

void EEClusterTask::cleanup(void){

  if ( ! enableCleanup_ ) return;

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEClusterTask");

    if ( meIslBEne_ ) dbe_->removeElement( meIslBEne_->getName() );
    meIslBEne_ = 0;

    if ( meIslBNum_ ) dbe_->removeElement( meIslBNum_->getName() );
    meIslBNum_ = 0;

    if ( meIslBCry_ ) dbe_->removeElement( meIslBCry_->getName() );
    meIslBCry_ = 0;

    if ( meIslBEneMap_ ) dbe_->removeElement( meIslBEneMap_->getName() );
    meIslBEneMap_ = 0;

    if ( meIslBNumMap_ ) dbe_->removeElement( meIslBNumMap_->getName() );
    meIslBNumMap_ = 0;

    if ( meIslBETMap_ ) dbe_->removeElement( meIslBETMap_->getName() );
    meIslBETMap_ = 0;

    if ( meIslBCryMap_ ) dbe_->removeElement( meIslBCryMap_->getName() );
    meIslBCryMap_ = 0;

    if ( meIslSEne_ ) dbe_->removeElement( meIslSEne_->getName() );
    meIslSEne_ = 0;

    if ( meIslSNum_ ) dbe_->removeElement( meIslSNum_->getName() );
    meIslSNum_ = 0;

    if ( meIslSSiz_ ) dbe_->removeElement( meIslSSiz_->getName() );
    meIslSSiz_ = 0;

    if ( meIslSEneMap_ ) dbe_->removeElement( meIslSEneMap_->getName() );
    meIslSEneMap_ = 0;

    if ( meIslSNumMap_ ) dbe_->removeElement( meIslSNumMap_->getName() );
    meIslSNumMap_ = 0;

    if ( meIslSETMap_ ) dbe_->removeElement( meIslSETMap_->getName() );
    meIslSETMap_ = 0;

    if ( meIslSSizMap_ ) dbe_->removeElement( meIslSSizMap_->getName() );
    meIslSSizMap_ = 0;

    if ( meHybS1toE_ ) dbe_->removeElement( meHybS1toE_->getName() );
    meHybS1toE_ = 0;

    if ( meInvMass_ ) dbe_->removeElement( meInvMass_->getName() );
    meInvMass_ = 0;

  }

  init_ = false;

}

void EEClusterTask::endJob(void){

  LogInfo("EEClusterTask") << "analyzed " << ievt_ << " events";

  if ( init_ ) this->cleanup();

}

void EEClusterTask::analyze(const Event& e, const EventSetup& c){

  if ( ! init_ ) this->setup();

  ievt_++;

  // --- Get the Basic Clusters from Island and Hybrid Algorithms ---

  // --- Barrel "Island" Basic Clusters ---
  try {

    Handle<BasicClusterCollection> pIslandBarrelBasicClusters;
    e.getByLabel(islandBarrelBasicClusterCollection_, pIslandBarrelBasicClusters);

    Int_t nbcc = pIslandBarrelBasicClusters->size();

    meIslBNum_->Fill(float(nbcc));

    for ( BasicClusterCollection::const_iterator bclusterItr = pIslandBarrelBasicClusters->begin(); bclusterItr != pIslandBarrelBasicClusters->end(); ++bclusterItr ) {

      BasicCluster bcluster = *(bclusterItr);

      meIslBEne_->Fill(bcluster.energy());
      meIslBCry_->Fill(float(bcluster.getHitsByDetId().size()));

      meIslBEneMap_->Fill(bcluster.eta(), bcluster.phi(), bcluster.energy());
      meIslBNumMap_->Fill(bcluster.eta(), bcluster.phi() );
      meIslBCryMap_->Fill(bcluster.eta(), bcluster.phi(), float(bcluster.getHitsByDetId().size()) );
      meIslBETMap_->Fill(bcluster.eta(), bcluster.phi(), float(bcluster.energy()) * sin(bcluster.position().theta()) );

    }

  } catch ( exception& ex ) {
    LogWarning("EEClusterTask") << " BasicClusterCollection: " << islandBarrelBasicClusterCollection_ << " not in event.";
  }

  // --- Barrel "Island" Super Clusters ----
  try {

    Handle<SuperClusterCollection> pIslandBarrelSuperClusters;
    e.getByLabel(islandBarrelSuperClusterCollection_, pIslandBarrelSuperClusters);

    Int_t nscc = pIslandBarrelSuperClusters->size();

    meIslSNum_->Fill(float(nscc));

    for ( SuperClusterCollection::const_iterator sclusterItr = pIslandBarrelSuperClusters->begin(); sclusterItr != pIslandBarrelSuperClusters->end(); ++sclusterItr ) {

      SuperCluster scluster = *(sclusterItr);

      meIslSEne_->Fill(scluster.energy());
      meIslSSiz_->Fill(float(scluster.clustersSize()));

      meIslSEneMap_->Fill(scluster.eta(), scluster.phi(), scluster.energy());
      meIslSNumMap_->Fill(scluster.eta(), scluster.phi() );
      meIslSETMap_->Fill(scluster.eta(), scluster.phi(), float(scluster.energy()) * sin(scluster.position().theta()) );
      meIslSSizMap_->Fill(scluster.eta(), scluster.phi(), float(scluster.clustersSize()) );

    }

  } catch ( exception& ex ) {
    LogWarning("EEClusterTask") << " SuperClusterCollection: " << islandBarrelSuperClusterCollection_ << " not in event.";
  }

  // --- Barrel "Hybrid" Super Clusters ---
  try {

    Handle<SuperClusterCollection> pHybridSuperClusters;
    e.getByLabel(hybridSuperClusterCollection_, pHybridSuperClusters);
    Int_t nscc = pHybridSuperClusters->size();

    Handle<BasicClusterShapeAssociationCollection> pHybridBarrelClusterShapeAssociation;
    try	{
      e.getByLabel(hybridBarrelClusterShapeAssociation_, pHybridBarrelClusterShapeAssociation);
    }	catch ( cms::Exception& ex )	{
      LogWarning("EEClusterTask") << "Can't get collection with label "   << hybridBarrelClusterShapeAssociation_.label();
    }

    //    meHybSNum_->Fill(float(nscc));

    TLorentzVector sc1_p(0,0,0,0);
    TLorentzVector sc2_p(0,0,0,0);

    for(  SuperClusterCollection::const_iterator sCluster = pHybridSuperClusters->begin(); sCluster != pHybridSuperClusters->end(); ++sCluster ) {

      // seed
      const ClusterShapeRef& tempClusterShape = pHybridBarrelClusterShapeAssociation->find(sCluster->seed())->val;
      //       meHybS1toS9_->Fill(tempClusterShape->eMax()/tempClusterShape->e3x3());
      meHybS1toE_->Fill(tempClusterShape->eMax()/sCluster->energy());


//       // for each basic cluster evaluate the distance from the seed
//       if (sCluster->clustersSize()>1) {
// 	basicCluster_iterator bc;
// 	for(bc = sCluster->clustersBegin(); bc!=sCluster->clustersEnd(); bc++) {
// 	  Float_t dtheta=fabs( (*bc)->position().theta() - sCluster->seed()->position().theta() );
// 	  Float_t dphi=fabs( (*bc)->position().phi() - sCluster->seed()->position().phi() );
// 	  // exclude the seed...
// 	  if (dtheta!=0 && dphi!=0) {
// 	    meHybDTheta_->Fill( dtheta );
// 	    meHybDPhi_->Fill( dphi );

// 	    meHybEneVsDTheta_->Fill( dtheta, (*bc)->energy() );
// 	    meHybEneVsDPhi_->Fill( dphi, (*bc)->energy() );
// 	  }
// 	}
//       }

      // look for the two most energetic super clusters
      if (nscc>1) {
	if (sCluster->energy()>sc1_p.Energy()) {
	  sc2_p=sc1_p;
	  sc1_p.SetPtEtaPhiE(sCluster->energy()*sin(sCluster->position().theta()),
			     sCluster->eta(), sCluster->phi(), sCluster->energy());
	}
	else if (sCluster->energy()>sc2_p.Energy()) {
	  sc2_p.SetPtEtaPhiE(sCluster->energy()*sin(sCluster->position().theta()),
			     sCluster->eta(), sCluster->phi(), sCluster->energy());
	}
      }
    }
    // Get the invariant mass of the two most energetic super clusters
    if (nscc>1) {
      TLorentzVector sum = sc1_p+sc2_p;
      meInvMass_->Fill(sum.M());
    }

  } catch ( exception& ex ) {
    LogWarning("EEClusterTask") << " SuperClusterCollection: not in event.";
  }

}

