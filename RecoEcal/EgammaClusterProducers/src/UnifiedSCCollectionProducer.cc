// C/C++ headers
#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
// Reconstruction Classes
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

// Class header file
#include "RecoEcal/EgammaClusterProducers/interface/UnifiedSCCollectionProducer.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"


/*
  UnifiedSCCollectionProducer:
  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  Takes  as  input the cleaned  and the  uncleaned collection of SC
  and  produces two collections of SC: one with the clean SC, but flagged
  such that with the algoID value one can identify the SC that are also
  in the unclean collection and a collection with the unclean only SC.
  This collection has the algoID enumeration of the SC altered
  such that:
  flags = 0   (cleanedOnly)     cluster is only in the cleaned collection
  flags = 100 (common)          cluster is common in both collections
  flags = 200 (uncleanedOnly)   cluster is only in the uncleaned collection

  In that way the user can get hold of objects from the
  -  cleaned   collection only if they choose flags <  200
  -  uncleaned collection only if they choose flags >= 100

  18 Aug 2010
  Nikolaos Rompotis and Chris Seez  - Imperial College London
  many thanks to David Wardrope, Shahram Rahatlou and Federico Ferri
*/


UnifiedSCCollectionProducer::UnifiedSCCollectionProducer(const edm::ParameterSet& ps)
{
  //
  // The debug level
  std::string debugString = ps.getParameter<std::string>("debugLevel");
  if      (debugString == "DEBUG")   debugL = HybridClusterAlgo::pDEBUG;
  else if (debugString == "INFO")    debugL = HybridClusterAlgo::pINFO;
  else                               debugL = HybridClusterAlgo::pERROR;

  // get the parameters
  // the cleaned collection:
  cleanBcCollection_ = ps.getParameter<edm::InputTag>("cleanBcCollection");
  cleanScCollection_ = ps.getParameter<edm::InputTag>("cleanScCollection");
  // the uncleaned collection
  uncleanBcCollection_ = ps.getParameter<edm::InputTag>("uncleanBcCollection");
  uncleanScCollection_ = ps.getParameter<edm::InputTag>("uncleanScCollection");
  // the names of the products to be produced:
  //
  // the clean collection: this is as it was before, but labeled
  bcCollection_ = ps.getParameter<std::string>("bcCollection");
  scCollection_ = ps.getParameter<std::string>("scCollection");
  // the unclean only collection: SC unique to the unclean collection
  bcCollectionUncleanOnly_ = ps.getParameter<std::string>("bcCollectionUncleanOnly");
  scCollectionUncleanOnly_ = ps.getParameter<std::string>("scCollectionUncleanOnly");
  // the products:
  produces< reco::BasicClusterCollection >(bcCollection_);
  produces< reco::SuperClusterCollection >(scCollection_);
  produces< reco::BasicClusterCollection >(bcCollectionUncleanOnly_);
  produces< reco::SuperClusterCollection >(scCollectionUncleanOnly_);
  
}


UnifiedSCCollectionProducer::~UnifiedSCCollectionProducer() {;}


void UnifiedSCCollectionProducer::produce(edm::Event& evt, 
				      const edm::EventSetup& es)
{
  if (debugL <= HybridClusterAlgo::pINFO)
    edm::LogInfo("UnifiedSC")<< ">>>>> Entering UnifiedSCCollectionProducer <<<<<";
  // get the input collections
  // __________________________________________________________________________
  //
  // cluster collections:
  edm::Handle<reco::BasicClusterCollection> pCleanBC;
  edm::Handle<reco::SuperClusterCollection> pCleanSC;
  //
  edm::Handle<reco::BasicClusterCollection> pUncleanBC;
  edm::Handle<reco::SuperClusterCollection> pUncleanSC;
  // clean collections ________________________________________________________
  //evt.getByLabel(cleanBcCollection_, pCleanBC);
  //if (!(pCleanBC.isValid())) 
  //  {
  //    if (debugL <= HybridClusterAlgo::pINFO)
  //      edm::LogInfo("UnifiedSC") << "could not handle clean basic clusters";
  //    return;
  //  }
  //const  reco::BasicClusterCollection cleanBS = *(pCleanBC.product());
  //
  evt.getByLabel(cleanScCollection_, pCleanSC);
  if (!(pCleanSC.isValid())) 
    {
      if (debugL <= HybridClusterAlgo::pINFO)
	std::cout << "could not handle clean super clusters" << std::endl;
      return;
    }
  //const  reco::SuperClusterCollection cleanSC = *(pCleanSC.product());

  // unclean collections ______________________________________________________
  //evt.getByLabel(uncleanBcCollection_, pUncleanBC);
  //if (!(pUncleanBC.isValid())) 
  //  {
  //    if (debugL <= HybridClusterAlgo::pINFO)
  //      edm::LogInfo("UnifiedSC")<<"could not handle unclean Basic Clusters!";
  //    return;
  //  }
  //const  reco::BasicClusterCollection uncleanBC = *(pUncleanBC.product());
  //
  evt.getByLabel(uncleanScCollection_, pUncleanSC);
  if (!(pUncleanSC.isValid())) 
    {
      if (debugL <= HybridClusterAlgo::pINFO)
	edm::LogInfo("UnifiedSC")<< "could not handle unclean super clusters!" ;
      return;
    }
  //const  reco::SuperClusterCollection uncleanSC = *(pUncleanSC.product());
  // collections are all taken now ____________________________________________
  //
  //
  // the collections to be produced ___________________________________________
  reco::BasicClusterCollection basicClusters;
  reco::SuperClusterCollection superClusters;
  //
  reco::BasicClusterCollection basicClustersUncleanOnly;
  reco::SuperClusterCollection superClustersUncleanOnly;
  //
  // run over the uncleaned SC and check how many of them are matched to 
  // the cleaned ones
  // if you find a matched one, then keep the info that it is matched 
  //    along with which clean SC was matched + its basic clusters
  // if you find an unmatched one, keep the info and store its basic clusters
  //
  // 
  int uncleanSize = (int) pUncleanSC->size();
  int cleanSize = (int) pCleanSC->size();
  if (debugL <= HybridClusterAlgo::pDEBUG)
    LogDebug("UnifiedSC") << "Size of Clean Collection: " << cleanSize 
			  << ", uncleanSize: " << uncleanSize;


  // keep the indices
  std::vector<int> inUncleanOnlyInd;      // counting the unclean
  std::vector<int> inCleanInd;            // counting the unclean
  std::vector<int> inCleanOnlyInd;        // counting the clean
  std::vector<DetId> scUncleanSeedDetId;  // counting the unclean
  std::vector<DetId> scCleanSeedDetId;    // counting the clean
  // ontains the index of the SC that owns that BS
  // first basic cluster index, second: 0 for unclean and 1 for clean
  std::vector< std::pair<int, int> > basicClusterOwner; 
  std::vector< std::pair<int, int> > basicClusterOwnerUncleanOnly; 
  // if this basic cluster is a seed it is 1
  std::vector<int> uncleanBasicClusterIsSeed;

  // loop over unclean SC _____________________________________________________
  for (int isc =0; isc< uncleanSize; ++isc) {
    reco::SuperClusterRef unscRef( pUncleanSC, isc);    
    const std::vector< std::pair<DetId, float> > & uhits = unscRef->hitsAndFractions();
    int uhitsSize = (int) uhits.size();
    bool foundTheSame = false;
    for (int jsc=0; jsc < cleanSize; ++jsc) { // loop over the cleaned SC
      reco::SuperClusterRef cscRef( pCleanSC, jsc );
      const std::vector<std::pair<DetId,float> > & chits = cscRef->hitsAndFractions();
      int chitsSize = (int) chits.size();
      foundTheSame = false;
      if (unscRef->seed()->seed()==cscRef->seed()->seed() && chitsSize == uhitsSize) { 
	// if the clusters are exactly the same then because the clustering
	// algorithm works in a deterministic way, the order of the rechits
	// will be the same
	for (int i=0; i< chitsSize; ++i) {
	  if (uhits[i].first != chits[i].first ) { break;}
	}
	foundTheSame = true;
      }
      if (foundTheSame) { // ok you have found it:
	// this supercluster belongs to both collections
	inUncleanOnlyInd.push_back(0);
	inCleanInd.push_back(jsc); // keeps the index of the clean SC
	scUncleanSeedDetId.push_back(unscRef->seed()->seed());
	//
	// keep its basic clusters:
	for (reco::CaloCluster_iterator bciter = unscRef->clustersBegin(); bciter != unscRef->clustersEnd(); ++bciter) {
	  // the basic clusters
	  //reco::CaloClusterPtr myclusterptr = *bciter;
	  //reco::CaloCluster mycluster = *myclusterptr;
	  basicClusters.push_back(**bciter);
	  // index of the unclean SC
	  basicClusterOwner.push_back( std::make_pair(isc,0) ); 
	} 
	break; // break the loop over unclean sc
      }
    }
    if (not foundTheSame) { // this SC is only in the unclean collection
      // mark it as unique in the uncleaned
      inUncleanOnlyInd.push_back(1);
      scUncleanSeedDetId.push_back(unscRef->seed()->seed());
      // keep all its basic clusters
      for (reco::CaloCluster_iterator bciter = unscRef->clustersBegin(); bciter != unscRef->clustersEnd(); ++bciter) {
	// the basic clusters
	//reco::CaloClusterPtr myclusterptr = *bciter;
	//reco::CaloCluster mycluster = *myclusterptr;
	basicClustersUncleanOnly.push_back(**bciter);
	basicClusterOwnerUncleanOnly.push_back( std::make_pair(isc,0) );
      }
    }
  } // loop over the unclean SC _______________________________________________
  //
  int inCleanSize = (int) inCleanInd.size();
  //
  // loop over the clean SC, check that are not in common with the unclean
  // ones and then store their SC as before ___________________________________
  for (int jsc =0; jsc< cleanSize; ++jsc) {
    // check whether this index is already in the common collection
    bool takenAlready = false;
    for (int j=0; j< inCleanSize; ++j) {
      if (jsc == inCleanInd[j]) { takenAlready = true ;break;}
    }
    if (takenAlready) {
      inCleanOnlyInd.push_back(0);
      scCleanSeedDetId.push_back(DetId(0));
      continue;
    }
    inCleanOnlyInd.push_back(1);
    reco::SuperClusterRef cscRef( pCleanSC, jsc );
    scCleanSeedDetId.push_back(cscRef->seed()->seed());
    for (reco::CaloCluster_iterator bciter = cscRef->clustersBegin(); bciter != cscRef->clustersEnd(); ++bciter) {
      // the basic clusters
      //reco::CaloClusterPtr myclusterptr = *bciter;
      //reco::CaloCluster mycluster = *myclusterptr;
      basicClusters.push_back(**bciter);
      basicClusterOwner.push_back( std::make_pair(jsc,1) );
    }
  } // end loop over clean SC _________________________________________________
  //
  //
  // at this point we have the basic cluster collection ready
  // 
  int bcSize = (int) basicClusters.size();
  int bcSizeUncleanOnly = (int) basicClustersUncleanOnly.size();
  if (debugL == HybridClusterAlgo::pDEBUG)
    LogDebug("UnifiedSC") << "Found cleaned SC: " << cleanSize 
			  <<  " uncleaned SC: "   << uncleanSize ;
  //
  // export the clusters to the event from the clean clusters
  std::auto_ptr< reco::BasicClusterCollection> 
    basicClusters_p(new reco::BasicClusterCollection);
  basicClusters_p->assign(basicClusters.begin(), basicClusters.end());
  edm::OrphanHandle<reco::BasicClusterCollection> bccHandle =  
    evt.put(basicClusters_p, bcCollection_);
  if (!(bccHandle.isValid())) {
    if (debugL <= HybridClusterAlgo::pINFO)
      edm::LogInfo("UnifiedSC")<< "could not handle the new BasicClusters!";
    return;
  }
  reco::BasicClusterCollection basicClustersProd = *bccHandle;
  if (debugL == HybridClusterAlgo::pDEBUG)
    edm::LogInfo("UnifiedSC")<< "Got the BasicClusters from the event again" << std::endl;
  //
  // export the clusters to the event: from the unclean only clusters
  std::auto_ptr< reco::BasicClusterCollection> 
    basicClustersUncleanOnly_p(new reco::BasicClusterCollection);
  basicClustersUncleanOnly_p->assign(basicClustersUncleanOnly.begin(), 
				     basicClustersUncleanOnly.end());
  edm::OrphanHandle<reco::BasicClusterCollection> bccHandleUncleanOnly =  
    evt.put(basicClustersUncleanOnly_p, bcCollectionUncleanOnly_);
  if (!(bccHandleUncleanOnly.isValid())) {
    if (debugL <= HybridClusterAlgo::pINFO)
      edm::LogInfo("UnifiedSC")<< "could not handle the new BasicClusters (Unclean Only)!" << std::endl;
    return;
  }
  reco::BasicClusterCollection basicClustersUncleanOnlyProd = *bccHandleUncleanOnly;
  if (debugL == HybridClusterAlgo::pDEBUG)
    edm::LogInfo("UnifiedSC")<< "Got the BasicClusters from the event again  (Unclean Only)" << std::endl;
  //

  // now we can build the SC collection
  //
  // start again from the unclean collection
  // all the unclean SC will become members of the new collection
  // with different algoIDs ___________________________________________________
  for (int isc=0; isc< uncleanSize; ++isc) {
    //std::cout << "working in ucl #" << isc << std::endl;
    reco::CaloClusterPtrVector clusterPtrVector;
    // the seed is the basic cluster with the highest energy
    reco::CaloClusterPtr seed; 
    if (inUncleanOnlyInd[isc] == 1) { // unclean SC Unique in Unclean
      for (int jbc=0; jbc< bcSizeUncleanOnly; ++jbc) {
	std::pair<int, int> theBcOwner = basicClusterOwnerUncleanOnly[jbc];
	if (theBcOwner.first == isc && theBcOwner.second == 0) {
	  reco::CaloClusterPtr currentClu=reco::CaloClusterPtr(bccHandleUncleanOnly,jbc);
	  clusterPtrVector.push_back(currentClu);
	  if (scUncleanSeedDetId[isc] == currentClu->seed()) {
	    seed = currentClu;
	  }
	}
      }

    }
    else { // unclean SC common in clean and unclean
      for (int jbc=0; jbc< bcSize; ++jbc) {
	std::pair<int, int> theBcOwner = basicClusterOwner[jbc];
	if (theBcOwner.first == isc && theBcOwner.second == 0) {
	  reco::CaloClusterPtr currentClu=reco::CaloClusterPtr(bccHandle,jbc);
	  clusterPtrVector.push_back(currentClu);
	  if (scUncleanSeedDetId[isc] == currentClu->seed()) {
	    seed = currentClu;
	  }
	}
      }
    }
    //std::cout << "before getting the uncl" << std::endl;
    reco::SuperClusterRef unscRef( pUncleanSC, isc ); 
    reco::SuperCluster newSC(unscRef->energy(), unscRef->position(), 
			     seed, clusterPtrVector );
    // now set the algoID for this SC again
    if (inUncleanOnlyInd[isc] == 1) {
            // set up the quality to unclean only .............
            newSC.setFlags(reco::CaloCluster::uncleanOnly);
            superClustersUncleanOnly.push_back(newSC);
    }
    else {
            // set up the quality to common  .............
            newSC.setFlags(reco::CaloCluster::common);
            superClusters.push_back(newSC);
    }
    // now you can store your SC

  } // end loop over unclean SC _______________________________________________
  //  flags numbering scheme
  //  flags =   0 = cleanedOnly     cluster is only in the cleaned collection
  //  flags = 100 = common          cluster is common in both collections
  //  flags = 200 = uncleanedOnly   cluster is only in the uncleaned collection

  // now loop over the clean SC and do the same but now you have to avoid the
  // the duplicated ones ______________________________________________________
  for (int jsc=0; jsc< cleanSize; ++jsc) {
    //std::cout << "working in cl #" << jsc << std::endl;
    // check that the SC is not in the unclean collection
    if (inCleanOnlyInd[jsc] == 0) continue;
    reco::CaloClusterPtrVector clusterPtrVector;
    // the seed is the basic cluster with the highest energy
    reco::CaloClusterPtr seed; 
    for (int jbc=0; jbc< bcSize; ++jbc) {
      std::pair<int, int> theBcOwner = basicClusterOwner[jbc];
      if (theBcOwner.first == jsc && theBcOwner.second == 1) {
	reco::CaloClusterPtr currentClu=reco::CaloClusterPtr(bccHandle,jbc);
	clusterPtrVector.push_back(currentClu);
	if (scCleanSeedDetId[jsc] == currentClu->seed()) {
	  seed = currentClu;
	}
      }
    }
    reco::SuperClusterRef cscRef( pCleanSC, jsc ); 
    reco::SuperCluster newSC(cscRef->energy(), cscRef->position(),
			     seed, clusterPtrVector );
    newSC.setFlags(reco::CaloCluster::cleanOnly);

    // add it to the collection:
    superClusters.push_back(newSC);

  } // end loop over clean SC _________________________________________________

  if (debugL == HybridClusterAlgo::pDEBUG)
    LogDebug("UnifiedSC")<< "New SC collection was created";

  std::auto_ptr< reco::SuperClusterCollection> 
    superClusters_p(new reco::SuperClusterCollection);
  superClusters_p->assign(superClusters.begin(), superClusters.end());

  evt.put(superClusters_p, scCollection_);
  if (debugL == HybridClusterAlgo::pDEBUG)
    LogDebug("UnifiedSC") << "Clusters (Basic/Super) added to the Event! :-)";

  std::auto_ptr< reco::SuperClusterCollection> 
    superClustersUncleanOnly_p(new reco::SuperClusterCollection);
  superClustersUncleanOnly_p->assign(superClustersUncleanOnly.begin(), 
				     superClustersUncleanOnly.end());

  evt.put(superClustersUncleanOnly_p, scCollectionUncleanOnly_);

  // ----- debugging ----------
  // print the new collection SC quantities
  if (debugL == HybridClusterAlgo::pDEBUG) {
    // print out the clean collection SC
    LogDebug("UnifiedSC") << "Clean Collection SC ";
    for (int i=0; i < cleanSize; ++i) {
      reco::SuperClusterRef cscRef( pCleanSC, i );
      std::cout << " >>> clean    #" << i << "; Energy: " << cscRef->energy()
		<< " eta: " << cscRef->eta() 
		<< " sc seed detid: " << cscRef->seed()->seed().rawId()
		<< std::endl;
    }
    // the unclean SC
    LogDebug("UnifiedSC") << "Unclean Collection SC ";
    for (int i=0; i < uncleanSize; ++i) {
      reco::SuperClusterRef uscRef( pUncleanSC, i );
      LogDebug("UnifiedSC") << " >>> unclean  #" << i << "; Energy: " << uscRef->energy()
			    << " eta: " << uscRef->eta() 
			    << " sc seed detid: " << uscRef->seed()->seed().rawId();
    }
    // the new collection
    LogDebug("UnifiedSC")<< "The new SC clean collection with size "<< superClusters.size()  << std::endl;
    
    int new_unclean = 0, new_clean=0;
    for (int i=0; i < (int) superClusters.size(); ++i) {
      const reco::SuperCluster & nsc = superClusters[i];
      LogDebug("UnifiedSC") << "SC was got" << std::endl
       << " ---> energy: " << nsc.energy() << std::endl
       << " ---> eta: " << nsc.eta() << std::endl
       << " ---> inClean: " << nsc.isInClean() << std::endl
       << " ---> id: " << nsc.seed()->seed().rawId() << std::endl
       << " >>> newSC    #" << i << "; Energy: " << nsc.energy()
       << " eta: " << nsc.eta()  << " isClean=" 
       << nsc.isInClean() << " isUnclean=" << nsc.isInUnclean()
       << " sc seed detid: " << nsc.seed()->seed().rawId();

      if (nsc.isInUnclean()) ++new_unclean;
      if (nsc.isInClean()) ++new_clean;
    }
    LogDebug("UnifiedSC")<< "The new SC unclean only collection with size "<< superClustersUncleanOnly.size();
    for (int i=0; i < (int) superClustersUncleanOnly.size(); ++i) {
      const reco::SuperCluster nsc = superClustersUncleanOnly[i];
      LogDebug ("UnifiedSC") << " >>> newSC    #" << i << "; Energy: " << nsc.energy()
		<< " eta: " << nsc.eta()  << " isClean=" 
		<< nsc.isInClean() << " isUnclean=" << nsc.isInUnclean()
		<< " sc seed detid: " << nsc.seed()->seed().rawId();
      if (nsc.isInUnclean()) ++new_unclean;
      if (nsc.isInClean()) ++new_clean;      
    }
    if ( (new_unclean != uncleanSize) || (new_clean != cleanSize) ) {
      LogDebug("UnifiedSC") << ">>>>!!!!!! MISMATCH: new unclean/ old unclean= " 
		<< new_unclean << " / " << uncleanSize 
	        << ", new clean/ old clean" << new_clean << " / " << cleanSize;
    }
  }
}




