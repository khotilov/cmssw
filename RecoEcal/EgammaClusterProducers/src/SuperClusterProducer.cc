// C/C++ headers
#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Reconstruction Classes
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// Class header file
#include "RecoEcal/EgammaClusterProducers/interface/SuperClusterProducer.h"


SuperClusterProducer::SuperClusterProducer(const edm::ParameterSet& ps)
{
  // The verbosity level
  std::string verbosityString = ps.getParameter<std::string>("VerbosityLevel");
  if      (verbosityString == "DEBUG")   verbosity = BremRecoveryClusterAlgo::DEBUG;
  else if (verbosityString == "WARNING") verbosity = BremRecoveryClusterAlgo::WARNING;
  else if (verbosityString == "INFO")    verbosity = BremRecoveryClusterAlgo::INFO;
  else                                   verbosity = BremRecoveryClusterAlgo::ERROR;

  endcapClusterProducer_ = ps.getParameter<std::string>("endcapClusterProducer");
  barrelClusterProducer_ = ps.getParameter<std::string>("barrelClusterProducer");

  endcapClusterCollection_ = ps.getParameter<std::string>("endcapClusterCollection");
  barrelClusterCollection_ = ps.getParameter<std::string>("barrelClusterCollection");

  endcapSuperclusterCollection_ = ps.getParameter<std::string>("endcapSuperclusterCollection");
  barrelSuperclusterCollection_ = ps.getParameter<std::string>("barrelSuperclusterCollection");

  barrelEtaSearchRoad_ = ps.getParameter<double>("barrelEtaSearchRoad");
  barrelPhiSearchRoad_ = ps.getParameter<double>("barrelPhiSearchRoad");
  endcapEtaSearchRoad_ = ps.getParameter<double>("endcapEtaSearchRoad");
  endcapPhiSearchRoad_ = ps.getParameter<double>("endcapPhiSearchRoad");
  seedEnergyThreshold_ = ps.getParameter<double>("seedEnergyThreshold");

  bremAlgo_p = new BremRecoveryClusterAlgo(barrelEtaSearchRoad_, barrelPhiSearchRoad_, 
					 endcapEtaSearchRoad_, endcapPhiSearchRoad_, 
					 seedEnergyThreshold_, verbosity);

  produces< reco::SuperClusterCollection >(endcapSuperclusterCollection_);
  produces< reco::SuperClusterCollection >(barrelSuperclusterCollection_);

  totalE = 0;
  noSuperClusters = 0;
  nEvt_ = 0;
}


SuperClusterProducer::~SuperClusterProducer()
{
  double averEnergy = totalE / noSuperClusters;
  edm::LogInfo("SuperClusterProducerInfo") << "-------------------------------------------------------";
  edm::LogInfo("SuperClusterProducerInfo") << "-------------------------------------------------------";
  edm::LogInfo("SuperClusterProducerInfo") << "average SuperCluster energy = " << averEnergy;
  edm::LogInfo("SuperClusterProducerInfo") << "-------------------------------------------------------";
  edm::LogInfo("SuperClusterProducerInfo") << "-------------------------------------------------------";
  delete bremAlgo_p;
}


void SuperClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  produceSuperclustersForECALPart(evt, endcapClusterProducer_, endcapClusterCollection_, endcapSuperclusterCollection_);
  produceSuperclustersForECALPart(evt, barrelClusterProducer_, barrelClusterCollection_, barrelSuperclusterCollection_);

  nEvt_++;
}

void SuperClusterProducer::produceSuperclustersForECALPart(edm::Event& evt, 
							   std::string clusterProducer, 
							   std::string clusterCollection,
							   std::string superclusterCollection)
{
  // get the cluster collection out and turn it to a BasicClusterRefVector:
  reco::BasicClusterRefVector *clusterRefVector = getClusterRefVector(evt, clusterProducer, clusterCollection);

  // run the brem recovery and get the SC collections
  reco::SuperClusterCollection *
    superclusterCollection_p = new reco::SuperClusterCollection(bremAlgo_p->makeSuperClusters(*clusterRefVector));

  std::auto_ptr<reco::SuperClusterCollection> superclusters_ap(new reco::SuperClusterCollection);
  superclusters_ap->assign(superclusterCollection_p->begin(), superclusterCollection_p->end());
  evt.put(superclusters_ap, superclusterCollection);

  reco::SuperClusterCollection::iterator it;
  for (it = superclusterCollection_p->begin(); it != superclusterCollection_p->end(); it++)
    {
      totalE += it->energy();
      noSuperClusters++;
    }
}


reco::BasicClusterRefVector *
SuperClusterProducer::getClusterRefVector(edm::Event& evt, std::string clusterProducer_, std::string clusterCollection_)
{  
  edm::Handle<reco::BasicClusterCollection> bccHandle;
  try
    {
      evt.getByLabel(clusterProducer_, clusterCollection_, bccHandle);
      if (!(bccHandle.isValid()))
	{
	  edm::LogError("SuperClusterProducerError") << "could not get a handle on the BasicCluster Collection!";
	  return 0;
	}
    } 
  catch ( cms::Exception& ex )
    {
      edm::LogError("SuperClusterProducerError") << "Error! can't get the product " << clusterCollection_.c_str(); 
      return 0;
    }

  const reco::BasicClusterCollection *clusterCollection_p = bccHandle.product();
  reco::BasicClusterRefVector *clusterRefVector_p = new reco::BasicClusterRefVector;
  for (unsigned int i = 0; i < clusterCollection_p->size(); i++)
    {
      clusterRefVector_p->push_back(reco::BasicClusterRef(bccHandle, i));
    }

  return clusterRefVector_p;
}                               




