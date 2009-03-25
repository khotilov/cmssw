// $Id: SuperCluster.cc,v 1.14 2009/03/24 10:47:55 ferriff Exp $
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
using namespace reco;

SuperCluster::SuperCluster( double energy, const math::XYZPoint& position ) :
  CaloCluster( energy, position ), rawEnergy_(-1.) {
}



SuperCluster::SuperCluster( double energy, const math::XYZPoint& position,
                            const CaloClusterPtr & seed,
                            const CaloClusterPtrVector& clusters,
			    double Epreshower, double phiWidth, double etaWidth) :
  CaloCluster(energy,position), rawEnergy_(-1.)
{
  phiWidth_ = phiWidth;
  etaWidth_ = etaWidth;
  seed_ = seed;
  preshowerEnergy_ = Epreshower;

  // set references to constituent basic clusters and update list of rechits
  for(CaloClusterPtrVector::const_iterator bcit  = clusters.begin();
                                            bcit != clusters.end();
                                          ++bcit) {
    clusters_.push_back( (*bcit) );

    // updated list of used hits
    const std::vector< std::pair<DetId, float> > & v1 = (*bcit)->hitsAndFractions();
    for( std::vector< std::pair<DetId, float> >::const_iterator diIt = v1.begin();
                                            diIt != v1.end();
                                           ++diIt ) {
      hitsAndFractions_.push_back( (*diIt) );
    } // loop over rechits
  } // loop over basic clusters

}



SuperCluster::SuperCluster( double energy, const math::XYZPoint& position,
                            const CaloClusterPtr & seed,
                            const CaloClusterPtrVector& clusters,
                            const CaloClusterPtrVector& preshowerClusters,
			    double Epreshower, double phiWidth, double etaWidth) :
  CaloCluster(energy,position), rawEnergy_(-1.)
{
  phiWidth_ = phiWidth;
  etaWidth_ = etaWidth;
  seed_ = seed;
  preshowerEnergy_ = Epreshower;

  // set references to constituent basic clusters and update list of rechits
  for(CaloClusterPtrVector::const_iterator bcit  = clusters.begin();
                                            bcit != clusters.end();
                                          ++bcit) {
    clusters_.push_back( (*bcit) );

    // updated list of used hits
    const std::vector< std::pair<DetId, float> > & v1 = (*bcit)->hitsAndFractions();
    for( std::vector< std::pair<DetId, float> >::const_iterator diIt = v1.begin();
                                            diIt != v1.end();
                                           ++diIt ) {
      hitsAndFractions_.push_back( (*diIt) );
    } // loop over rechits
  } // loop over basic clusters

  // set references to preshower clusters
  for(CaloClusterPtrVector::const_iterator pcit  = preshowerClusters.begin();
                                            pcit != preshowerClusters.end();
                                          ++pcit) {
    preshowerClusters_.push_back( (*pcit) );
  }
}



double SuperCluster::rawEnergy() const
{
  if (rawEnergy_<0) {
    rawEnergy_ = 0.;
    for(CaloClusterPtrVector::const_iterator bcItr = clustersBegin(); bcItr != clustersEnd(); bcItr++)
      {
	rawEnergy_ += (*bcItr)->energy();
      }
  }
  return rawEnergy_;
}
