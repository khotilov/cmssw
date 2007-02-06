// $Id: SuperCluster.cc,v 1.6 2006/09/12 23:53:40 rahatlou Exp $
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/EgammaReco/interface/ClusterPi0Discriminator.h"
using namespace reco;

SuperCluster::SuperCluster( double energy, const math::XYZPoint& position ) :
  EcalCluster( energy, position ) {
}

SuperCluster::SuperCluster( double energy, const math::XYZPoint& position,
                            const BasicClusterRef & seed,
                            const BasicClusterRefVector& clusters ) :
  EcalCluster(energy,position)
{

  seed_ = seed;

  // set references to constituent basic clusters and update list of rechits
  for(BasicClusterRefVector::const_iterator bcit  = clusters.begin();
                                            bcit != clusters.end();
                                          ++bcit) {
    clusters_.push_back( (*bcit) );

    // updated list of used hits
    const std::vector<DetId> & v1 = (*bcit)->getHitsByDetId();
    for( std::vector<DetId>::const_iterator diIt = v1.begin();
                                            diIt != v1.end();
                                           ++diIt ) {
      usedHits_.push_back( (*diIt) );
    } // loop over rechits
  } // loop over basic clusters

}
