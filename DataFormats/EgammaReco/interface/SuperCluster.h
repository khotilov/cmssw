#ifndef EgammaReco_SuperCluster_h
#define EgammaReco_SuperCluster_h
/** \class reco::SuperCluster SuperCluster.h DataFormats/EgammaReco/interface/SuperCluster.h
 *  
 * A SuperCluster reconstructed in the Electromagnetic Calorimeter
 * contains references to seed and constituent BasicClusters
 *
 * \author Luca Lista, INFN
 *
 * \version $Id: SuperCluster.h,v 1.18 2009/03/24 10:47:55 ferriff Exp $
 *
 */
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <Rtypes.h>

namespace reco {
  class SuperCluster : public CaloCluster {
  public:

    typedef math::XYZPoint Point;

    /// default constructor
    SuperCluster() : CaloCluster(0., Point(0.,0.,0.)), rawEnergy_(-1.) {}

    /// constructor defined by CaloCluster - will have to use setSeed and add() separately
    SuperCluster( double energy, const Point& position );

    SuperCluster( double energy, const Point& position,
                  const CaloClusterPtr & seed,
                  const CaloClusterPtrVector& clusters,
		  double Epreshower=0.,
		  double phiWidth=0., double etaWidth=0. );

    // to be merged in the previous one? -- FIXME
    SuperCluster( double energy, const Point& position,
                  const CaloClusterPtr & seed,
                  const CaloClusterPtrVector& clusters,
                  const CaloClusterPtrVector& preshowerClusters,
		  double Epreshower=0.,
		  double phiWidth=0., double etaWidth=0. );

    /// raw uncorrected energy (sum of energies of component BasicClusters)
    double rawEnergy() const;

    /// energy deposited in preshower 
    double preshowerEnergy() const { return preshowerEnergy_; }

    /// obtain phi and eta width of the Super Cluster
    double phiWidth() const { return phiWidth_; }
    double etaWidth() const { return etaWidth_; }

    //Assign new variables to supercluster
    void setPreshowerEnergy( double preshowerEnergy ) { preshowerEnergy_ = preshowerEnergy; };
    void setPhiWidth( double pw ) { phiWidth_ = pw; }
    void setEtaWidth( double ew ) { etaWidth_ = ew; }

    /// seed BasicCluster
    const CaloClusterPtr & seed() const { return seed_; }

    /// fist iterator over BasicCluster constituents
    CaloCluster_iterator clustersBegin() const { return clusters_.begin(); }

    /// last iterator over BasicCluster constituents
    CaloCluster_iterator clustersEnd() const { return clusters_.end(); }

    /// fist iterator over PreshowerCluster constituents
    CaloCluster_iterator preshowerClustersBegin() const { return preshowerClusters_.begin(); }

    /// last iterator over PreshowerCluster constituents
    CaloCluster_iterator preshowerClustersEnd() const { return preshowerClusters_.end(); }

    /// number of BasicCluster constituents
    size_t clustersSize() const { return clusters_.size(); }

    /// list of used xtals by DetId // now inherited by CaloCluster
    //std::vector<DetId> getHitsByDetId() const { return usedHits_; }

    /// set reference to seed BasicCluster
    void setSeed( const CaloClusterPtr & r ) { seed_ = r; }

    /// add reference to constituent BasicCluster
    void addCluster( const CaloClusterPtr & r ) { clusters_.push_back( r ); }

    /// add reference to constituent BasicCluster
    void addPreshowerCluster( const CaloClusterPtr & r ) { preshowerClusters_.push_back( r ); }

  private:

    /// reference to BasicCluster seed
    CaloClusterPtr seed_;

    /// references to BasicCluster constitunets
    CaloClusterPtrVector clusters_;

    /// references to BasicCluster constitunets
    CaloClusterPtrVector preshowerClusters_;

    /// used hits by detId - retrieved from BC constituents -- now inherited from CaloCluster
    //std::vector<DetId> usedHits_;

    double preshowerEnergy_;

    mutable double rawEnergy_;
    
    Double32_t phiWidth_;
    Double32_t etaWidth_;

  };

}
#endif
