#ifndef DataFormats_ParticleFlowReco_PFCluster_h
#define DataFormats_ParticleFlowReco_PFCluster_h

#include "Math/GenVector/PositionVector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include <iostream>
#include <vector>



class PFClusterAlgo;

namespace reco {

  /**\class PFCluster
     \brief Particle flow cluster, see clustering algorithm in PFClusterAlgo
          
     \author Colin Bernet
     \date   July 2006
  */
  class PFCluster {
  public:

    /// type definition
    enum Type {
      TYPE_TOPOLOGICAL = 1, 
      TYPE_PF = 2 
    };


    typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;
  
    /// default constructor
    PFCluster();
  
    /// constructor
    PFCluster(unsigned id, int type);

    /// constructor
    PFCluster(unsigned id, int type, int layer, double energy,
	      double x, double y, double z );

    /// copy constructor
    PFCluster(const PFCluster& other);

    /// destructor
    ~PFCluster();

    /// resets clusters parameters
    void reset();
    
    /// add a given fraction of the rechit
    void addRecHitFraction( const reco::PFRecHitFraction& frac);
						
    //C this function will be moved to PFClusterAlgo
    /// \brief updates cluster info from rechit
    /// 
    /// algo = POSCALC_LIN (POSCALC_LOG) for linear (logarithmic) weighting 
    ///
    /// if linear, the rechit position is weighted by the rechit energy E
    ///
    /// if logarithmic, the rechit position is weighted by log(E/p1). 
    ///
    /// if p1 = -1, it is determined automatically 
    /// auto determination of p1 works only for ECAL and HCAL
    ///
    /// ncrystal is the number of crystals around the seed used in the 
    /// calculation. can be -1 (all), 5, or 9. 
/*     void calculatePosition( int algo,  */
/* 			    double p1 = 0,  */
/* 			    bool depcor = true,  */
/* 			    int  ncrystals = -1); */

    /// vector of rechit fractions
    const std::vector< reco::PFRecHitFraction >& recHitFractions() const 
      { return rechits_; }

    /// set cluster id
    void          setId(unsigned id) {id_ = id;} 
  
    /// cluster id
    unsigned      id() const {return id_;}
  
    /// cluster type
    int           type() const {return type_;}

    /// cluster layer, see PFLayer.h in this directory
    int           layer() const {return layer_;}          

    /// cluster energy
    double        energy() const {return energy_;}

    /// cluster position: cartesian
    const math::XYZPoint& positionXYZ() const {return posxyz_;}

    /// cluster position: rho, eta, phi
    const REPPoint&       positionREP() const {return posrep_;}

    /// calculates posrep_ once and for all
    void calculatePositionREP() {
      posrep_.SetCoordinates( posxyz_.Rho(), posxyz_.Eta(), posxyz_.Phi() ); 
    }


    /// \todo move to PFClusterTools
    static double getDepthCorrection(double energy, bool isBelowPS = false,
				     bool isHadron = false);

    void         setColor(int color) {color_ = color;}

    int          color() const {return color_;}
  
    //C remove this
/*     PFCluster& operator+=(const PFCluster&); */

    PFCluster& operator=(const PFCluster&);

    friend    std::ostream& operator<<(std::ostream& out, 
				       const PFCluster& cluster);
    /// counter
    static unsigned     instanceCounter_;

    /// \todo move to PFClusterTools
    static void setDepthCorParameters(int mode, 
				      double a, double b, 
				      double ap, double bp ) {
      depthCorMode_ = mode;
      depthCorA_ = a; 
      depthCorB_ = b; 
      depthCorAp_ = ap; 
      depthCorBp_ = bp; 
    } 
    

  private:
  
    /// vector of rechit fractions (transient)
    std::vector< reco::PFRecHitFraction >  rechits_;

    /// cluster id
    unsigned      id_;

    /// cluster type
    int           type_;

    /// cluster layer, see PFClusterLayer.h
    int           layer_;          

    /// cluster energy
    double        energy_;

    /// cluster position: cartesian
    math::XYZPoint      posxyz_;

    /// cluster position: rho, eta, phi (transient)
    REPPoint            posrep_;


    /// \todo move to PFClusterTools
    static int    depthCorMode_;
    /// \todo move to PFClusterTools
    static double depthCorA_;
    /// \todo move to PFClusterTools
    static double depthCorB_ ;
    /// \todo move to PFClusterTools
    static double depthCorAp_;
    /// \todo move to PFClusterTools
    static double depthCorBp_;


    /// color (transient)
    int                 color_;
    
    friend class PFClusterAlgo;
  };
}

#endif
