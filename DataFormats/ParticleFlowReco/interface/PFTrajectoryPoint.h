#ifndef DataFormats_ParticleFlowReco_PFTrajectoryPoint_h
#define DataFormats_ParticleFlowReco_PFTrajectoryPoint_h
/** 
 */
#include <vector>
#include <map>
#include <iostream>

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Math/GenVector/PositionVector3D.h"

namespace reco {

  /**\class PFTrajectoryPoint
     \brief Trajectory point (a PFRecTrack holds several trajectory points)
     
     \todo   deal with origin and end vertices of PFSimParticles
     \todo   remove HCAL exit
     \author Renaud Bruneliere
     \date   July 2006
  */
  class PFTrajectoryPoint {

  public:
    typedef ROOT::Math::PositionVector3D<ROOT::Math::CylindricalEta3D<Double32_t> > REPPoint;

    /// Define the different layers where the track can be propagated
    enum LayerType_t {
       /// Point of closest approach from beam axis (initial point in the case of PFSimParticle)
      ClosestApproach = 0,
      BeamPipe = 1,       
      /// Preshower layer 1
      PS1 = 2,             
      /// Preshower layer 2
      PS2 = 3,             
      /// ECAL front face
      ECALEntrance = 4,  
      /// expected maximum of the shower in ECAL, for an e/gamma particle
      /// \todo: add an ECALShowerMaxHadrons
      ECALShowerMax = 5,   
      /// HCAL front face
      HCALEntrance = 6,
      /// HCAL exit
      HCALExit = 7,
      NLayers = 8
    };

    /// default constructor. Set variables at default dummy values
    PFTrajectoryPoint();

    /// constructor from values. 
    /// set detId to -1 if this point is not from a tracker layer
    /// set layer to -1 if this point is not from a tracker layer
    PFTrajectoryPoint(int detId,
		      int layer,
		      const math::XYZPoint& posxyz, 
		      const math::XYZTLorentzVector& momentum); 

    /// copy
    PFTrajectoryPoint(const PFTrajectoryPoint& other);

    /// destructor
    virtual ~PFTrajectoryPoint();


    /// measurement detId
    int detId() const    { return detId_; }

    /// trajectory point layer
    int layer() const    { return layer_; }

    /// is this point valid ? 
    bool     isValid() const {
      if( layer_ == -1 && detId_ == -1 ) return false;
      else return true;
    }

    /// is this point corresponding to an intersection with a tracker layer ?
    bool isTrackerLayer() const {
      if(detId_ >= 0 ) return true; 
      else return false;
    }

    /// cartesian position (x, y, z)
    const math::XYZPoint& positionXYZ() const { return posxyz_; }

    /// trajectory position in (rho, eta, phi) base
    const REPPoint& positionREP() const { return posrep_; }

    /// calculate posrep_ once and for all
    void calculatePositionREP() {
      posrep_.SetCoordinates( posxyz_.Rho(), posxyz_.Eta(), posxyz_.Phi() );
    }

    /// 4-momenta quadrivector
    const math::XYZTLorentzVector& momentum() const    { return momentum_; }

    bool   operator==(const reco::PFTrajectoryPoint& other) const;

    friend std::ostream& operator<<(std::ostream& out, const reco::PFTrajectoryPoint& trajPoint);

  private:

    /// Is the measurement corresponding to a tracker layer or propagated to
    /// a layer external to the tracker (preshower, ECAL, HCAL)
    bool isTrackerLayer_;

    /// detid if measurement is corresponding to a tracker layer
    int detId_;             

    /// propagated layer
    int layer_;

    /// cartesian position (x, y, z)
    math::XYZPoint          posxyz_;

    /// position in (rho, eta, phi) base (transient)
    REPPoint                posrep_;

    /// momentum quadrivector
    math::XYZTLorentzVector momentum_;

  };
  
}

#endif
