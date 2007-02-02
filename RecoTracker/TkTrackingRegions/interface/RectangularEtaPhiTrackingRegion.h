#ifndef RectangularEtaPhiTrackingRegion_H
#define RectangularEtaPhiTrackingRegion_H

/** \class RectangularEtaPhiTrackingRegion
 * A concrete implementation of TrackingRegion. 
 * Apart of vertex constraint from TrackingRegion in this implementation
 * the region of interest is further constrainted in phi and eta 
 * around the direction of the region 
 */

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionBase.h"
#include "RecoTracker/TkTrackingRegions/interface/TkTrackingRegionsMargin.h"
//#include "CommonDet/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "RecoTracker/TkTrackingRegions/interface/HitRZConstraint.h"
#include "RecoTracker/TkTrackingRegions/interface/OuterHitPhiPrediction.h"
#include "FWCore/Framework/interface/EventSetup.h"
class OuterEstimator;
class BarrelDetLayer;
class ForwardDetLayer;

class RectangularEtaPhiTrackingRegion : public TrackingRegionBase {
public:

  typedef TkTrackingRegionsMargin<float> Margin;

 /// dummy constructor
 RectangularEtaPhiTrackingRegion() { }

 /** constructor (symmetric eta and phi margins). <BR>
  * dir        - the direction around which region is constructed <BR>
  *              the initial direction of the momentum of the particle 
  *              should be in the range <BR> 
  *              phi of dir +- deltaPhi <BR>
  *              eta of dir +- deltaEta <BR> 
  *              
  * vertexPos  - the position of the vertex (origin) of the of the region.<BR>
  *              It is a centre of cylinder constraind with rVertex, zVertex.
  *              The track of the particle should cross the cylinder <BR>
  *              WARNING: in the current implementaion the vertexPos is
  *              supposed to be placed on the beam line, i.e. to be of the form
  *              (0,0,float)
  *
  * ptMin      - minimal pt of interest <BR>
  * rVertex    - radius of the cylinder around beam line where the tracks
  *              of interest should point to. <BR>
  * zVertex    - half height of the cylinder around the beam line
  *              where the tracks of interest should point to.   <BR>
  * deltaEta   - allowed deviation of the initial direction of particle
  *              in eta in respect to direction of the region <BR>
  *  deltaPhi  - allowed deviation of the initial direction of particle
  *              in phi in respect to direction of the region 
  */
  RectangularEtaPhiTrackingRegion( const GlobalVector & dir, 
				   const GlobalPoint & vertexPos,
                                   float ptMin, float rVertex, float zVertex,
                                   float deltaEta, float deltaPhi,
                                   float dummy = 0.,
                                   bool precise = true) 
    : TrackingRegionBase( dir, vertexPos, Range( -1/ptMin, 1/ptMin), 
			  rVertex, zVertex),
   thePhiMargin( Margin( fabs(deltaPhi),fabs(deltaPhi))),
   thePrecise(precise)
   { initEtaRange(dir, Margin( fabs(deltaEta),fabs(deltaEta))); }
 
 /** constructor (asymmetrinc eta and phi margins). <BR>
  * non equal left-right eta and phi bounds around direction are
  * possible. The ranges are defined using \c Margin.
  * the meaning of other arguments is the same as in the case of 
  * symmetring bounds to direction of the region.
  */
  RectangularEtaPhiTrackingRegion( const GlobalVector & dir, 
		                       const GlobalPoint & vertexPos,
                                   float ptMin, float rVertex, float zVertex,
                                   Margin etaMargin,
                                   Margin phiMargin,
                                   float dummy = 0., bool precise = true) 
    : TrackingRegionBase( dir, vertexPos, Range( -1/ptMin, 1/ptMin), 
      rVertex, zVertex), thePhiMargin( phiMargin), thePrecise(precise)
    { initEtaRange(dir, etaMargin); }

 /** constructor (explicit pt range, asymmetrinc eta and phi margins). <BR>
  * the meaning of other arguments is the same as in the case of 
  * symmetring bounds to direction of the region.
  */
  RectangularEtaPhiTrackingRegion( const GlobalVector & dir, 
		                       const GlobalPoint & vertexPos,
                                   Range invPtRange, 
                                   float rVertex, float zVertex,
                                   Margin etaMargin,
                                   Margin phiMargin,
                                   float dummy = 0.,
                                   bool precise = true) 
    : TrackingRegionBase( dir, vertexPos, invPtRange, rVertex, zVertex),
      thePhiMargin( phiMargin), thePrecise(precise)
    { initEtaRange(dir, etaMargin); }


  /// allowed eta range [eta_min, eta_max] interval
  const Range & etaRange() const { return theEtaRange; }

  /// defined phi range around phi0, margin is [phi_left,phi_right]. 
  /// region is defined in a range: [phi0-phi_left, phi0+phi_right]
  const Margin & phiMargin() const { return thePhiMargin; }

  /// is precise error calculation switched on 
  bool  isPrecise() const { return thePrecise; }


  //  virtual vector<RecHit> hits(const DetLayer* layer) const;

  virtual HitRZCompatibility* checkRZ(const DetLayer* layer, 
				      const TrackingRecHit*  outerHit,
				      const edm::EventSetup& iSetup) const;

  virtual RectangularEtaPhiTrackingRegion* clone() const { 
    return new RectangularEtaPhiTrackingRegion(*this);
  }

  static std::string const & name() 
    { static std::string local("RectangularEtaPhiTrackingRegion"); return local; }
  virtual std::string const & getName() const {return name();}


private:

  OuterEstimator * estimator(const BarrelDetLayer* layer,const edm::EventSetup& iSetup) const;
  OuterEstimator * estimator(const ForwardDetLayer* layer,const edm::EventSetup& iSetup) const;

  OuterHitPhiPrediction phiWindow(const edm::EventSetup& iSetup) const;
  HitRZConstraint rzConstraint() const;

  void  initEtaRange( const GlobalVector & dir, const Margin& margin);

private:

  Range theEtaRange;
  Margin thePhiMargin;
  bool thePrecise;

};

#endif
