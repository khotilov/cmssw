#ifndef GenericTransientTrackingRecHit_H
#define GenericTransientTrackingRecHit_H

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

class GenericTransientTrackingRecHit: public TransientTrackingRecHit{
public:

//RC   virtual GenericTransientTrackingRecHit * clone() const {
//RC     return new GenericTransientTrackingRecHit(*this);
//RC   }

//RC  virtual RecHitPointer clone( const TrajectoryStateOnSurface&) const {
//RC    return clone();
//RC  }

  virtual ~GenericTransientTrackingRecHit() {delete trackingRecHit_;}

  virtual AlgebraicVector parameters() const {return trackingRecHit_->parameters();}
  virtual AlgebraicSymMatrix parametersError() const {return trackingRecHit_->parametersError();}
  virtual DetId geographicalId() const {return trackingRecHit_->geographicalId();}
  virtual AlgebraicMatrix projectionMatrix() const {return trackingRecHit_->projectionMatrix();}
  virtual int dimension() const {return trackingRecHit_->dimension();}

  virtual LocalPoint localPosition() const {return trackingRecHit_->localPosition();}
  virtual LocalError localPositionError() const {return trackingRecHit_->localPositionError();}

  virtual const GeomDetUnit * detUnit() const;

  virtual bool canImproveWithTrack() const {return false;}

  virtual const TrackingRecHit * hit() const {return trackingRecHit_;};
  
  virtual bool isValid() const{return trackingRecHit_->isValid();}

  virtual std::vector<const TrackingRecHit*> recHits() const {
    return ((const TrackingRecHit *)(trackingRecHit_))->recHits();
  }
  virtual std::vector<TrackingRecHit*> recHits() {
    return trackingRecHit_->recHits();
  }

  static RecHitPointer build( const GeomDet * geom, const TrackingRecHit * rh) {
    return RecHitPointer( new GenericTransientTrackingRecHit( geom, rh));
  }

protected:

  // private constructors enforce usage of builders
  GenericTransientTrackingRecHit(const GeomDet * geom, const TrackingRecHit * rh) :
    TransientTrackingRecHit(geom) {
    trackingRecHit_ = rh->clone();
  }
  GenericTransientTrackingRecHit( const GenericTransientTrackingRecHit & other ) :
    TransientTrackingRecHit( other.det()) {
    trackingRecHit_ = other.hit()->clone();
  }

private:

  TrackingRecHit * trackingRecHit_;

  // should not have assignment operator (?)
  GenericTransientTrackingRecHit & operator= (const GenericTransientTrackingRecHit & t) {
    trackingRecHit_ = t.hit()->clone();
    return *(this);
  }

  // hide the clone method for ReferenceCounted. Warning: this method is still 
  // accessible via the bas class TrackingRecHit interface!
   virtual GenericTransientTrackingRecHit * clone() const {
     return new GenericTransientTrackingRecHit(*this);
   }

};

#endif

