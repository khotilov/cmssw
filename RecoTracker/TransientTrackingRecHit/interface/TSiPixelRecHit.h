
#ifndef RECOTRACKER_TRANSIENTRACKINGRECHIT_TSiPixelRecHit_H
#define RECOTRACKER_TRANSIENTRACKINGRECHIT_TSiPixelRecHit_H

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/HelpertRecHit2DLocalPos.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


// class GeomDetUnit;

class TSiPixelRecHit : public TransientTrackingRecHit {
public:

  typedef SiPixelRecHit::ClusterRef clusterRef;


  virtual ~TSiPixelRecHit() {}

  virtual AlgebraicVector parameters() const {return theHitData.parameters();}

  virtual AlgebraicSymMatrix parametersError() const {
    return HelpertRecHit2DLocalPos().parError( theHitData.localPositionError(), *det()); 
    // return theHitData->parametersError();
  }

  virtual AlgebraicMatrix projectionMatrix() const {return theHitData.projectionMatrix();}
  virtual int dimension() const {return theHitData.dimension();}

  virtual LocalPoint localPosition() const {return theHitData.localPosition();}
  virtual LocalError localPositionError() const {return theHitData.localPositionError();}

  virtual void getKfComponents( KfComponentsHolder & holder ) const {
      HelpertRecHit2DLocalPos().getKfComponents(holder, theHitData, *det()); 
  }

  virtual const TrackingRecHit * hit() const {return &theHitData;};
  
  virtual std::vector<const TrackingRecHit*> recHits() const {
    return hit()->recHits();
  }
  virtual std::vector<TrackingRecHit*> recHits() {
    return theHitData.recHits();
  }

  virtual bool canImproveWithTrack() const {return true;}

  //RC  virtual TSiPixelRecHit* clone (const TrajectoryStateOnSurface& ts) const;
  virtual RecHitPointer clone (const TrajectoryStateOnSurface& ts) const;

  virtual const GeomDetUnit* detUnit() const;

  // Extension of the TransientTrackingRecHit interface

  const SiPixelRecHit* specificHit() const {return &theHitData;};
  const PixelClusterParameterEstimator* cpe() const {return theCPE;}

  static RecHitPointer build( const GeomDet * geom, const SiPixelRecHit* rh, 
			      const PixelClusterParameterEstimator* cpe,
			      float weight=1., float annealing=1.) {
  return RecHitPointer( new TSiPixelRecHit( geom, rh, cpe,weight, annealing));
  }

  static RecHitPointer build( const LocalPoint& pos, const LocalError& err,
			      const GeomDet* det, 
			      clusterRef cluster,
			      const PixelClusterParameterEstimator* cpe,
			      float weight=1., float annealing=1.) {
    return RecHitPointer( new TSiPixelRecHit( pos, err, det, cluster, cpe,weight, annealing));
  }


private:

  SiPixelRecHit                         theHitData;
  const PixelClusterParameterEstimator* theCPE;

  /// This constructor copy the TrackingRecHit, it should be used when the 
  /// TrackingRecHit exist already in some collection
  TSiPixelRecHit(const GeomDet * geom, const SiPixelRecHit* rh, 
		 const PixelClusterParameterEstimator* cpe,
		 float weight, float annealing) : 
    TransientTrackingRecHit(geom, *rh, weight, annealing), theHitData(*rh), theCPE(cpe) {}

  /// Creates the TrackingRecHit internally, avoids redundent cloning
  TSiPixelRecHit( const LocalPoint& pos, const LocalError& err,
				const GeomDet* det, 
		  //				const SiPixelCluster& clust,
		  clusterRef clust,
		  const PixelClusterParameterEstimator* cpe,
		  float weight, float annealing) :
    TransientTrackingRecHit(det,weight, annealing), 
    theHitData( pos, err, det->geographicalId(), clust),
    theCPE(cpe)
  {}
  
  //  TSiPixelRecHit( const TSiPixelRecHit& other ) :
  //  TransientTrackingRecHit( other.det()), 
  //  theHitData( other.specificHit()->clone()),
  //  theCPE( other.cpe())  {}

  virtual TSiPixelRecHit * clone() const {
    return new TSiPixelRecHit(*this);
  }

};



#endif
