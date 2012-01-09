#ifndef RECOTRACKER_TRANSIENTRACKINGRECHIT_TSiStripRecHit1D_H
#define RECOTRACKER_TRANSIENTRACKINGRECHIT_TSiStripRecHit1D_H

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/HelpertRecHit2DLocalPos.h"
#include "DataFormats/Common/interface/RefGetter.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include<limits>

class TSiStripRecHit1D : public TransientTrackingRecHit{
public:

  typedef SiStripRecHit1D::ClusterRef SiStripClusterRef;
  
  typedef edm::LazyGetter<SiStripCluster>::value_ref  SiStripRegionalClusterRef;

  virtual ~TSiStripRecHit1D() {}

  
  virtual void getKfComponents( KfComponentsHolder & holder ) const {
      HelpertRecHit2DLocalPos().getKfComponents(holder, theHitData, *det()); 
  }
 

  virtual AlgebraicVector parameters() const {return theHitData.parameters();}

  
  virtual AlgebraicSymMatrix parametersError() const {
    return HelpertRecHit2DLocalPos().parError( theHitData.localPositionError(), *det()); 
    //    return theHitData->parametersError();
  }
  

  virtual AlgebraicMatrix projectionMatrix() const {return theHitData.projectionMatrix();}
  virtual int dimension() const {return theHitData.dimension();}

  virtual LocalPoint localPosition() const {return theHitData.localPosition();}
  virtual LocalError localPositionError() const {return theHitData.localPositionError();}

  virtual const TrackingRecHit * hit() const {return &theHitData;};
  
  virtual std::vector<const TrackingRecHit*> recHits() const {
    return hit()->recHits();
  }
  virtual std::vector<TrackingRecHit*> recHits() {
    return theHitData.recHits();
  }

  virtual const GeomDetUnit* detUnit() const;

  virtual bool canImproveWithTrack() const {return true;}

  //RC virtual TSiStripRecHit2DLocalPos* clone(const TrajectoryStateOnSurface& ts) const;
  virtual RecHitPointer clone(const TrajectoryStateOnSurface& ts) const;

  // Extension of the TransientTrackingRecHit interface

  const SiStripRecHit1D* specificHit() const {return &theHitData;};
  const StripClusterParameterEstimator* cpe() const {return theCPE;}

  static RecHitPointer build( const GeomDet * geom, const SiStripRecHit1D* rh,
			      const StripClusterParameterEstimator* cpe,
			      float weight=1., float annealing=1.,
			      bool computeCoarseLocalPosition=false) {
    return RecHitPointer( new TSiStripRecHit1D( geom, rh, cpe, weight, annealing,computeCoarseLocalPosition));
  }

  static RecHitPointer build( const LocalPoint& pos, const LocalError& err,
			      const GeomDet* det,
			      const SiStripClusterRef clust,
			      const StripClusterParameterEstimator* cpe,
			      float weight=1., float annealing=1.) {
    return RecHitPointer( new TSiStripRecHit1D( pos, err, det, clust, cpe, weight, annealing));
  }

  static RecHitPointer build( const LocalPoint& pos, const LocalError& err,
			      const GeomDet* det,
			      const SiStripRegionalClusterRef clust,
			      const StripClusterParameterEstimator* cpe,
			      float weight=1., float annealing=1.) {
    return RecHitPointer( new TSiStripRecHit1D( pos, err, det, clust, cpe, weight, annealing));
  }



private:

  SiStripRecHit1D              theHitData;
  const StripClusterParameterEstimator* theCPE;

  TSiStripRecHit1D (const GeomDet * geom, const SiStripRecHit1D* rh,
		    const StripClusterParameterEstimator* cpe,
		    float weight, float annealing,
		    bool computeCoarseLocalPosition) : 
    TransientTrackingRecHit(geom, weight, annealing), theCPE(cpe) 
    {
      if (rh->hasPositionAndError() || !computeCoarseLocalPosition)
	theHitData = SiStripRecHit1D(*rh);
      else{
      const GeomDetUnit* gdu = dynamic_cast<const GeomDetUnit*>(geom);
      LogDebug("TSiStripRecHit2DLocalPos")<<"calculating coarse position/error.";
      if (gdu){
	if (rh->cluster().isNonnull()){
	  StripClusterParameterEstimator::LocalValues lval= theCPE->localParameters(*rh->cluster(), *gdu);
	  LocalError le(lval.second.xx(),0.,std::numeric_limits<float>::max()); //Correct??
	  theHitData = SiStripRecHit1D(lval.first, le, geom->geographicalId(),rh->cluster());
	}else{
	  StripClusterParameterEstimator::LocalValues lval= theCPE->localParameters(*rh->cluster_regional(), *gdu);
	  LocalError le(lval.second.xx(),0.,std::numeric_limits<float>::max()); //Correct??
	  theHitData = SiStripRecHit1D(lval.first, le, geom->geographicalId(),rh->cluster_regional());
	}
      }else{
	edm::LogError("TSiStripRecHit2DLocalPos")<<" geomdet does not cast into geomdet unit. cannot create strip local parameters.";
	theHitData = SiStripRecHit1D(*rh);
      }
      }
    }

  /// Creates the TrackingRecHit internally, avoids redundent cloning
  TSiStripRecHit1D( const LocalPoint& pos, const LocalError& err,
		    const GeomDet* det,
		    const SiStripClusterRef clust,
		    const StripClusterParameterEstimator* cpe,
		    float weight, float annealing) :
    TransientTrackingRecHit(det, weight, annealing), theHitData(pos, err, det->geographicalId(), clust), 
    theCPE(cpe){} 

  //  TSiStripRecHit2DLocalPos( const TSiStripRecHit2DLocalPos& other ) :
  //  TransientTrackingRecHit( other.det()), 
  //  theHitData( other.specificHit()->clone()),
  //  theCPE( other.cpe()) {}

  TSiStripRecHit1D( const LocalPoint& pos, const LocalError& err,
		    const GeomDet* det,
		    const SiStripRegionalClusterRef clust,			    
		    const StripClusterParameterEstimator* cpe,
		    float weight, float annealing) :
    TransientTrackingRecHit(det, weight, annealing), theHitData(pos, err, det->geographicalId(), clust), 
    theCPE(cpe){} 
  
  

  virtual TSiStripRecHit1D* clone() const {
    return new TSiStripRecHit1D(*this);
  }

};

#endif
