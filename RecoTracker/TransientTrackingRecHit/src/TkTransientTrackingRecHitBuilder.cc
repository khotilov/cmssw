#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
//
// include all the concrete ones
//
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripMatchedRecHit.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiPixelRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/InvalidTransientRecHit.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DMatchedLocalPos.h"


TkTransientTrackingRecHitBuilder::TkTransientTrackingRecHitBuilder( const TrackingGeometry* trackingGeometry, 
								    const PixelClusterParameterEstimator * p,
								    const StripClusterParameterEstimator * s):  
  tGeometry_(trackingGeometry),
  pixelCPE(p),
  stripCPE(s) {}
  
TransientTrackingRecHit* TkTransientTrackingRecHitBuilder::build (const TrackingRecHit * p) {

  if ( const SiStripRecHit2DLocalPos* sh = dynamic_cast<const SiStripRecHit2DLocalPos*>(p)) { 
    return ( new TSiStripRecHit2DLocalPos(tGeometry_->idToDet(p->geographicalId()), sh, stripCPE ) ); 
  } else if ( const SiStripRecHit2DMatchedLocalPos* mh = dynamic_cast<const SiStripRecHit2DMatchedLocalPos*>(p)) {
    return ( new TSiStripMatchedRecHit(tGeometry_->idToDet(p->geographicalId()), mh )); 
  } else if ( const SiPixelRecHit* ph = dynamic_cast<const SiPixelRecHit*>(p)) {
    return ( new TSiPixelRecHit(tGeometry_->idToDet(p->geographicalId()), ph, pixelCPE) ); 
  }else if (dynamic_cast<const InvalidTrackingRecHit*>(p)){
    return ( new InvalidTransientRecHit((p->geographicalId().rawId() == 0 ? 0 : 
					tGeometry_->idToDet(p->geographicalId())
					) )); 
  }
  return 0;
}
