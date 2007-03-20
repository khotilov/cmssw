#ifndef RECOTRACKER_TRANSIENTRECHITBUILDER_H
#define RECOTRACKER_TRANSIENTRECHITBUILDER_H

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitMatcher.h"

class TkTransientTrackingRecHitBuilder : public TransientTrackingRecHitBuilder {
  
 public:
  TkTransientTrackingRecHitBuilder (const TrackingGeometry* trackingGeometry, 
				    const PixelClusterParameterEstimator * ,
				    const StripClusterParameterEstimator * ,
                                    const SiStripRecHitMatcher           *);
  TransientTrackingRecHit::RecHitPointer build (const TrackingRecHit * p) const ;
  const PixelClusterParameterEstimator * pixelClusterParameterEstimator(){return pixelCPE;}
  const StripClusterParameterEstimator * stripClusterParameterEstimator(){return stripCPE;}
  const SiStripRecHitMatcher           * siStripRecHitMatcher(){return theMatcher;}
    


 private:
  const TrackingGeometry* tGeometry_;
  const PixelClusterParameterEstimator * pixelCPE;
  const StripClusterParameterEstimator * stripCPE;
  const SiStripRecHitMatcher           * theMatcher;
  
};


#endif
