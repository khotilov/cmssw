#include "CLHEP/Units/PhysicalConstants.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionForwardEstimator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h" 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "TrackingTools/DetLayers/interface/PhiLess.h"
#include "TrackingTools/DetLayers/interface/rangesIntersect.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalDetRangeZPhi.h"




  // zero value indicates incompatible ts - hit pair
std::pair<bool,double> ConversionForwardEstimator::estimate( const TrajectoryStateOnSurface& ts, 
							const TransientTrackingRecHit& hit) const {
  LogDebug("ConversionForwardEstimator") << "::estimate( const TrajectoryStateOnSurface& ts ...) " << "\n";
  
  std::pair<bool,double> result;
  
  float tsPhi = ts.globalParameters().position().phi();
  GlobalPoint gp = hit.globalPosition();
  float rhPhi = gp.phi();
  float rhR = gp.perp();

  // allow an r fudge of 1.5 * times the sigma
  // nodt used float dr = 1.5 * hit.localPositionError().yy();
//   cout << " err " << hit.globalPositionError().phierr(gp) 
//        << " "     << hit.globalPositionError().rerr(gp) << endl;

  // not used float zLayer = ts.globalParameters().position().z();
  float rLayer = ts.globalParameters().position().perp();

  float newdr = sqrt(pow(dr_,2)+4.*hit.localPositionError().yy());
  float rMin = rLayer - newdr;
  float rMax = rLayer + newdr;
  float phiDiff = tsPhi - rhPhi;
  if (phiDiff > pi) phiDiff -= twopi;
  if (phiDiff < -pi) phiDiff += twopi; 

   LogDebug("ConversionForwardEstimator") << " ForwardEstimator: RecHit at " << gp << "\n";
   LogDebug("ConversionForwardEstimator") << "                   rMin = " << rMin << ", rMax = " << rMax << ", rHit = " << rhR << "\n";
   LogDebug("ConversionForwardEstimator") << "                   thePhiRangeMin = " << thePhiRangeMin << ", thePhiRangeMax = " << thePhiRangeMax << ", phiDiff = " << phiDiff << "\n";

    /*
    if(ptype > 1) {
    cout << "endcap " << gp.perp() <<" " << gp.phi() <<" " << gp.z() <<" " << ptype << " " << phiDiff <<" " << rhR-rLayer << " dr " << dr_ << " " << newdr << endl;
    }
  */
  
  if ( phiDiff < thePhiRangeMax && phiDiff > thePhiRangeMin && 
       rhR < rMax && rhR > rMin) {
  
    /*
    cout << "      estimator returns 1 with phiDiff " << thePhiRangeMin << " < " << phiDiff << " < "
	 << thePhiRangeMax << " and rhR " << rMin << " < " << rhR << " < " << rMax << endl;
    cout << " YES " << phiDiff << " " <<rLayer-rhR << endl;
    cout << "                  => RECHIT ACCEPTED " << endl;
    */
    result.first= true;
    result.second=phiDiff;
  } else {
    /*
    cout << "      estimator returns 0 with phiDiff " << thePhiRangeMin << " < " << phiDiff << " < "
     << thePhiRangeMax << " and  rhR " << rMin << " < " << rhR << " < " << rMax << endl;
    */
    result.first= false;
    result.second=0;    
    
  }
  
  return result;
  
}

bool ConversionForwardEstimator::estimate( const TrajectoryStateOnSurface& ts, 
			   const BoundPlane& plane) const {

   LogDebug("ConversionForwardEstimator") << "::estimate( const TrajectoryStateOnSurface& ts, const BoundPlane& plane) " << "\n";  
  // this method should return one if a detector ring is close enough
  //     to the hit, zero otherwise.
  //     Now time is wasted looking for hits in the rings which are anyhow
  //     too far from the prediction   
  return true ;

}



MeasurementEstimator::Local2DVector
ConversionForwardEstimator::maximalLocalDisplacement( const TrajectoryStateOnSurface& ts,
                                                        const BoundPlane& plane) const
{
  
   LogDebug("ConversionForwardEstimator") << "::maximalLocalDisplacement  " << "\n";
  
  
  if ( ts.hasError() ) {
    LocalError le = ts.localError().positionError();
    LogDebug("ConversionForwardEstimator") << "::maximalLocalDisplacent local error " << le.xx() << " " << le.yy() << "\n";
    //return Local2DVector( sqrt(le.xx())*nSigmaCut(), sqrt(le.yy())*nSigmaCut());
    return Local2DVector( sqrt(le.xx()), sqrt(le.yy()) );
  }
  else return Local2DVector(0,0);
 

}


