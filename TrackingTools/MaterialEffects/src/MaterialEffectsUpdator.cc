#include "TrackingTools/TrajectoryState/interface/SurfaceSideDefinition.h"
#include "TrackingTools/MaterialEffects/interface/MaterialEffectsUpdator.h"

using namespace SurfaceSideDefinition;

// static initialization
AlgebraicSymMatrix55  MaterialEffectsUpdator::theNullMatrix;

//
// Update of the trajectory state (implemented in base class since general for
//   all classes returning deltaP and deltaCov.
//
bool MaterialEffectsUpdator::updateStateInPlace (TrajectoryStateOnSurface& TSoS, 
				                 const PropagationDirection propDir) const {
  //
  // Check if 
  // - material is associated to surface
  // - propagation direction is not anyDirection
  // - side of surface is not atCenterOfSurface (could be handled with 50% material?)
  //
  const Surface& surface = TSoS.surface();
  if ( !surface.mediumProperties() || propDir==anyDirection || 
       TSoS.surfaceSide()==atCenterOfSurface )  return true;
  //
  // Check, if already on right side of surface
  //
  if ( (propDir==alongMomentum && TSoS.surfaceSide()==afterSurface ) ||
       (propDir==oppositeToMomentum && TSoS.surfaceSide()==beforeSurface ) )  return true;
  //
  // Update momentum. In case of failure: return invalid state
  //
  LocalTrajectoryParameters lp = TSoS.localParameters();
  if ( !lp.updateP(deltaP(TSoS,propDir)) )  
    return false;
  //
  // Update covariance matrix?
  //
  SurfaceSide side = propDir==alongMomentum ? afterSurface : beforeSurface;
  if ( TSoS.hasError() ) {
    AlgebraicSymMatrix55 eloc = TSoS.localError().matrix() + deltaLocalError(TSoS,propDir);
    //TSoS = TrajectoryStateOnSurface(lp,LocalTrajectoryError(eloc),surface, &(TSoS.globalParameters().magneticField()),side);
    TSoS.update(lp,LocalTrajectoryError(eloc),surface,
                &(TSoS.globalParameters().magneticField()),side);
  }
  else {
    TSoS.update(lp,surface,&(TSoS.globalParameters().magneticField()),side);
    //TSoS = TrajectoryStateOnSurface(lp,surface,&(TSoS.globalParameters().magneticField()),side);
  }
  return true;
}

bool newArguments (const TrajectoryStateOnSurface & TSoS, PropagationDirection  propDir) const {
  // check that track as same momentum and direction, surface has same radLen
  // it optimize also against multiple evaluations on different "surfaces" 
  // belonging to contigous detectors with same radLem 
  bool ok = 
    theLastOverP != TSoS.localParameters().qbp() ||
    theLastDxdz != TSoS.localParameters().dxdz() || 
    theLastRL    != TSoS.surface().mediumProperties()->radLen() ||
    theLastPropDir != propDir;
  if (ok) {
    theLastOverP = TSoS.localParameters().qbp();
    theLastDxdz  = TSoS.localParameters().dxdz(); 
    theLastRL  = TSoS.surface().mediumProperties()->radLen();
    theLastPropDir = propDir;
  }
  return ok;
}
  
