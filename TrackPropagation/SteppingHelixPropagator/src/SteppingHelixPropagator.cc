/** \class SteppingHelixPropagator
 *  Propagator implementation using steps along a helix.
 *  Minimal geometry navigation.
 *  Material effects (multiple scattering and energy loss) are based on tuning
 *  to MC and (eventually) data. 
 *  Implementation file contents follow.
 *
 *  $Date: 2007/05/09 21:48:10 $
 *  $Revision: 1.38 $
 *  \author Vyacheslav Krutelyov (slava77)
 */

//
// Original Author:  Vyacheslav Krutelyov
//         Created:  Fri Mar  3 16:01:24 CST 2006
// $Id: SteppingHelixPropagator.cc,v 1.38 2007/05/09 21:48:10 slava77 Exp $
//
//


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include "MagneticField/VolumeGeometry/interface/MagVolume.h"
#include "Utilities/Timing/interface/TimingReport.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Cone.h"

#include "TrackingTools/GeomPropagators/interface/PropagationExceptions.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <sstream>
#include <typeinfo>

SteppingHelixPropagator::SteppingHelixPropagator() :
  Propagator(anyDirection)
{
  field_ = 0;
}

SteppingHelixPropagator::SteppingHelixPropagator(const MagneticField* field, 
						 PropagationDirection dir):
  Propagator(dir),
  unit66_(AlgebraicMatrixID())
{
  field_ = field;
  vbField_ = dynamic_cast<const VolumeBasedMagneticField*>(field_);
  covRot_ = AlgebraicMatrix66();
  dCTransform_ = unit66_;
  debug_ = false;
  noMaterialMode_ = false;
  noErrorPropagation_ = false;
  applyRadX0Correction_ = false;
  useMagVolumes_ = true;
  useMatVolumes_ = true;
  returnTangentPlane_ = true;
  sendLogWarning_ = false;
  for (int i = 0; i <= MAX_POINTS; i++){
    svBuf_[i].cov = AlgebraicSymMatrix66();
    svBuf_[i].matDCov = AlgebraicSymMatrix66();
    svBuf_[i].isComplete = true;
    svBuf_[i].isValid_ = true;
    svBuf_[i].hasErrorPropagated_ = !noErrorPropagation_;
  }
  defaultStep_ = 5.;
}

TrajectoryStateOnSurface 
SteppingHelixPropagator::propagate(const FreeTrajectoryState& ftsStart, const Plane& pDest) const {
  return propagateWithPath(ftsStart, pDest).first;
}

TrajectoryStateOnSurface 
SteppingHelixPropagator::propagate(const FreeTrajectoryState& ftsStart, const Cylinder& cDest) const
{
  return propagateWithPath(ftsStart, cDest).first;
}

FreeTrajectoryState
SteppingHelixPropagator::propagate(const FreeTrajectoryState& ftsStart, const GlobalPoint& pDest) const
{
  return propagateWithPath(ftsStart, pDest).first;
}

FreeTrajectoryState
SteppingHelixPropagator::propagate(const FreeTrajectoryState& ftsStart, 
				   const GlobalPoint& pDest1, const GlobalPoint& pDest2) const
{
  return propagateWithPath(ftsStart, pDest1, pDest2).first;
}


std::pair<TrajectoryStateOnSurface, double> 
SteppingHelixPropagator::propagateWithPath(const FreeTrajectoryState& ftsStart, 
					   const Plane& pDest) const {

  setIState(SteppingHelixStateInfo(ftsStart));

  const StateInfo& svCurrent = propagate(svBuf_[0], pDest);

  return TsosPP(svCurrent.getStateOnSurface(pDest), svCurrent.path());
}

std::pair<TrajectoryStateOnSurface, double> 
SteppingHelixPropagator::propagateWithPath(const FreeTrajectoryState& ftsStart, 
					   const Cylinder& cDest) const {

  setIState(SteppingHelixStateInfo(ftsStart));

  const StateInfo& svCurrent = propagate(svBuf_[0], cDest);

  return TsosPP(svCurrent.getStateOnSurface(cDest, returnTangentPlane_), svCurrent.path());
}


std::pair<FreeTrajectoryState, double> 
SteppingHelixPropagator::propagateWithPath(const FreeTrajectoryState& ftsStart, 
					   const GlobalPoint& pDest) const {
  setIState(SteppingHelixStateInfo(ftsStart));

  const StateInfo& svCurrent = propagate(svBuf_[0], pDest);

  FreeTrajectoryState ftsDest;
  svCurrent.getFreeState(ftsDest);

  return FtsPP(ftsDest, svCurrent.path());
}

std::pair<FreeTrajectoryState, double> 
SteppingHelixPropagator::propagateWithPath(const FreeTrajectoryState& ftsStart, 
					   const GlobalPoint& pDest1, const GlobalPoint& pDest2) const {

  if ((pDest1-pDest2).mag() < 1e-10){
    if (sendLogWarning_){
      edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: the points should be at a bigger distance"
						<<std::endl;
    }
    return FtsPP();
  }
  setIState(SteppingHelixStateInfo(ftsStart));
  
  const StateInfo& svCurrent = propagate(svBuf_[0], pDest1, pDest2);

  FreeTrajectoryState ftsDest;
  svCurrent.getFreeState(ftsDest);

  return FtsPP(ftsDest, svCurrent.path());
}


const SteppingHelixStateInfo&
SteppingHelixPropagator::propagate(const SteppingHelixStateInfo& sStart, 
				   const Surface& sDest) const {
  
  if (! sStart.isValid()){
    if (sendLogWarning_){
      edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: invalid input state"
						<<std::endl;
    }
    return invalidState_;
  }

  const Plane* pDest = dynamic_cast<const Plane*>(&sDest);
  if (pDest != 0) return propagate(sStart, *pDest);

  const Cylinder* cDest = dynamic_cast<const Cylinder*>(&sDest);
  if (cDest != 0) return propagate(sStart, *cDest);

  throw PropagationException("The surface is neither Cylinder nor Plane");

}

const SteppingHelixStateInfo&
SteppingHelixPropagator::propagate(const SteppingHelixStateInfo& sStart, 
				   const Plane& pDest) const {
  
  if (! sStart.isValid()){
    if (sendLogWarning_){
      edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: invalid input state"
						<<std::endl;
    }    
    return invalidState_;
  }
  setIState(sStart);
  
  GlobalPoint rPlane = pDest.position();
  GlobalVector nPlane(pDest.rotation().zx(), pDest.rotation().zy(), pDest.rotation().zz());

  double pars[6] = { pDest.position().x(), pDest.position().y(), pDest.position().z(),
		     pDest.rotation().zx(), pDest.rotation().zy(), pDest.rotation().zz() };
  
  propagate(PLANE_DT, pars);
  
  return svBuf_[cIndex_(nPoints_-1)];
}

const SteppingHelixStateInfo&
SteppingHelixPropagator::propagate(const SteppingHelixStateInfo& sStart, 
				   const Cylinder& cDest) const {
  
  if (! sStart.isValid()){
    if (sendLogWarning_){
      edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: invalid input state"
						<<std::endl;
    }    
    return invalidState_;
  }
  setIState(sStart);
  
  double pars[6];
  pars[RADIUS_P] = cDest.radius();

  
  propagate(RADIUS_DT, pars);
  
  return svBuf_[cIndex_(nPoints_-1)];
}

const SteppingHelixStateInfo&
SteppingHelixPropagator::propagate(const SteppingHelixStateInfo& sStart, 
				   const GlobalPoint& pDest) const {
  
  if (! sStart.isValid()){
    if (sendLogWarning_){
      edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: invalid input state"
						<<std::endl;
    }    
    return invalidState_;
  }
  setIState(sStart);
  
  double pars[6] = {pDest.x(), pDest.y(), pDest.z(), 0, 0, 0};

  
  propagate(POINT_PCA_DT, pars);
  
  return svBuf_[cIndex_(nPoints_-1)];
}

const SteppingHelixStateInfo&
SteppingHelixPropagator::propagate(const SteppingHelixStateInfo& sStart, 
				   const GlobalPoint& pDest1, const GlobalPoint& pDest2) const {
  
  if ((pDest1-pDest2).mag() < 1e-10 || !sStart.isValid()){
    if (sendLogWarning_){
      if ((pDest1-pDest2).mag() < 1e-10)
	edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: points are too close"
						  <<std::endl;
      if (!sStart.isValid())
	edm::LogWarning("SteppingHelixPropagator")<<"Can't propagate: invalid input state"
						  <<std::endl;
    }
    return invalidState_;
  }
  setIState(sStart);
  
  double pars[6] = {pDest1.x(), pDest1.y(), pDest1.z(),
		    pDest2.x(), pDest2.y(), pDest2.z()};
  
  propagate(LINE_PCA_DT, pars);
  
  return svBuf_[cIndex_(nPoints_-1)];
}

void SteppingHelixPropagator::setIState(const SteppingHelixStateInfo& sStart) const {
  nPoints_ = 0;
  if (sStart.isComplete ) {
    svBuf_[cIndex_(nPoints_)] = sStart;
    nPoints_++;
  } else {
    setIState(sStart.p3, sStart.r3, sStart.q, sStart.cov,
	      propagationDirection());
  }
  svBuf_[cIndex_(0)].hasErrorPropagated_ = sStart.hasErrorPropagated_ & !noErrorPropagation_;
}

void SteppingHelixPropagator::setIState(const SteppingHelixPropagator::Vector& p3, 
					const SteppingHelixPropagator::Point& r3, int charge, 
					const AlgebraicSymMatrix66& cov, PropagationDirection dir) const {
  nPoints_ = 0;
  loadState(svBuf_[cIndex_(nPoints_)], p3, r3, charge, cov, dir);
  nPoints_++;
}

SteppingHelixPropagator::Result 
SteppingHelixPropagator::propagate(SteppingHelixPropagator::DestType type, 
				   const double pars[6], double epsilon)  const{

  static const std::string metname = "SteppingHelixPropagator";
  StateInfo* svCurrent = &svBuf_[cIndex_(nPoints_-1)];

  //check if it's going to work at all
  double tanDist = 0;
  double dist = 0;
  PropagationDirection refDirection = anyDirection;
  Result result = refToDest(type, (*svCurrent), pars, dist, tanDist, refDirection);

  if (result != SteppingHelixStateInfo::OK ){
    svCurrent->status_ = result;
    svCurrent->isValid_ = result == SteppingHelixStateInfo::OK;
    svCurrent->field = field_;
    if (sendLogWarning_){
      edm::LogWarning(metname)<<" Failed after first refToDest check with status "
			      <<SteppingHelixStateInfo::ResultName[result]
			      <<std::endl;
    }
    return result;
  }

  result = SteppingHelixStateInfo::UNDEFINED;
  bool makeNextStep = true;
  double dStep = defaultStep_;
  PropagationDirection dir,oldDir;
  dir = propagationDirection(); 
  oldDir = dir;
  int nOsc = 0;

  double distMag = 1e12;
  double tanDistMag = 1e12;

  double distMat = 1e12;
  double tanDistMat = 1e12;

  double tanDistNextCheck = -0.1;//just need a negative start val
  double tanDistMagNextCheck = -0.1;
  double tanDistMatNextCheck = -0.1;
  double oldDStep = 0;
  PropagationDirection oldRefDirection = propagationDirection();

  Result resultToMat = SteppingHelixStateInfo::UNDEFINED;
  Result resultToMag = SteppingHelixStateInfo::UNDEFINED;

  bool isFirstStep = true;

  while (makeNextStep){
    dStep = defaultStep_;
    svCurrent = &svBuf_[cIndex_(nPoints_-1)];
    double curZ = svCurrent->r3.z();
    double curR = svCurrent->r3.perp();
    if ( fabs(curZ) < 440 && curR < 260) dStep = defaultStep_*2;
    //    refDirection = propagationDirection();

    tanDistNextCheck -= oldDStep;
    tanDistMagNextCheck -= oldDStep;
    tanDistMatNextCheck -= oldDStep;
    
    if (tanDistNextCheck < 0){
      //use pre-computed values if it's the first step
      if (! isFirstStep) refToDest(type, (*svCurrent), pars, dist, tanDist, refDirection);
      tanDistNextCheck = fabs(tanDist)*0.5 - 0.5; //need a better guess (to-do)
      //reasonable limit
      if (tanDistNextCheck >  defaultStep_*20. ) tanDistNextCheck = defaultStep_*20.;
      oldRefDirection = refDirection;
    } else {
      tanDist  = tanDist > 0. ? tanDist - oldDStep : tanDist + oldDStep; 
      refDirection = oldRefDirection;
      if (debug_) LogTrace(metname)<<"Skipped refToDest: guess tanDist = "<<tanDist
				   <<" next check at "<<tanDistNextCheck<<std::endl;
    }

    if (propagationDirection() == anyDirection){
      dir = refDirection;
    } else {
      dir = propagationDirection();
    }

    if (useMagVolumes_ && ! (fabs(dist) < fabs(epsilon))){//need to know the general direction
      if (tanDistMagNextCheck < 0){
	resultToMag = refToMagVolume((*svCurrent), dir, distMag, tanDistMag, fabs(dist));
	tanDistMagNextCheck = fabs(tanDistMag)*0.5-0.5; //need a better guess (to-do)
	//reasonable limit; "turn off" checking if bounds are further than the destination
	if (tanDistMagNextCheck >  defaultStep_*20. 
	    || fabs(dist) < fabs(distMag)
	    || resultToMag ==SteppingHelixStateInfo::INACC) tanDistMagNextCheck = defaultStep_*20.;	
	if (resultToMag != SteppingHelixStateInfo::INACC 
	    && resultToMag != SteppingHelixStateInfo::OK) tanDistMagNextCheck = -1;
      } else {
	//	resultToMag = SteppingHelixStateInfo::OK;
	tanDistMag  = tanDistMag > 0. ? tanDistMag - oldDStep : tanDistMag + oldDStep; 
	if (debug_) LogTrace(metname)<<"Skipped refToMag: guess tanDistMag = "<<tanDistMag<<std::endl;
      }
    }

    if (useMatVolumes_ && ! (fabs(dist) < fabs(epsilon))){//need to know the general direction
      if (tanDistMatNextCheck < 0){
	resultToMat = refToMatVolume((*svCurrent), dir, distMat, tanDistMat, fabs(dist));
	tanDistMatNextCheck = fabs(tanDistMat)*0.5-0.5; //need a better guess (to-do)
	//reasonable limit; "turn off" checking if bounds are further than the destination
	if (tanDistMatNextCheck >  defaultStep_*20. 
	    || fabs(dist) < fabs(distMat)
	    || resultToMat ==SteppingHelixStateInfo::INACC ) tanDistMatNextCheck = defaultStep_*20.;
	if (resultToMat != SteppingHelixStateInfo::INACC 
	    && resultToMat != SteppingHelixStateInfo::OK) tanDistMatNextCheck = -1;
      } else {
	//	resultToMat = SteppingHelixStateInfo::OK;
	tanDistMat  = tanDistMat > 0. ? tanDistMat - oldDStep : tanDistMat + oldDStep; 
	if (debug_) LogTrace(metname)<<"Skipped refToMat: guess tanDistMat = "<<tanDistMat<<std::endl;
      }
    }

    double rDotP = svCurrent->r3.dot(svCurrent->p3);
    if ((fabs(curZ) > 1.5e3 || curR >800.) 
	&& ((dir == alongMomentum && rDotP > 0) 
	    || (dir == oppositeToMomentum && rDotP < 0) )
	){
      dStep = fabs(tanDist) -1e-12;
    }
    double tanDistMin = fabs(tanDist);
    if (tanDistMin > fabs(tanDistMag)+0.05 && resultToMag == SteppingHelixStateInfo::OK){
      tanDistMin = fabs(tanDistMag)+0.05;     //try to step into the next volume
    }
    if (tanDistMin > fabs(tanDistMat)+0.05 && resultToMat == SteppingHelixStateInfo::OK){
      tanDistMin = fabs(tanDistMat)+0.05;     //try to step into the next volume
    }
    if (fabs(tanDistMin) < dStep){
      dStep = fabs(tanDistMin); 
      if (type == POINT_PCA_DT || type == LINE_PCA_DT){
	//being lazy here; the best is to take into account the curvature
	dStep = fabs(tanDistMin)*0.5; 
      }
    }
    //keep this path length for the next step
    oldDStep = dStep;

    if (dStep > 1e-10 && ! (fabs(dist) < fabs(epsilon))){
      StateInfo* svNext = &svBuf_[cIndex_(nPoints_)];
      makeAtomStep((*svCurrent), (*svNext), dStep, dir, HEL_AS_F);
      nPoints_++;    svCurrent = &svBuf_[cIndex_(nPoints_-1)];
      if (oldDir != dir){
	nOsc++;
	tanDistNextCheck = -1;//check dist after osc
	tanDistMagNextCheck = -1;
	tanDistMatNextCheck = -1;
      }
      oldDir = dir;
    }

    if (nOsc>1 && fabs(dStep)>epsilon){
      if (debug_) LogTrace(metname)<<"Ooops"<<std::endl;
    }

    if (fabs(dist) < fabs(epsilon)  ) result = SteppingHelixStateInfo::OK;

    if ((type == POINT_PCA_DT || type == LINE_PCA_DT )
	&& fabs(dStep) < fabs(epsilon)  ){
      //now check if it's not a branch point (peek ahead at 1 cm)
      double nextDist = 0;
      double nextTanDist = 0;
      PropagationDirection nextRefDirection = anyDirection;
      StateInfo* svNext = &svBuf_[cIndex_(nPoints_)];
      makeAtomStep((*svCurrent), (*svNext), 1., dir, HEL_AS_F);
      nPoints_++;     svCurrent = &svBuf_[cIndex_(nPoints_-1)];
      refToDest(type, (*svCurrent), pars, nextDist, nextTanDist, nextRefDirection);
      if ( fabs(nextDist) > fabs(dist)){
	nPoints_--;      svCurrent = &svBuf_[cIndex_(nPoints_-1)];
	result = SteppingHelixStateInfo::OK;
	if (debug_){
	  LogTrace(metname)<<"Found real local minimum in PCA"<<std::endl;
	}
      } else {
	//keep this trial point and continue
	dStep = defaultStep_;
	if (debug_){
	  LogTrace(metname)<<"Found branch point in PCA"<<std::endl;
	}
      }
    }

    if (nPoints_ > MAX_STEPS*1./defaultStep_ || nOsc > 6) result = SteppingHelixStateInfo::FAULT;

    if (svCurrent->p3.mag() < 0.1 ) result = SteppingHelixStateInfo::RANGEOUT;

    if ( curR > 20000 || fabs(curZ) > 20000 ) result = SteppingHelixStateInfo::INACC;

    makeNextStep = result == SteppingHelixStateInfo::UNDEFINED;
    svCurrent->status_ = result;
    svCurrent->isValid_ = result == SteppingHelixStateInfo::OK;
    svCurrent->field = field_;

    isFirstStep = false;
  }

  if (sendLogWarning_ && result != SteppingHelixStateInfo::OK){
    edm::LogWarning(metname)<<" Propagation failed with status "
			    <<SteppingHelixStateInfo::ResultName[result]
			    <<std::endl;
    if (result == SteppingHelixStateInfo::RANGEOUT)
      edm::LogWarning(metname)<<" Momentum at last point is too low (<0.1) p_last = "
			      <<svCurrent->p3.mag()
			      <<std::endl;
    if (result == SteppingHelixStateInfo::INACC)
      edm::LogWarning(metname)<<" Went too far: the last point is at "<<svCurrent->r3
			      <<std::endl;
    if (result == SteppingHelixStateInfo::FAULT && nOsc > 6)
      edm::LogWarning(metname)<<" Infinite loop condidtion detected: going in cycles. Break after 6 cycles"
			      <<std::endl;
    if (result == SteppingHelixStateInfo::FAULT && nPoints_ > MAX_STEPS*1./defaultStep_)
      edm::LogWarning(metname)<<" Tired to go farther. Made too many steps: more than "
			      <<MAX_STEPS*1./defaultStep_
			      <<std::endl;
    
  }

  if (debug_){
    switch (type) {
    case RADIUS_DT:
      LogTrace(metname)<<"going to radius "<<pars[RADIUS_P]<<std::endl;
      break;
    case Z_DT:
      LogTrace(metname)<<"going to z "<<pars[Z_P]<<std::endl;
      break;
    case PATHL_DT:
      LogTrace(metname)<<"going to pathL "<<pars[PATHL_P]<<std::endl;
      break;
    case PLANE_DT:
      {
	Point rPlane(pars[0], pars[1], pars[2]);
	Vector nPlane(pars[3], pars[4], pars[5]);
	LogTrace(metname)<<"going to plane r0:"<<rPlane<<" n:"<<nPlane<<std::endl;
      }
      break;
    case POINT_PCA_DT:
      {
	Point rDest(pars[0], pars[1], pars[2]);
	LogTrace(metname)<<"going to PCA to point "<<rDest<<std::endl;
      }
      break;
    case LINE_PCA_DT:
      {
	Point rDest1(pars[0], pars[1], pars[2]);
	Point rDest2(pars[3], pars[4], pars[5]);
	LogTrace(metname)<<"going to PCA to line "<<rDest1<<" - "<<rDest2<<std::endl;
      }
      break;
    default:
      LogTrace(metname)<<"going to NOT IMPLEMENTED"<<std::endl;
      break;
    }
    LogTrace(metname)<<"Made "<<nPoints_-1<<" steps and stopped at(cur step) "<<svCurrent->r3<<" nOsc "<<nOsc<<std::endl;
  }
  
  return result;
}
  
void SteppingHelixPropagator::loadState(SteppingHelixPropagator::StateInfo& svCurrent, 
					const SteppingHelixPropagator::Vector& p3, 
					const SteppingHelixPropagator::Point& r3, int charge,
					const AlgebraicSymMatrix66& cov, PropagationDirection dir) const{
  static const std::string metname = "SteppingHelixPropagator";
  svCurrent.q = charge;
  svCurrent.p3 = p3;
  svCurrent.r3 = r3;
  svCurrent.dir = dir == alongMomentum ? 1.: -1.;

  svCurrent.path_ = 0; // this could've held the initial path
  svCurrent.radPath = 0;

  GlobalPoint gPoint(r3.x(), r3.y(), r3.z());
  GlobalPoint gPointNegZ(svCurrent.r3.x(), svCurrent.r3.y(), -fabs(svCurrent.r3.z()));

  GlobalVector bf;
  // = field_->inTesla(gPoint);
  if (useMagVolumes_){
    if (vbField_ ){
      svCurrent.magVol = vbField_->findVolume(gPointNegZ);
    } else {
      edm::LogWarning(metname)<<"Failed to cast into VolumeBasedMagneticField: fall back to the default behavior"<<std::endl;
      svCurrent.magVol = 0;
    }
    if (debug_){
      LogTrace(metname)<<"Got volume at "<<svCurrent.magVol<<std::endl;
    }
  }
  
  if (useMagVolumes_ && svCurrent.magVol != 0){
    bf = svCurrent.magVol->inTesla(gPointNegZ);
    if (gPoint.z() > 0){
      svCurrent.bf.set(-bf.x(), -bf.y(), bf.z());
    } else {
      svCurrent.bf.set(bf.x(), bf.y(), bf.z());
    }
  } else {
    bf = field_->inTesla(gPoint);
    svCurrent.bf.set(bf.x(), bf.y(), bf.z());
  }
  if (svCurrent.bf.mag() < 1e-6) svCurrent.bf.set(0., 0., 1e-6);



  double dEdXPrime = 0;
  double dEdx = 0;
  double radX0 = 0;
  MatBounds rzLims;
  dEdx = getDeDx(svCurrent, dEdXPrime, radX0, rzLims);
  svCurrent.dEdx = dEdx;    svCurrent.dEdXPrime = dEdXPrime;
  svCurrent.radX0 = radX0;
  svCurrent.rzLims = rzLims;

  svCurrent.cov =cov;

  svCurrent.isComplete = true;

  if (debug_){
    LogTrace(metname)<<"Loaded at  path: "<<svCurrent.path_<<" radPath: "<<svCurrent.radPath
		     <<" p3 "<<" pt: "<<svCurrent.p3.perp()<<" phi: "<<svCurrent.p3.phi()
		     <<" eta: "<<svCurrent.p3.eta()
		     <<" "<<svCurrent.p3
		     <<" r3: "<<svCurrent.r3
		     <<" bField: "<<svCurrent.bf.mag()
		     <<std::endl;
    LogTrace(metname)<<"Input Covariance in Global RF "<<cov<<std::endl;
  }
}

void SteppingHelixPropagator::getNextState(const SteppingHelixPropagator::StateInfo& svPrevious, 
					   SteppingHelixPropagator::StateInfo& svNext,
					   double dP, const SteppingHelixPropagator::Vector& tau,
					   const SteppingHelixPropagator::Vector& drVec, double dS, double dX0,
					   const AlgebraicMatrix66& dCovTransform) const{
  static const std::string metname = "SteppingHelixPropagator";
  svNext.q = svPrevious.q;
  svNext.dir = dS > 0.0 ? 1.: -1.; 
  svNext.p3 = tau;  svNext.p3*=(svPrevious.p3.mag() - svNext.dir*fabs(dP));

  svNext.r3 = svPrevious.r3; svNext.r3 += drVec;

  svNext.path_ = svPrevious.path_ + dS;
  svNext.radPath = svPrevious.radPath + dX0;

  GlobalPoint gPointNegZ(svNext.r3.x(), svNext.r3.y(), -fabs(svNext.r3.z()));

  GlobalVector bf; 

  if (useMagVolumes_){
    if (vbField_ != 0){
      svNext.magVol = vbField_->findVolume(gPointNegZ);
    } else {
      LogTrace(metname)<<"Failed to cast into VolumeBasedMagneticField"<<std::endl;
      svNext.magVol = 0;
    }
    if (debug_){
      LogTrace(metname)<<"Got volume at "<<svNext.magVol<<std::endl;
    }
  }

  if (useMagVolumes_ && svNext.magVol != 0){
    bf = svNext.magVol->inTesla(gPointNegZ);
    if (svNext.r3.z() > 0){
      svNext.bf.set(-bf.x(), -bf.y(), bf.z());
    } else {
      svNext.bf.set(bf.x(), bf.y(), bf.z());
    }
  } else {
    GlobalPoint gPoint(svNext.r3.x(), svNext.r3.y(), svNext.r3.z());
    bf = field_->inTesla(gPoint);
    svNext.bf.set(bf.x(), bf.y(), bf.z());
  }
  if (svNext.bf.mag() < 1e-6) svNext.bf.set(0., 0., 1e-6);
  
  
  double dEdXPrime = 0;
  double dEdx = 0;
  double radX0 = 0;
  MatBounds rzLims;
  dEdx = getDeDx(svNext, dEdXPrime, radX0, rzLims);
  svNext.dEdx = dEdx;    svNext.dEdXPrime = dEdXPrime;
  svNext.radX0 = radX0;
  svNext.rzLims = rzLims;

  //update Emat only if it's valid
  svNext.hasErrorPropagated_ = svPrevious.hasErrorPropagated_;
  if (svPrevious.hasErrorPropagated_){
    //    svNext.cov = ROOT::Math::Similarity(dCovTransform, svPrevious.cov);
    AlgebraicMatrix66 tmp = dCovTransform*svPrevious.cov;
    ROOT::Math::AssignSym::Evaluate(svNext.cov, tmp*ROOT::Math::Transpose(dCovTransform));

    svNext.cov += svPrevious.matDCov;
  } else {
    //could skip dragging along the unprop. cov later .. now
    // svNext.cov = svPrevious.cov;
  }

  if (debug_){
    LogTrace(metname)<<"Now at  path: "<<svNext.path_<<" radPath: "<<svNext.radPath
		     <<" p3 "<<" pt: "<<svNext.p3.perp()<<" phi: "<<svNext.p3.phi()
		     <<" eta: "<<svNext.p3.eta()
		     <<" "<<svNext.p3
		     <<" r3: "<<svNext.r3
		     <<" dPhi: "<<acos(svNext.p3.unit().dot(svPrevious.p3.unit()))
		     <<" bField: "<<svNext.bf.mag()
		     <<" dedx: "<<svNext.dEdx
		     <<std::endl;
    LogTrace(metname)<<"New Covariance "<<svNext.cov<<std::endl;
    LogTrace(metname)<<"Transf by dCovTransform "<<dCovTransform<<std::endl;
  }
}

void SteppingHelixPropagator::setRep(SteppingHelixPropagator::Basis& rep, 
				     const SteppingHelixPropagator::Vector& tau) const{
  Vector zRep(0., 0., 1.);
  rep.lX = tau;
  rep.lY = zRep.cross(tau); rep.lY *= 1./tau.perp();
  rep.lZ = rep.lX.cross(rep.lY);
}

bool SteppingHelixPropagator::makeAtomStep(SteppingHelixPropagator::StateInfo& svCurrent,
					   SteppingHelixPropagator::StateInfo& svNext,
					   double dS, 
					   PropagationDirection dir, 
					   SteppingHelixPropagator::Fancy fancy) const{
  static const std::string metname = "SteppingHelixPropagator";
  if (debug_){
    LogTrace(metname)<<"Make atom step "<<svCurrent.path_<<" with step "<<dS<<" in direction "<<dir<<std::endl;
  }

  double dP = 0;
  Vector tau = svCurrent.p3; tau *= 1./tau.mag();
  Vector tauNext(tau);
  Vector drVec;

  dS = dir == alongMomentum ? fabs(dS) : -fabs(dS);


  double radX0 = 1e24;

  switch (fancy){
  case HEL_AS_F:
  case HEL_ALL_F:{
    double p0 = svCurrent.p3.mag();
    double b0 = svCurrent.bf.mag();

    //get to the mid-point first
    double phi = 0.0029979*svCurrent.q*b0/p0*dS/2.;
    bool phiSmall = fabs(phi) < 3e-8;

    double cosPhi = cos(phi);
    double sinPhi = sin(phi);

    double oneLessCosPhi = 1.-cosPhi;
    double oneLessCosPhiOPhi = oneLessCosPhi/phi;
    double sinPhiOPhi = sinPhi/phi;
    double phiLessSinPhiOPhi = 1 - sinPhiOPhi;
    if (phiSmall){
      oneLessCosPhi = 0.5*phi*phi;//*(1.- phi*phi/12.);
      oneLessCosPhiOPhi = 0.5*phi;//*(1.- phi*phi/12.);
      sinPhiOPhi = 1. - phi*phi/6.;
      phiLessSinPhiOPhi = phi*phi/6.;//*(1. - phi*phi/20.);
    }

    Vector bHat = svCurrent.bf; bHat *= 1./bHat.mag();
    Vector btVec(bHat.cross(tau));
    Vector bbtVec(bHat.cross(btVec));

    //don't need it here    tauNext = tau + bbtVec*oneLessCosPhi - btVec*sinPhi;
    drVec = tau; drVec += bbtVec*phiLessSinPhiOPhi; drVec -= btVec*oneLessCosPhiOPhi;
    drVec *= dS/2.;

    double dEdx = svCurrent.dEdx;
    double dEdXPrime = svCurrent.dEdXPrime;
    radX0 = svCurrent.radX0;
    dP = dEdx*dS;

    //improve with above values:
    drVec += svCurrent.r3;
    GlobalVector bfGV;
    Vector bf; //(bfGV.x(), bfGV.y(), bfGV.z());
    // = svCurrent.magVol->inTesla(GlobalPoint(drVec.x(), drVec.y(), -fabs(drVec.z())));
    if (useMagVolumes_ && svCurrent.magVol != 0){
      // this negative-z business will break at some point
      bfGV = svCurrent.magVol->inTesla(GlobalPoint(drVec.x(), drVec.y(), -fabs(drVec.z())));
      if (drVec.z() > 0){
	bf.set(-bfGV.x(), -bfGV.y(), bfGV.z());
      } else {
	bf.set(bfGV.x(), bfGV.y(), bfGV.z());
      }
    } else {
      bfGV = field_->inTesla(GlobalPoint(drVec.x(), drVec.y(), drVec.z()));
      bf.set(bfGV.x(), bfGV.y(), bfGV.z());
    }
    b0 = bf.mag();
    if (b0 < 1e-6) {
      b0 = 1e-6;
      bf.set(0., 0., 1e-6);
    }
    if (debug_){
      LogTrace(metname)<<"Improved b "<<b0
		       <<" at r3 "<<drVec<<std::endl;
    }

    if (fabs((b0-svCurrent.bf.mag())*dS) > 1){
      //missed the mag volume boundary?
      if (debug_){
	LogTrace(metname)<<"Large bf*dS change "<<fabs((b0-svCurrent.bf.mag())*dS)
			 <<" --> recalc dedx"<<std::endl;
      }
      svNext.r3 = drVec;
      svNext.bf = bf;
      svNext.p3 = svCurrent.p3;
      MatBounds rzTmp;
      dEdx = getDeDx(svNext, dEdXPrime, radX0, rzTmp);
      dP = dEdx*dS;      
    }
    //p0 is mid-way and b0 from mid-point
    p0 += dP/2.; p0 = p0 < 1e-2 ? 1e-2 : p0;

    phi = 0.0029979*svCurrent.q*b0/p0*dS;
    phiSmall = fabs(phi) < 3e-8;

    if (phiSmall){
      sinPhi = phi;
      cosPhi = 1. -phi*phi/2;
      oneLessCosPhi = 0.5*phi*phi;//*(1.- phi*phi/12.); //<-- ~below double-precision for phi<3e-8
      oneLessCosPhiOPhi = 0.5*phi;//*(1.- phi*phi/12.);
      sinPhiOPhi = 1. - phi*phi/6.;
      phiLessSinPhiOPhi = phi*phi/6.;//*(1. - phi*phi/20.);
    }else {
      cosPhi = cos(phi); 
      sinPhi = sin(phi);
      oneLessCosPhi = 1.-cosPhi;
      oneLessCosPhiOPhi = oneLessCosPhi/phi;
      sinPhiOPhi = sinPhi/phi;
      phiLessSinPhiOPhi = 1. - sinPhiOPhi;
    }

    bHat = bf; bHat *= 1./bHat.mag();
    btVec = bHat.cross(tau);
    bbtVec = bHat.cross(btVec);

    tauNext = tau; tauNext += bbtVec*oneLessCosPhi; tauNext -= btVec*sinPhi;
    drVec = tau; drVec += bbtVec*phiLessSinPhiOPhi; drVec -= btVec*oneLessCosPhiOPhi;
    drVec *= dS;
    
    
    if (svCurrent.hasErrorPropagated_){
      double theta02 = 14.e-3/p0*sqrt(fabs(dS)/radX0); // .. drop log term (this is non-additive)
      theta02 *=theta02;
      if (applyRadX0Correction_){
	// this provides the integrand for theta^2
	// if summed up along the path, should result in 
	// theta_total^2 = Int_0^x0{ f(x)dX} = (13.6/p0)^2*x0*(1+0.036*ln(x0+1))
	// x0+1 above is to make the result infrared safe.
	double x0 = fabs(svCurrent.radPath);
	double dX0 = fabs(dS)/radX0;
	double alphaX0 = 13.6e-3/p0; alphaX0 *= alphaX0;
	double betaX0 = 0.038;
	theta02 = dX0*alphaX0*(1+betaX0*log(x0+1))*(1 + betaX0*log(x0+1) + 2.*betaX0*x0/(x0+1) );
      }
      
      double epsilonP0 = 1.+ dP/p0;
      double omegaP0 = -dP/p0 + dS*dEdXPrime;      
      

      double dsp = dS/p0;

      Vector tbtVec(tau.cross(btVec));

      dCTransform_ = unit66_;
      //make everything in global
      //case I: no "spatial" derivatives |--> dCtr({1,2,3,4,5,6}{1,2,3}) = 0    
      dCTransform_(0,3) = dsp*(bHat.x()*tbtVec.x() 
			       + cosPhi*tau.x()*bbtVec.x()
			       + ((1.-bHat.x()*bHat.x()) + phi*tau.x()*btVec.x())*sinPhiOPhi);

      dCTransform_(0,4) = dsp*(bHat.z()*oneLessCosPhiOPhi + bHat.x()*tbtVec.y()
			       + cosPhi*tau.y()*bbtVec.x() 
			       + (-bHat.x()*bHat.y() + phi*tau.y()*btVec.x())*sinPhiOPhi);
      dCTransform_(0,5) = dsp*(-bHat.y()*oneLessCosPhiOPhi + bHat.x()*tbtVec.z()
			       + cosPhi*tau.z()*bbtVec.x()
			       + (-bHat.x()*bHat.z() + phi*tau.z()*btVec.x())*sinPhiOPhi);

      dCTransform_(1,3) = dsp*(-bHat.z()*oneLessCosPhiOPhi + bHat.y()*tbtVec.x()
			       + cosPhi*tau.x()*bbtVec.y()
			       + (-bHat.x()*bHat.y() + phi*tau.x()*btVec.y())*sinPhiOPhi);
      dCTransform_(1,4) = dsp*(bHat.y()*tbtVec.y() 
			       + cosPhi*tau.y()*bbtVec.y()
			       + ((1.-bHat.y()*bHat.y()) + phi*tau.y()*btVec.y())*sinPhiOPhi);
      dCTransform_(1,5) = dsp*(bHat.x()*oneLessCosPhiOPhi + bHat.y()*tbtVec.z()
			       + cosPhi*tau.z()*bbtVec.y() 
			       + (-bHat.y()*bHat.z() + phi*tau.z()*btVec.y())*sinPhiOPhi);

      dCTransform_(2,3) = dsp*(bHat.y()*oneLessCosPhiOPhi + bHat.z()*tbtVec.x()
			       + cosPhi*tau.x()*bbtVec.z() 
			       + (-bHat.x()*bHat.z() + phi*tau.x()*btVec.z())*sinPhiOPhi);
      dCTransform_(2,4) = dsp*(-bHat.x()*oneLessCosPhiOPhi + bHat.z()*tbtVec.y()
			       + cosPhi*tau.y()*bbtVec.z()
			       + (-bHat.y()*bHat.z() + phi*tau.y()*btVec.z())*sinPhiOPhi);
      dCTransform_(2,5) = dsp*(bHat.z()*tbtVec.z() 
			       + cosPhi*tau.z()*bbtVec.z()
			       + ((1.-bHat.z()*bHat.z()) + phi*tau.z()*btVec.z())*sinPhiOPhi);


      dCTransform_(3,3) = epsilonP0*(1. - oneLessCosPhi*(1.-bHat.x()*bHat.x())
				     + phi*tau.x()*(cosPhi*btVec.x() - sinPhi*bbtVec.x()))
	+ omegaP0*tau.x()*tauNext.x();
      dCTransform_(3,4) = epsilonP0*(bHat.x()*bHat.y()*oneLessCosPhi + bHat.z()*sinPhi
				     + phi*tau.y()*(cosPhi*btVec.x() - sinPhi*bbtVec.x()))
	+ omegaP0*tau.y()*tauNext.x();
      dCTransform_(3,5) = epsilonP0*(bHat.x()*bHat.z()*oneLessCosPhi - bHat.y()*sinPhi
				     + phi*tau.z()*(cosPhi*btVec.x() - sinPhi*bbtVec.x()))
	+ omegaP0*tau.z()*tauNext.x();

      dCTransform_(4,3) = epsilonP0*(bHat.x()*bHat.y()*oneLessCosPhi - bHat.z()*sinPhi
				     + phi*tau.x()*(cosPhi*btVec.y() - sinPhi*bbtVec.y()))
	+ omegaP0*tau.x()*tauNext.y();
      dCTransform_(4,4) = epsilonP0*(1. - oneLessCosPhi*(1.-bHat.y()*bHat.y())
				     + phi*tau.y()*(cosPhi*btVec.y() - sinPhi*bbtVec.y()))
	+ omegaP0*tau.y()*tauNext.y();
      dCTransform_(4,5) = epsilonP0*(bHat.y()*bHat.z()*oneLessCosPhi + bHat.x()*sinPhi
				     + phi*tau.z()*(cosPhi*btVec.y() - sinPhi*bbtVec.y()))
	+ omegaP0*tau.z()*tauNext.y();
    
      dCTransform_(5,3) = epsilonP0*(bHat.x()*bHat.z()*oneLessCosPhi + bHat.y()*sinPhi
				     + phi*tau.x()*(cosPhi*btVec.z() - sinPhi*bbtVec.z()))
	+ omegaP0*tau.x()*tauNext.z();
      dCTransform_(5,4) = epsilonP0*(bHat.y()*bHat.z()*oneLessCosPhi - bHat.x()*sinPhi
				     + phi*tau.y()*(cosPhi*btVec.z() - sinPhi*bbtVec.z()))
	+ omegaP0*tau.y()*tauNext.z();
      dCTransform_(5,5) = epsilonP0*(1. - oneLessCosPhi*(1.-bHat.z()*bHat.z())
				     + phi*tau.z()*(cosPhi*btVec.z() - sinPhi*bbtVec.z()))
	+ omegaP0*tau.z()*tauNext.z();
    

      Basis rep; setRep(rep, tauNext);
      //mind the sign of dS and dP (dS*dP < 0 allways)
      //covariance should grow no matter which direction you propagate
      //==> take abs values.
      //reset not needed: fill all below  svCurrent.matDCov *= 0.;
      double mulRR = theta02*dS*dS/3.;
      double mulRP = theta02*fabs(dS)*p0/2.;
      double mulPP = theta02*p0*p0;
      double losPP = dP*dP*1.6/fabs(dS)*(1.0 + p0*1e-3);
      //another guess .. makes sense for 1 cm steps 2./dS == 2 [cm] / dS [cm] at low pt
      //double it by 1TeV
      //not gaussian anyways
      // derived from the fact that sigma_p/eLoss ~ 0.08 after ~ 200 steps

      //symmetric RR part
      svCurrent.matDCov(0,0) = mulRR*(rep.lY.x()*rep.lY.x() + rep.lZ.x()*rep.lZ.x());
      svCurrent.matDCov(0,1) = mulRR*(rep.lY.x()*rep.lY.y() + rep.lZ.x()*rep.lZ.y());
      svCurrent.matDCov(0,2) = mulRR*(rep.lY.x()*rep.lY.z() + rep.lZ.x()*rep.lZ.z());
      svCurrent.matDCov(1,1) = mulRR*(rep.lY.y()*rep.lY.y() + rep.lZ.y()*rep.lZ.y());
      svCurrent.matDCov(1,2) = mulRR*(rep.lY.y()*rep.lY.z() + rep.lZ.y()*rep.lZ.z());
      svCurrent.matDCov(2,2) = mulRR*(rep.lY.z()*rep.lY.z() + rep.lZ.z()*rep.lZ.z());

      //symmetric PP part
      svCurrent.matDCov(3,3) = mulPP*(rep.lY.x()*rep.lY.x() + rep.lZ.x()*rep.lZ.x())
	+ losPP*rep.lX.x()*rep.lX.x();
      svCurrent.matDCov(3,4) = mulPP*(rep.lY.x()*rep.lY.y() + rep.lZ.x()*rep.lZ.y())
	+ losPP*rep.lX.x()*rep.lX.y();
      svCurrent.matDCov(3,5) = mulPP*(rep.lY.x()*rep.lY.z() + rep.lZ.x()*rep.lZ.z())
	+ losPP*rep.lX.x()*rep.lX.z();
      svCurrent.matDCov(4,4) = mulPP*(rep.lY.y()*rep.lY.y() + rep.lZ.y()*rep.lZ.y())
	+ losPP*rep.lX.y()*rep.lX.y();
      svCurrent.matDCov(4,5) = mulPP*(rep.lY.y()*rep.lY.z() + rep.lZ.y()*rep.lZ.z())
	+ losPP*rep.lX.y()*rep.lX.z();
      svCurrent.matDCov(5,5) = mulPP*(rep.lY.z()*rep.lY.z() + rep.lZ.z()*rep.lZ.z())
	+ losPP*rep.lX.z()*rep.lX.z();

      //still symmetric but fill 9 elements
      svCurrent.matDCov(0,3) = mulRP*(rep.lY.x()*rep.lY.x() + rep.lZ.x()*rep.lZ.x());
      svCurrent.matDCov(0,4) = mulRP*(rep.lY.x()*rep.lY.y() + rep.lZ.x()*rep.lZ.y());
      svCurrent.matDCov(0,5) = mulRP*(rep.lY.x()*rep.lY.z() + rep.lZ.x()*rep.lZ.z());
      svCurrent.matDCov(1,3) = mulRP*(rep.lY.x()*rep.lY.y() + rep.lZ.x()*rep.lZ.y());
      svCurrent.matDCov(1,4) = mulRP*(rep.lY.y()*rep.lY.y() + rep.lZ.y()*rep.lZ.y());
      svCurrent.matDCov(1,5) = mulRP*(rep.lY.y()*rep.lY.z() + rep.lZ.y()*rep.lZ.z());
      svCurrent.matDCov(2,3) = mulRP*(rep.lY.x()*rep.lY.z() + rep.lZ.x()*rep.lZ.z());
      svCurrent.matDCov(2,4) = mulRP*(rep.lY.y()*rep.lY.z() + rep.lZ.y()*rep.lZ.z());
      svCurrent.matDCov(2,5) = mulRP*(rep.lY.z()*rep.lY.z() + rep.lZ.z()*rep.lZ.z());
      
    }
    break;
  }
    //   case POL_1_F:
    //   case POL_2_F:
    //   case POL_M_F:
    //     break;
  default:
    break;
  }

  double pMag = svCurrent.p3.mag();

  if (dir == oppositeToMomentum) dP = -fabs(dP);
  dP = dP > pMag ? pMag-1e-5 : dP;
  getNextState(svCurrent, svNext, dP, tauNext, drVec, dS, dS/radX0,
	       dCTransform_);
  return true;
}

double SteppingHelixPropagator::getDeDx(const SteppingHelixPropagator::StateInfo& sv, 
					double& dEdXPrime, double& radX0,
					MatBounds& rzLims) const{
  radX0 = 1.e24;
  dEdXPrime = 0.;
  rzLims = MatBounds();
  if (noMaterialMode_) return 0;

  double dEdx = 0.;

  double lR = sv.r3.perp();
  double lZ = fabs(sv.r3.z());

  //assume "Iron" .. seems to be quite the same for brass/iron/PbW04
  //good for Fe within 3% for 0.2 GeV to 10PeV
  double p0 = sv.p3.mag();

  //0.065 (PDG) --> 0.044 to better match with MPV
  double dEdX_mat = -(11.4 + 0.96*fabs(log(p0*2.8)) + 0.033*p0*(1.0 - pow(p0, -0.33)) )*1e-3; 
  //in GeV/cm .. 0.8 to get closer to the median or MPV
  double dEdX_HCal = 0.95*dEdX_mat; //extracted from sim
  double dEdX_ECal = 0.45*dEdX_mat;
  double dEdX_coil = 0.35*dEdX_mat; //extracted from sim .. closer to 40% in fact
  double dEdX_Fe =   dEdX_mat;
  double dEdX_MCh =  0.053*dEdX_mat; //chambers on average
  double dEdX_Trk =  0.0114*dEdX_mat;
  double dEdX_Vac =  0.0;

  double radX0_HCal = 1.44/0.8; //guessing
  double radX0_ECal = 0.89/0.7;
  double radX0_coil = 4.; //
  double radX0_Fe =   1.76;
  double radX0_MCh =  1e3; //
  double radX0_Trk =  320.;
  double radX0_Air =  3.e4;
  double radX0_Vac =  3.e9; //"big" number for vacuum


  //not all the boundaries are set below: this will be a default
  if (! (lR < 380 && lZ < 785)){
    if (lZ > 785 ) rzLims = MatBounds(0, 1e4, 785, 1e4);
    if (lZ < 785 ) rzLims = MatBounds(380, 1e4, 0, 785);
  }

  //this should roughly figure out where things are 
  //(numbers taken from Fig1.1.2 TDR and from geom xmls)
  if (lR < 2.9){ //inside beampipe
    dEdx = dEdX_Vac; radX0 = radX0_Vac;
    rzLims = MatBounds(0, 2.9, 0, 785);
  }
  else if (lR < 129){
    if (lZ < 294){ 
      rzLims = MatBounds(2.9, 129, 0, 294);
      dEdx = dEdx = dEdX_Trk; radX0 = radX0_Trk; 
      //somewhat empirical formula that ~ matches the average if going from 0,0,0
      //assuming "uniform" tracker material
      //doesn't really track material layer to layer
      double lEtaDet = fabs(sv.r3.eta());
      double scaleRadX = lEtaDet > 1.5 ? 0.7724 : sin(2.*atan(exp(-0.5*lEtaDet)));
      scaleRadX *= scaleRadX;
      if (lEtaDet > 2 && lZ > 20) scaleRadX *= (lEtaDet-1.);
      if (lEtaDet > 2.5 && lZ > 20) scaleRadX *= (lEtaDet-1.);
      radX0 *= scaleRadX;
    }
    else if (lZ < 372){ 
      rzLims = MatBounds(2.9, 129, 294, 372);
      dEdx = dEdX_ECal; radX0 = radX0_ECal; 
    }//EE averaged out over a larger space
    else if (lZ < 398){
      rzLims = MatBounds(2.9, 129, 372, 398);
      dEdx = dEdX_HCal*0.05; radX0 = radX0_Air; 
    }//betw EE and HE
    else if (lZ < 555){ 
      rzLims = MatBounds(2.9, 129, 398, 555);
      dEdx = dEdX_HCal*0.96; radX0 = radX0_HCal/0.96; 
    } //HE calor abit less dense
    else {
      rzLims = MatBounds(2.9, 129, 555, 785);
      //iron .. don't care about no material in front of HF (too forward)
      if (! (lZ > 568 && lZ < 625 && lR > 85 ) // HE support 
	  && ! (lZ > 785 && lZ < 850 && lR > 118)) {dEdx = dEdX_Fe; radX0 = radX0_Fe; }
      else  { dEdx = dEdX_MCh; radX0 = radX0_MCh; } //ME at eta > 2.2
    }
  }
  else if (lR < 287){
    if (lZ < 372 && lR < 177){ 
      if (lZ < 304) rzLims = MatBounds(129, 177, 0, 304);
      else if (lZ < 343){
	if (lR < 135 ) rzLims = MatBounds(129, 135, 304, 343);
	else if (lR < 172 ){
	  if (lZ < 311 ) rzLims = MatBounds(135, 172, 304, 311);
	  else rzLims = MatBounds(135, 172, 311, 343);
	} else {
	  if (lZ < 328) rzLims = MatBounds(172, 177, 304, 328);
	  else rzLims = MatBounds(172, 177, 328, 343);
	}
      } else {
	if (lR < 156 ) rzLims = MatBounds(129, 156, 343, 372);
	else if ( (lZ - 343) > (lR - 156)*1.38 ) 
	  rzLims = MatBounds(156, 177, 127.73, 372, 0.943726, 0);
	else rzLims = MatBounds(156, 177, 343, 127.73, 0, 0.943726);
      }

      if (!(lR > 135 && lZ <343 && lZ > 304 )
	  && ! (lR > 156 && lZ < 372 && lZ > 343 && ((lZ-343.)< (lR-156.)*1.38)))
	{
	  //the crystals are the same length, but they are not 100% of material
	  double cosThetaEquiv = 0.8/sqrt(1.+lZ*lZ/lR/lR) + 0.2;
	  if (lZ > 343) cosThetaEquiv = 1.;
	  dEdx = dEdX_ECal*cosThetaEquiv; radX0 = radX0_ECal/cosThetaEquiv; 
	} //EB
      else { 
	if ( (lZ > 304 && lZ < 328 && lR < 177 && lR > 135) 
	     && ! (lZ > 311 && lR < 172) ) {dEdx = dEdX_Fe; radX0 = radX0_Fe; } //Tk_Support
	else {dEdx = dEdX_ECal*0.2; radX0 = radX0_Air;} //cables go here
      }
    }
    else if (lZ < 554){ 
      if (lR < 177){
	if ( lZ > 372 && lZ < 398 )rzLims = MatBounds(129, 177, 372, 398);
	else if (lZ < 548) rzLims = MatBounds(129, 177, 398, 548);
	else rzLims = MatBounds(129, 177, 548, 554);
      }
      else if (lR < 193){
	if ((lZ - 307) < (lR - 177.)*1.739) rzLims = MatBounds(177, 193, 0, -0.803, 0, 1.04893);
	else if ( lZ < 389)  rzLims = MatBounds(177, 193, -0.803, 389, 1.04893, 0.);
	else if ( lZ < 548 ) rzLims = MatBounds(177, 193, 389, 548);
	else rzLims = MatBounds(177, 193, 548, 554);
      }
      else if (lR < 264){
	if ( (lZ - 375.7278) < (lR - 193.)/1.327) rzLims = MatBounds(193, 264, 0, 230.287, 0, 0.645788);
	else if ( (lZ - 392.7278) < (lR - 193.)/1.327) 
	  rzLims = MatBounds(193, 264, 230.287, 247.287, 0.645788, 0.645788);
	else if ( lZ < 517 ) rzLims = MatBounds(193, 264, 247.287, 517, 0.645788, 0);
	else if (lZ < 548){
	  if (lR < 246 ) rzLims = MatBounds(193, 246, 517, 548);
	  else rzLims = MatBounds(246, 264, 517, 548);
	}
	else rzLims = MatBounds(193, 264, 548, 554);
      }
      else if ( lR < 275){
	if (lZ < 433) rzLims = MatBounds(264, 275, 0, 433);
	else rzLims = MatBounds(264, 275, 433, 554);
      }
      else {
	if (lZ < 402) rzLims = MatBounds(275, 287, 0, 402);
	else rzLims = MatBounds(275, 287, 402, 554);
      }

      if ((lZ < 433 || lR < 264) && (lZ < 402 || lR < 275) && (lZ < 517 || lR < 246) //notches
	  //I should've made HE and HF different .. now need to shorten HE to match
	  && lZ < 548
	  && ! (lZ < 389 && lZ > 335 && lR < 193 ) //not a gap
	  && ! (lZ > 307 && lZ < 335 && lR < 193 && ((lZ - 307) > (lR - 177.)*1.739)) //not a gap
	  && ! (lR < 177 && lZ < 398) //under the HE nose
	  && ! (lR < 264 && lR > 175 && fabs(441.5 - lZ + (lR - 269.)/1.327) < 8.5) ) //not a gap
	{ dEdx = dEdX_HCal; radX0 = radX0_HCal; }//hcal
      else {dEdx = dEdX_HCal*0.05; radX0 = radX0_Air; }//endcap gap
    }
    else if (lZ < 564){
      if (lR < 251) {
	rzLims = MatBounds(129, 251, 554, 564);      
	dEdx = dEdX_Fe; radX0 = radX0_Fe; 
      }//HE support
      else { 
	rzLims = MatBounds(251, 287, 554, 564);      
	dEdx = dEdX_MCh; radX0 = radX0_MCh; 
      }
    }
    else if (lZ < 625){ 
      rzLims = MatBounds(129, 287, 564, 625);      
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 785){
      if (lR < 275) rzLims = MatBounds(129, 275, 625, 785);
      else {
	if (lZ < 720) rzLims = MatBounds(275, 287, 625, 720);
	else rzLims = MatBounds(275, 287, 720, 785);
      }
      if (! (lR > 275 && lZ < 720)) { dEdx = dEdX_Fe; radX0 = radX0_Fe; }//iron
      else { dEdx = dEdX_MCh; radX0 = radX0_MCh; }
    }
    else if (lZ < 850){
      rzLims = MatBounds(129, 287, 785, 850);
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 910){
      rzLims = MatBounds(129, 287, 850, 910);
      dEdx = dEdX_Fe; radX0 = radX0_Fe; 
    }//iron
    else if (lZ < 975){
      rzLims = MatBounds(129, 287, 910, 975);
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 1000){
      rzLims = MatBounds(129, 287, 975, 1000);
      dEdx = dEdX_Fe; radX0 = radX0_Fe; 
    }//iron
    else { dEdx = 0; radX0 = radX0_Air;}
  }
  else if (lR <380 && lZ < 667){
    if (lZ < 630) rzLims = MatBounds(287, 380, 0, 630);      
    else rzLims = MatBounds(287, 380, 630, 667);  

    if (lZ < 630) { dEdx = dEdX_coil; radX0 = radX0_coil; }//a guess for the solenoid average
    else {dEdx = 0; radX0 = radX0_Air; }//endcap gap
  }
  else {
    if (lZ < 667) {
      double bMag = sv.bf.mag();
      if (lR < 850){
	if (bMag > 0.75 && ! (lZ > 500 && lR <500 && bMag < 1.15)
	    && ! (lZ < 450 && lR > 420 && bMag < 1.15 ) )
	  { dEdx = dEdX_Fe; radX0 = radX0_Fe; }//iron
	else { dEdx = dEdX_MCh; radX0 = radX0_MCh; }
      } else {dEdx = 0; radX0 = radX0_Air; }
    } 
    else if (lR > 750 ){dEdx = 0; radX0 = radX0_Air; }
    else if (lZ < 724){
      if (lR < 380 ) rzLims = MatBounds(287, 380, 667, 724); 
      else rzLims = MatBounds(380, 750, 667, 724); 
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 785){
      if (lR < 380 ) rzLims = MatBounds(287, 380, 724, 785); 
      else rzLims = MatBounds(380, 750, 724, 785); 
      dEdx = dEdX_Fe; radX0 = radX0_Fe; 
    }//iron
    else if (lZ < 850){
      rzLims = MatBounds(287, 750, 785, 850); 
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 910){
      rzLims = MatBounds(287, 750, 850, 910); 
      dEdx = dEdX_Fe; radX0 = radX0_Fe; 
    }//iron
    else if (lZ < 975){
      rzLims = MatBounds(287, 750, 910, 975); 
      dEdx = dEdX_MCh; radX0 = radX0_MCh; 
    }
    else if (lZ < 1000){
      rzLims = MatBounds(287, 750, 975, 1000); 
      dEdx = dEdX_Fe; radX0 = radX0_Fe; 
    }//iron
    else {dEdx = 0; radX0 = radX0_Air; }//air
  }
  
  dEdXPrime = dEdx == 0 ? 0 : -dEdx/dEdX_mat*(0.96/p0 + 0.033 - 0.022*pow(p0, -0.33))*1e-3; //== d(dEdX)/dp

  return dEdx;
}


int SteppingHelixPropagator::cIndex_(int ind) const{
  int result = ind%MAX_POINTS;  
  if (ind != 0 && result == 0){
    result = MAX_POINTS;
  }
  return result;
}

SteppingHelixPropagator::Result
SteppingHelixPropagator::refToDest(SteppingHelixPropagator::DestType dest, 
				   const SteppingHelixPropagator::StateInfo& sv,
				   const double pars[6], 
				   double& dist, double& tanDist, 
				   PropagationDirection& refDirection,
				   double fastSkipDist) const{
  static const std::string metname = "SteppingHelixPropagator";
  Result result = SteppingHelixStateInfo::NOT_IMPLEMENTED;
  double curZ = sv.r3.z();
  double curR = sv.r3.perp();

  switch (dest){
  case RADIUS_DT:
    {
      double cosDPhiPR = cos((sv.r3.deltaPhi(sv.p3)));
      dist = pars[RADIUS_P] - curR;
      if (fabs(dist) > fastSkipDist){
	result = SteppingHelixStateInfo::INACC;
	break;
      }
      tanDist = dist/sv.p3.perp()*sv.p3.mag();
      refDirection = dist*cosDPhiPR > 0 ?
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  case Z_DT:
    {
      dist = pars[Z_P] - curZ;
      if (fabs(dist) > fastSkipDist){
	result = SteppingHelixStateInfo::INACC;
	break;
      }
      tanDist = dist/sv.p3.z()*sv.p3.mag();
      refDirection = sv.p3.z()*dist > 0. ?
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  case PLANE_DT:
    {
      Point rPlane(pars[0], pars[1], pars[2]);
      Vector nPlane(pars[3], pars[4], pars[5]);
      
      double dRDotN = (sv.r3 - rPlane).dot(nPlane);
      
      dist = fabs(dRDotN);
      if (fabs(dist) > fastSkipDist){
	result = SteppingHelixStateInfo::INACC;
	break;
      }
      double p0 = sv.p3.mag();
      double b0 = sv.bf.mag();
      double tN = sv.p3.dot(nPlane)/p0;
      if (fabs(tN)>1e-24) tanDist = -dRDotN/tN;
      if (fabs(tanDist) > 1e4) tanDist = 1e4;
      if (b0>1.5e-6){
	double kVal = 0.0029979*sv.q/p0*b0;
	double aVal = tanDist*kVal;
	Vector lVec = sv.bf.cross(sv.p3); lVec *= 1./b0/p0;
	double bVal = lVec.dot(nPlane)/tN;
	if (fabs(aVal*bVal)< 0.3){
	  double cVal = - sv.bf.cross(lVec).dot(nPlane)/b0/tN; //1- bHat_n*bHat_tau/tau_n;
	  double tanDCorr = bVal/2. + (bVal*bVal/2. + cVal/6)*aVal; 
	  tanDCorr *= aVal;
	  //+ (-bVal/24. + 0.625*bVal*bVal*bVal + 5./12.*bVal*cVal)*aVal*aVal*aVal
	  if (debug_) LogTrace(metname)<<tanDist<<" vs "<<tanDist*(1.+tanDCorr)<<" corr "<<tanDist*tanDCorr<<std::endl;
	  tanDist *= (1.+tanDCorr);
	} else {
	  if (debug_) LogTrace(metname)<<"ABVal "<< fabs(aVal*bVal)
				       <<" = "<<aVal<<" * "<<bVal<<" too large:: will not converge"<<std::endl;
	}
      }
      refDirection = tN*dRDotN < 0. ?
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  case CONE_DT:
    {
      //assumes the cone axis/vertex is along z
      Point cVertex(pars[0], pars[1], pars[2]);
      Vector relV3 = sv.r3 - cVertex;
      double theta(pars[3]);
      if (cVertex.perp() < 1e-5){
	double sinDTheta = sin(theta-relV3.theta());
	double cosDTheta = cos(theta-relV3.theta());
	bool isInside = sin(theta) > sin(relV3.theta()) 
	  && cos(theta)*cos(relV3.theta()) > 0;
	dist = isInside || cosDTheta > 0 ? 
	  relV3.mag()*sinDTheta : relV3.mag();
	if (fabs(dist) > fastSkipDist){
	  result = SteppingHelixStateInfo::INACC;
	  break;
	}
	double normPhi = isInside ? 
	  Geom::pi() - relV3.phi() : relV3.phi();
	double normTheta = theta > Geom::pi()/2. ? 
	  (isInside ? 1.5*Geom::pi() - theta : theta - Geom::pi()/2.) 
	  : (isInside ? Geom::pi()/2 - theta : theta + Geom::pi()/2);
	//this is a normVector from the cone to the point
	Vector norm; norm.setRThetaPhi(fabs(dist), normTheta, normPhi);
	double cosDThetaP = cos(norm.theta() - sv.p3.theta());
	tanDist = dist/fabs(cosDThetaP);
	refDirection = norm.dot(sv.p3) > 0 ?
	  oppositeToMomentum : alongMomentum;
	if (debug_){
	  LogTrace(metname)<<"refToDest:toCone the point is "
			   <<(isInside? "in" : "out")<<"side the cone"
			   <<std::endl;
	}
      }
    }
    break;
    //   case CYLINDER_DT:
    //     break;
  case PATHL_DT:
    {
      double curS = fabs(sv.path_);
      dist = pars[PATHL_P] - curS;
      if (fabs(dist) > fastSkipDist){
	result = SteppingHelixStateInfo::INACC;
	break;
      }
      tanDist = dist;
      refDirection = pars[PATHL_P] > 0 ? 
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  case POINT_PCA_DT:
    {
      Point pDest(pars[0], pars[1], pars[2]);
      dist = (sv.r3 - pDest).mag()+ 1e-24;//add a small number to avoid 1/0
      if (fabs(dist) > fastSkipDist){
	result = SteppingHelixStateInfo::INACC;
	break;
      }
      tanDist = (sv.r3.dot(sv.p3) - pDest.dot(sv.p3))/sv.p3.mag();
      refDirection = tanDist < 0 ?
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  case LINE_PCA_DT:
    {
      Point rLine(pars[0], pars[1], pars[2]);
      Vector dLine(pars[3], pars[4], pars[5]);
      dLine = (dLine - rLine);
      dLine *= 1./dLine.mag();

      Vector dR = sv.r3 - rLine;
      Vector dRPerp = dR - dLine*(dR.dot(dLine));
      dist = dRPerp.mag() + 1e-24;//add a small number to avoid 1/0
      tanDist = dRPerp.dot(sv.p3)/sv.p3.mag();
      //angle wrt line
      double cosAlpha = dLine.dot(sv.p3)/sv.p3.mag();
      tanDist *= fabs(1./sqrt(fabs(1.-cosAlpha*cosAlpha)+1e-96));
      refDirection = tanDist < 0 ?
	alongMomentum : oppositeToMomentum;
      result = SteppingHelixStateInfo::OK;
    }
    break;
  default:
    {
      //some large number
      dist = 1e12;
      tanDist = 1e12;
      refDirection = anyDirection;
      result = SteppingHelixStateInfo::NOT_IMPLEMENTED;
    }
    break;
  }

  if (debug_){
    LogTrace(metname)<<"refToDest input: dest"<<dest<<" pars[]: ";
    for (int i = 0; i < 6; i++){
      LogTrace(metname)<<", "<<i<<" "<<pars[i];
    }
    LogTrace(metname)<<std::endl;
    LogTrace(metname)<<"refToDest output: "
		     <<"\t dist"<< dist
		     <<"\t tanDist"<< tanDist      
		     <<"\t refDirection"<< refDirection
		     <<std::endl;
  }

  return result;
}

SteppingHelixPropagator::Result
SteppingHelixPropagator::refToMagVolume(const SteppingHelixPropagator::StateInfo& sv,
					PropagationDirection dir,
					double& dist, double& tanDist,
					double fastSkipDist) const{

  static const std::string metname = "SteppingHelixPropagator";
  Result result = SteppingHelixStateInfo::NOT_IMPLEMENTED;
  const MagVolume* cVol = sv.magVol;

  if (cVol == 0) return result;
  std::vector<VolumeSide> cVolFaces(cVol->faces());

  double distToFace[6];
  double tanDistToFace[6];
  PropagationDirection refDirectionToFace[6];
  Result resultToFace[6];
  int iFDest = -1;
  
  uint iFDestSorted[6];
  uint nDestSorted =0;

  if (debug_){
    LogTrace(metname)<<"Trying volume "<<DDSolidShapesName::name(cVol->shapeType())
		     <<" with "<<cVolFaces.size()<<" faces"<<std::endl;
  }

  for (uint iFace = 0; iFace < cVolFaces.size(); iFace++){
    if (iFace > 5){
      LogTrace(metname)<<"Too many faces"<<std::endl;
    }
    if (debug_){
      LogTrace(metname)<<"Start with face "<<iFace<<std::endl;
    }
//     const Plane* cPlane = dynamic_cast<const Plane*>(&cVolFaces[iFace].surface());
//     const Cylinder* cCyl = dynamic_cast<const Cylinder*>(&cVolFaces[iFace].surface());
//     const Cone* cCone = dynamic_cast<const Cone*>(&cVolFaces[iFace].surface());
    const Surface* cPlane = 0; //only need to know the loc->glob transform
    const Cylinder* cCyl = 0;
    const Cone* cCone = 0;
    if (typeid(cVolFaces[iFace].surface()) == typeid(const Plane&)){
      cPlane = &cVolFaces[iFace].surface();
    } else if (typeid(cVolFaces[iFace].surface()) == typeid(const Cylinder&)){
      cCyl = reinterpret_cast<const Cylinder*>(&cVolFaces[iFace].surface());
    } else if (typeid(cVolFaces[iFace].surface()) == typeid(const Cone&)){
      cCone = reinterpret_cast<const Cone*>(&cVolFaces[iFace].surface());
    } else {
      edm::LogWarning(metname)<<"Could not cast a volume side surface to a known type"<<std::endl;
    }
    
    if (debug_){
      if (cPlane!=0) LogTrace(metname)<<"The face is a plane at "<<cPlane<<std::endl;
      if (cCyl!=0) LogTrace(metname)<<"The face is a cylinder at "<<cCyl<<std::endl;
    }

    double pars[6];
    DestType dType = UNDEFINED_DT;
    if (cPlane != 0){
      GlobalPoint rPlane = cPlane->position();
      // = cPlane->toGlobal(LocalVector(0,0,1.)); nPlane = nPlane.unit();
      GlobalVector nPlane(cPlane->rotation().zx(), cPlane->rotation().zy(), cPlane->rotation().zz());
      
      if (sv.r3.z() < 0){
	pars[0] = rPlane.x(); pars[1] = rPlane.y(); pars[2] = rPlane.z();
	pars[3] = nPlane.x(); pars[4] = nPlane.y(); pars[5] = nPlane.z();
      } else {
	pars[0] = rPlane.x(); pars[1] = rPlane.y(); pars[2] = -rPlane.z();
	pars[3] = nPlane.x(); pars[4] = nPlane.y(); pars[5] = -nPlane.z();
      }
      dType = PLANE_DT;
    } else if (cCyl != 0){
      if (debug_){
	LogTrace(metname)<<"Cylinder at "<<cCyl->position()
			 <<" rorated by "<<cCyl->rotation()
			 <<std::endl;
      }
      pars[RADIUS_P] = cCyl->radius();
      dType = RADIUS_DT;
    } else if (cCone != 0){
      if (debug_){
	LogTrace(metname)<<"Cone at "<<cCone->position()
			 <<" rorated by "<<cCone->rotation()
			 <<" vertex at "<<cCone->vertex()
			 <<" angle of "<<cCone->openingAngle()
			 <<std::endl;
      }
      if (sv.r3.z() < 0){
	pars[0] = cCone->vertex().x(); pars[1] = cCone->vertex().y(); 
	pars[2] = cCone->vertex().z();
	pars[3] = cCone->openingAngle();
      } else {
	pars[0] = cCone->vertex().x(); pars[1] = cCone->vertex().y(); 
	pars[2] = -cCone->vertex().z();
	pars[3] = Geom::pi() - cCone->openingAngle();
      }
      dType = CONE_DT;
    } else {
      LogTrace(metname)<<"Unknown surface"<<std::endl;
      resultToFace[iFace] = SteppingHelixStateInfo::UNDEFINED;
      continue;
    }
    resultToFace[iFace] = 
      refToDest(dType, sv, pars, 
		distToFace[iFace], tanDistToFace[iFace], refDirectionToFace[iFace], fastSkipDist);    
    
    if (resultToFace[iFace] != SteppingHelixStateInfo::OK){
      if (resultToFace[iFace] == SteppingHelixStateInfo::INACC) result = SteppingHelixStateInfo::INACC;
      continue;
    }
    
    //keep those in right direction for later use
    if (refDirectionToFace[iFace] == dir || fabs(distToFace[iFace]/tanDistToFace[iFace]) < 2e-2){
      iFDestSorted[nDestSorted] = iFace;
      nDestSorted++;
    }
    if (debug_) 
      LogTrace(metname)<<cVol<<" "<<iFace<<" "
		       <<tanDistToFace[iFace]<<" "<<distToFace[iFace]<<" "<<refDirectionToFace[iFace]<<" "<<dir<<std::endl;
  }
  
  for (uint i = 0;i<nDestSorted; ++i){
    uint iMax = nDestSorted-i-1;
    for (uint j=0;j<nDestSorted-i; ++j){
      if (fabs(tanDistToFace[iFDestSorted[j]]) > fabs(tanDistToFace[iFDestSorted[iMax]]) ){
	iMax = j;
      }
    }
    uint iTmp = iFDestSorted[nDestSorted-i-1];
    iFDestSorted[nDestSorted-i-1] = iFDestSorted[iMax];
    iFDestSorted[iMax] = iTmp;
  }
  if (debug_){
    for (uint i=0;i<nDestSorted;++i){
      LogTrace(metname)<<cVol<<" "<<i<<" "<<iFDestSorted[i]<<" "<<tanDistToFace[iFDestSorted[i]]<<std::endl;
    }
  }

  //now go from the shortest to the largest distance hoping to get a point in the volume.
  //other than in case of a near-parallel travel this should stop after the first try
  for (uint i=0; i<nDestSorted;++i){
    iFDest = iFDestSorted[i];

    double sign = dir == alongMomentum ? 1. : -1.;
    Point gPointEst(sv.r3);
    Vector lDelta(sv.p3); lDelta *= 1./sv.p3.mag()*sign*sqrt(fabs(distToFace[iFDest]*tanDistToFace[iFDest]));
    gPointEst += lDelta;
    if (debug_){
      LogTrace(metname)<<"Linear est point "<<gPointEst
		       <<" for iFace "<<iFDest<<std::endl;
    }
    GlobalPoint gPointEstNegZ(gPointEst.x(), gPointEst.y(), -fabs(gPointEst.z()));
    if ( cVol->inside(gPointEstNegZ) ){
      if (debug_){
	LogTrace(metname)<<"The point is inside the volume"<<std::endl;
      }
      
      result = SteppingHelixStateInfo::OK;
      dist = distToFace[iFDest];
      tanDist = tanDistToFace[iFDest];
      if (debug_){
	LogTrace(metname)<<"Got a point near closest boundary -- face "<<iFDest<<std::endl;
      }
      break;
    }
  }

  if (iFDest == -1){
    if (debug_) LogTrace(metname)<<"Failed to find a dest point inside the volume"<<std::endl;
  }

  return result;
}


SteppingHelixPropagator::Result
SteppingHelixPropagator::refToMatVolume(const SteppingHelixPropagator::StateInfo& sv,
					PropagationDirection dir,
					double& dist, double& tanDist,
					double fastSkipDist) const{

  static const std::string metname = "SteppingHelixPropagator";
  Result result = SteppingHelixStateInfo::NOT_IMPLEMENTED;

  double parLim[6] = {sv.rzLims.rMin, sv.rzLims.rMax, 
		      sv.rzLims.zMin, sv.rzLims.zMax, 
		      sv.rzLims.th1, sv.rzLims.th2 };

  double distToFace[4];
  double tanDistToFace[4];
  PropagationDirection refDirectionToFace[4];
  Result resultToFace[4];
  int iFDest = -1;
  
  for (uint iFace = 0; iFace < 4; iFace++){
    if (debug_){
      LogTrace(metname)<<"Start with mat face "<<iFace<<std::endl;
    }

    double pars[6];
    DestType dType = UNDEFINED_DT;
    if (iFace > 1){
      if (fabs(parLim[iFace+2])< 1e-6){//plane
	if (sv.r3.z() < 0){
	  pars[0] = 0; pars[1] = 0; pars[2] = -parLim[iFace];
	  pars[3] = 0; pars[4] = 0; pars[5] = 1;
	} else {
	  pars[0] = 0; pars[1] = 0; pars[2] = parLim[iFace];
	  pars[3] = 0; pars[4] = 0; pars[5] = 1;
	}
	dType = PLANE_DT;
      } else {
	if (sv.r3.z() > 0){
	  pars[0] = 0; pars[1] = 0; 
	  pars[2] = parLim[iFace];
	  pars[3] = Geom::pi()/2. - parLim[iFace+2];
	} else {
	  pars[0] = 0; pars[1] = 0; 
	  pars[2] = - parLim[iFace];
	  pars[3] = Geom::pi()/2. + parLim[iFace+2];
	}
	dType = CONE_DT;	
      }
    } else {
      pars[RADIUS_P] = parLim[iFace];
      dType = RADIUS_DT;
    }

    resultToFace[iFace] = 
      refToDest(dType, sv, pars, 
		distToFace[iFace], tanDistToFace[iFace], refDirectionToFace[iFace], fastSkipDist);
    
    if (resultToFace[iFace] != SteppingHelixStateInfo::OK){
      if (resultToFace[iFace] == SteppingHelixStateInfo::INACC) result = SteppingHelixStateInfo::INACC;
      continue;
    }
    if (refDirectionToFace[iFace] == dir || fabs(distToFace[iFace]/tanDistToFace[iFace]) < 2e-2){
      double sign = dir == alongMomentum ? 1. : -1.;
      Point gPointEst(sv.r3);
      Vector lDelta(sv.p3); lDelta *= sign*sqrt(fabs(distToFace[iFace]*tanDistToFace[iFace]))/sv.p3.mag();
      gPointEst += lDelta;
      if (debug_){
	LogTrace(metname)<<"Linear est point "<<gPointEst
			 <<std::endl;
      }
      double lZ = fabs(gPointEst.z());
      double lR = gPointEst.perp();
      if ( (lZ - parLim[2]) > lR*tan(parLim[4]) 
	   && (lZ - parLim[3]) < lR*tan(parLim[5])  
	   && lR > parLim[0] && lR < parLim[1]){
	if (debug_){
	  LogTrace(metname)<<"The point is inside the volume"<<std::endl;
	}
	//OK, guessed a point still inside the volume
	if (iFDest == -1){
	  iFDest = iFace;
	} else {
	  if (fabs(tanDistToFace[iFDest]) > fabs(tanDistToFace[iFace])){
	    iFDest = iFace;
	  }
	}
      } else {
	if (debug_){
	  LogTrace(metname)<<"The point is NOT inside the volume"<<std::endl;
	}
      }
    }

  }
  if (iFDest != -1){
    result = SteppingHelixStateInfo::OK;
    dist = distToFace[iFDest];
    tanDist = tanDistToFace[iFDest];
    if (debug_){
      LogTrace(metname)<<"Got a point near closest boundary -- face "<<iFDest<<std::endl;
    }
  } else {
    if (debug_){
      LogTrace(metname)<<"Failed to find a dest point inside the volume"<<std::endl;
    }
  }

  return result;
}


