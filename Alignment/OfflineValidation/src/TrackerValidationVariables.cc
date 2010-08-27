
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "CondFormats/Alignment/interface/Definitions.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"

#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"

#include "TMath.h"

TrackerValidationVariables::TrackerValidationVariables()
{

}

TrackerValidationVariables::TrackerValidationVariables(const edm::EventSetup& es, const edm::ParameterSet& iSetup) 
  : conf_(iSetup)
{
  es.get<TrackerDigiGeometryRecord>().get(tkGeom_);
  es.get<IdealMagneticFieldRecord>().get(magneticField_);
}

TrackerValidationVariables::~TrackerValidationVariables()
{

}

void
TrackerValidationVariables::fillHitQuantities(const Trajectory* trajectory, std::vector<AVHitStruct> & v_avhitout)
{
  TrajectoryStateCombiner tsoscomb;

  const std::vector<TrajectoryMeasurement> &tmColl = trajectory->measurements();
  for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin();
      itTraj != tmColl.end();
      ++itTraj) {
    
    if (!itTraj->updatedState().isValid()) continue;
    
    TrajectoryStateOnSurface tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState() );
    TransientTrackingRecHit::ConstRecHitPointer hit = itTraj->recHit();
    
    if(!hit->isValid() || hit->geographicalId().det() != DetId::Tracker) continue;
    
    AVHitStruct hitStruct;
    const DetId& hit_detId = hit->geographicalId();
    uint IntRawDetID = (hit_detId.rawId());	
    uint IntSubDetID = (hit_detId.subdetId());
    
    if(IntSubDetID == 0) continue;
    
    //first calculate residuals in cartesian coordinates in the local module coordinate system
    
    LocalPoint lPHit = hit->localPosition();
    LocalPoint lPTrk = tsos.localPosition();
    LocalVector lVTrk = tsos.localDirection();

    hitStruct.localAlpha = atan2(lVTrk.x(), lVTrk.z()); // wrt. normal tg(alpha)=x/z
    hitStruct.localBeta  = atan2(lVTrk.y(), lVTrk.z()); // wrt. normal tg(beta)= y/z

    //LocalError errHit = hit->localPositionError();
    // adding APE to hitError
    AlgebraicROOTObject<2>::SymMatrix mat = asSMatrix<2>(hit->parametersError());
    LocalError errHit = LocalError( mat(0,0),mat(0,1),mat(1,1) );
    LocalError errTrk = tsos.localError().positionError();
    
    //check for negative error values: track error can have negative value, if matrix inversion fails (very rare case)
    //hit error should always give positive values
    if(errHit.xx()<0. || errHit.yy()<0. || errTrk.xx()<0. || errTrk.yy()<0.){
      edm::LogError("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillHitQuantities"
						  << "One of the squared error methods gives negative result"
						  << "\n\terrHit.xx()\terrHit.yy()\terrTrk.xx()\terrTrk.yy()"
						  << "\n\t" << errHit.xx()
						  << "\t" << errHit.yy()
						  << "\t" << errTrk.xx()
						  << "\t" << errTrk.yy();
      continue;
    }
    
    align::LocalVector res = lPTrk - lPHit;
    
    float resXErr = std::sqrt( errHit.xx() + errTrk.xx() );
    float resYErr = std::sqrt( errHit.yy() + errTrk.yy() );
    
    hitStruct.resX = res.x();
    hitStruct.resY = res.y();
    hitStruct.resErrX = resXErr;
    hitStruct.resErrY = resYErr;

    // hitStruct.localX = lPhit.x();
    // hitStruct.localY = lPhit.y();
    // EM: use predictions for local coordinates
    hitStruct.localX = lPTrk.x();
    hitStruct.localY = lPTrk.y();

    // now calculate residuals taking global orientation of modules and radial topology in TID/TEC into account
    float resXprime(999.F), resYprime(999.F);
    float resXprimeErr(999.F), resYprimeErr(999.F);
    
    if(hit->detUnit()){ // is it a single physical module?
      const GeomDetUnit& detUnit = *(hit->detUnit());
      float uOrientation(-999.F), vOrientation(-999.F);
      float resXTopol(999.F), resYTopol(999.F);
      
      const Surface& surface = hit->detUnit()->surface();
      const BoundPlane& boundplane = hit->detUnit()->surface();
      const Bounds& bound = boundplane.bounds();
      
      float length = 0;
      float width = 0;
      LocalPoint lPModule(0.,0.,0.), lUDirection(1.,0.,0.), lVDirection(0.,1.,0.);
      GlobalPoint gPModule    = surface.toGlobal(lPModule),
	          gUDirection = surface.toGlobal(lUDirection),
	          gVDirection = surface.toGlobal(lVDirection);
      
      if (IntSubDetID == PixelSubdetector::PixelBarrel ||
	  IntSubDetID == StripSubdetector::TIB || 
	  IntSubDetID == StripSubdetector::TOB) {
	
	uOrientation = deltaPhi(gUDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	vOrientation = gVDirection.z() - gPModule.z() >= 0 ? +1.F : -1.F;
	resXTopol = res.x();
	resYTopol = res.y();
	resXprimeErr = resXErr;
	resYprimeErr = resYErr;

	const RectangularPlaneBounds *rectangularBound = dynamic_cast<const RectangularPlaneBounds*>(&bound);
	length = rectangularBound->length();
	width = rectangularBound->width();
	hitStruct.localXnorm = 2*hitStruct.localX/width;
	hitStruct.localYnorm = 2*hitStruct.localY/length;
	
      } else if (IntSubDetID == PixelSubdetector::PixelEndcap) {
	
	uOrientation = gUDirection.perp() - gPModule.perp() >= 0 ? +1.F : -1.F;
	vOrientation = deltaPhi(gVDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	resXTopol = res.x();
	resYTopol = res.y();
	resXprimeErr = resXErr;
	resYprimeErr = resYErr;
	
	const RectangularPlaneBounds *rectangularBound = dynamic_cast<const RectangularPlaneBounds*>(&bound);
	length = rectangularBound->length();
	width = rectangularBound->width();
	hitStruct.localXnorm = 2*hitStruct.localX/width;
	hitStruct.localYnorm = 2*hitStruct.localY/length;
	
      } else if (IntSubDetID == StripSubdetector::TID ||
		 IntSubDetID == StripSubdetector::TEC) {
	
	uOrientation = deltaPhi(gUDirection.phi(),gPModule.phi()) >= 0. ? +1.F : -1.F;
	vOrientation = gVDirection.perp() - gPModule.perp() >= 0. ? +1.F : -1.F;
	
	if (!dynamic_cast<const RadialStripTopology*>(&detUnit.topology()))continue;
	const RadialStripTopology& topol = dynamic_cast<const RadialStripTopology&>(detUnit.topology());
	
	MeasurementPoint measHitPos = topol.measurementPosition(lPHit);
	MeasurementPoint measTrkPos = topol.measurementPosition(lPTrk);
	
	MeasurementError measHitErr = topol.measurementError(lPHit,errHit);
	MeasurementError measTrkErr = topol.measurementError(lPTrk,errTrk);
	
	if (measHitErr.uu()<0. ||
	    measHitErr.vv()<0. ||
	    measTrkErr.uu()<0. ||
	    measTrkErr.vv()<0.){
	  edm::LogError("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillHitQuantities"
						      << "One of the squared error methods gives negative result"
						      << "\n\tmeasHitErr.uu()\tmeasHitErr.vv()\tmeasTrkErr.uu()\tmeasTrkErr.vv()"
						      << "\n\t" << measHitErr.uu()
						      << "\t" << measHitErr.vv()
						      << "\t" << measTrkErr.uu()
						      << "\t" << measTrkErr.vv();
	  continue;
	}
	  
	float localStripLengthHit = topol.localStripLength(lPHit);
	float localStripLengthTrk = topol.localStripLength(lPTrk);
	float phiHit = topol.stripAngle(measHitPos.x());
	float phiTrk = topol.stripAngle(measTrkPos.x());
	float r_0 = topol.originToIntersection();
	
	resXTopol = (phiTrk-phiHit)*r_0;
	//resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit;
	float cosPhiHit(cos(phiHit)), cosPhiTrk(cos(phiTrk)),
	      sinPhiHit(sin(phiHit)), sinPhiTrk(sin(phiTrk));
	float l_0 = r_0 - topol.detHeight()/2;
	resYTopol = measTrkPos.y()*localStripLengthTrk - measHitPos.y()*localStripLengthHit + l_0*(1/cosPhiTrk - 1/cosPhiHit);
	
	resXprimeErr = std::sqrt(measHitErr.uu()+measTrkErr.uu())*topol.angularWidth()*r_0;
	//resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk);
	float helpSummand = l_0*l_0*topol.angularWidth()*topol.angularWidth()*(sinPhiHit*sinPhiHit/pow(cosPhiHit,4)*measHitErr.uu()
									       + sinPhiTrk*sinPhiTrk/pow(cosPhiTrk,4)*measTrkErr.uu() );
	resYprimeErr = std::sqrt(measHitErr.vv()*localStripLengthHit*localStripLengthHit
				 + measTrkErr.vv()*localStripLengthTrk*localStripLengthTrk + helpSummand );  

	const TrapezoidalPlaneBounds *trapezoidalBound = dynamic_cast<const TrapezoidalPlaneBounds*>(&bound);
	length = trapezoidalBound->length();
	width = trapezoidalBound->widthAtHalfLength();
	hitStruct.localXnorm = 2*hitStruct.localX/width;
	hitStruct.localYnorm = 2*hitStruct.localY/length;

      } else {
	edm::LogWarning("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillHitQuantities" 
						      << "No valid tracker subdetector " << IntSubDetID;
	continue;
      }  
      
      resXprime = resXTopol*uOrientation;
      resYprime = resYTopol*vOrientation;
      
    } else { // not a detUnit, so must be a virtual 2D-Module
      // FIXME: at present only for det units residuals are calculated and filled in the hitStruct
      // But in principle this method should also be useable for the gluedDets (2D modules in TIB, TID, TOB, TEC)
      // In this case, only orientation should be taken into account for primeResiduals, but not the radial topology
      // At present, default values (999.F) are given out
    }
    
    hitStruct.resXprime = resXprime;
    hitStruct.resYprime = resYprime;
    hitStruct.resXprimeErr = resXprimeErr;
    hitStruct.resYprimeErr = resYprimeErr;
    
    hitStruct.rawDetId = IntRawDetID;
    hitStruct.phi = tsos.globalDirection().phi();
    hitStruct.eta = tsos.globalDirection().eta();
    
    // first try for overlap residuals
    // based on Code from Keith and Wolfgang
    if (itTraj+1 != tmColl.end()) {
      TransientTrackingRecHit::ConstRecHitPointer hit2 = (itTraj+1)->recHit();
      TrackerAlignableId ali1, ali2;
      if (hit2->isValid() && 
	  ali1.typeAndLayerFromDetId(hit->geographicalId()) == ali2.typeAndLayerFromDetId(hit2->geographicalId())  &&
	  hit2->geographicalId().rawId() !=  SiStripDetId(IntRawDetID).partnerDetId()) {	    
	  
	float overlapPath_;
	TrajectoryStateCombiner combiner_;
	AnalyticalPropagator propagator(&(*magneticField_));
	// forward and backward predicted states at module 1
	TrajectoryStateOnSurface fwdPred1 = (itTraj)->forwardPredictedState();
	TrajectoryStateOnSurface bwdPred1 = (itTraj)->backwardPredictedState();
	if ( !fwdPred1.isValid() || !bwdPred1.isValid() ) continue;
	// backward predicted state at module 2
	TrajectoryStateOnSurface bwdPred2 = (itTraj+1)->backwardPredictedState();
	TrajectoryStateOnSurface fwdPred2 = (itTraj+1)->forwardPredictedState();
	if ( !bwdPred2.isValid() ) continue;
	// extrapolation bwdPred2 to module 1
	TrajectoryStateOnSurface bwdPred2At1 = propagator.propagate(bwdPred2,fwdPred1.surface());
	if ( !bwdPred2At1.isValid() ) continue;
	// combination with fwdPred1 (ref. state, best estimate without hits 1 and 2)
	TrajectoryStateOnSurface comb1 = combiner_.combine(fwdPred1,bwdPred2At1);
	if ( !comb1.isValid() ) continue;
	  
	//
	// propagation of reference parameters to module 2
	//
	std::pair<TrajectoryStateOnSurface,double> tsosWithS =
	  propagator.propagateWithPath(comb1,bwdPred2.surface());
	TrajectoryStateOnSurface comb1At2 = tsosWithS.first;
	
	// Alternative possibility, not used at present
	// TrajectoryStateOnSurface comb1At2 = propagator.propagate(comb1,bwdPred2.surface());
	  
	if ( !comb1At2.isValid() ) continue;
	overlapPath_ = tsosWithS.second;
	  
	std::vector<GlobalPoint> predictedPositions;
	predictedPositions.push_back(comb1.globalPosition());
	predictedPositions.push_back(comb1At2.globalPosition());
	
	GlobalVector diff_pred = predictedPositions[0] - predictedPositions[1];
	
	TrajectoryStateOnSurface tsos2 = tsoscomb( (itTraj+1)->forwardPredictedState(), (itTraj+1)->backwardPredictedState() );
	align::LocalVector res2 = tsos2.localPosition() - hit2->localPosition();
	//float overlapresidual = res2.x() - res.x();
	float overlapresidual = diff_pred.x();
	
	hitStruct.overlapres = std::make_pair(hit2->geographicalId().rawId(),overlapresidual);
      }
    }
   
    v_avhitout.push_back(hitStruct);
  } 
}

void
TrackerValidationVariables::fillTrackQuantities(const edm::Event& event, std::vector<AVTrackStruct> & v_avtrackout)
{
  edm::Handle<TrajTrackAssociationCollection> TrajTracksMap;
  if (event.getByLabel(conf_.getParameter<edm::InputTag>("Tracks"), TrajTracksMap)) {
 
    const Trajectory* trajectory;
    const reco::Track* track;

    for ( TrajTrackAssociationCollection::const_iterator iPair = TrajTracksMap->begin();
          iPair != TrajTracksMap->end();
	  ++iPair) {

      trajectory = &(*(*iPair).key);
      track = &(*(*iPair).val);

      AVTrackStruct trackStruct;
      LogDebug("TrackerValidationVariables") << "track pt " << track->p();
     
      trackStruct.p = track->p();
      trackStruct.pt = track->pt();
      trackStruct.ptError = track->ptError();
      trackStruct.px = track->px();
      trackStruct.py = track->py();
      trackStruct.pz = track->pz();
      trackStruct.eta = track->eta();
      trackStruct.phi = track->phi();
      trackStruct.chi2 = track->chi2();
      trackStruct.chi2Prob= TMath::Prob(track->chi2(),track->ndof());
      trackStruct.normchi2 = track->normalizedChi2();
      GlobalPoint gPoint(track->vx(), track->vy(), track->vz());
      double theLocalMagFieldInInverseGeV = magneticField_->inInverseGeV(gPoint).z();
      trackStruct.kappa = -track->charge()*theLocalMagFieldInInverseGeV/track->pt();
      trackStruct.charge = track->charge();
      trackStruct.d0 = track->d0();
      trackStruct.dz = track->dz();
      trackStruct.numberOfValidHits = track->numberOfValidHits();
      trackStruct.numberOfLostHits = track->numberOfLostHits();

      fillHitQuantities(trajectory, trackStruct.hits);
      
      v_avtrackout.push_back(trackStruct);
    }
    
  } else {

    edm::LogError("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillQuantities" 
						<< "No track collection found: skipping event";
  }
}

void
TrackerValidationVariables::fillHitQuantities(const edm::Event& event, std::vector<AVHitStruct> & v_avhitout)
{
  edm::Handle<TrajTrackAssociationCollection> TrajTracksMap;
  if (event.getByLabel(conf_.getParameter<edm::InputTag>("Tracks"), TrajTracksMap)) {
 
    const Trajectory* trajectory;

    for ( TrajTrackAssociationCollection::const_iterator iPair = TrajTracksMap->begin();
          iPair != TrajTracksMap->end();
	  ++iPair) {

      trajectory = &(*(*iPair).key);
      
      fillHitQuantities(trajectory, v_avhitout);
    }
  }
}
