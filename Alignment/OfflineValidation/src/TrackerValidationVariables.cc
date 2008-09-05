
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "CLHEP/Vector/RotationInterfaces.h" 

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"
#include "RecoTracker/TrackProducer/interface/TrackingRecHitLessFromGlobalPosition.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementVector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/RadialStripTopology.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"  
#include "Alignment/TrackerAlignment/interface/TrackerAlignableId.h"
#include "Alignment/OfflineValidation/interface/TrackerValidationVariables.h"


TrackerValidationVariables::TrackerValidationVariables(){}


TrackerValidationVariables::TrackerValidationVariables(const edm::EventSetup& es, const edm::ParameterSet& iSetup) 
  : conf_(iSetup)
{
  es.get<TrackerDigiGeometryRecord>().get( tkGeom_ );
  es.get<IdealMagneticFieldRecord>().get(magneticField_);

}

TrackerValidationVariables::~TrackerValidationVariables() {}

void 
TrackerValidationVariables::fillHitQuantities(const edm::Event& iEvent, 
				      std::vector<AVHitStruct>& v_avhitout )
{
  edm::Handle<std::vector<Trajectory> > trajCollectionHandle;
  iEvent.getByLabel(conf_.getParameter<std::string>("trajectoryInput"),trajCollectionHandle);
  
  TrajectoryStateCombiner tsoscomb;
  LogDebug("TrackerValidationVariables") << "trajColl->size(): " << trajCollectionHandle->size() ;
  for(std::vector<Trajectory>::const_iterator it = trajCollectionHandle->begin(), itEnd = trajCollectionHandle->end(); 
      it!=itEnd;++it){
    std::vector<TrajectoryMeasurement> tmColl = it->measurements();
    for(std::vector<TrajectoryMeasurement>::const_iterator itTraj = tmColl.begin(), itTrajEnd = tmColl.end(); 
	itTraj != itTrajEnd; ++itTraj) {


      if(! itTraj->updatedState().isValid()) continue;
      
      
      TrajectoryStateOnSurface tsos = tsoscomb( itTraj->forwardPredictedState(), itTraj->backwardPredictedState() );
      TransientTrackingRecHit::ConstRecHitPointer hit = itTraj->recHit();
      if(! hit->isValid() || hit->geographicalId().det() != DetId::Tracker ) {
	continue; 
      } else {
	AVHitStruct hitStruct;
	const DetId & hit_detId = hit->geographicalId();
	uint IntRawDetID = (hit_detId.rawId());	
	uint IntSubDetID = (hit_detId.subdetId());
	
	if(IntSubDetID == 0 )
	  continue;
	
	align::LocalVector res = tsos.localPosition() - hit->localPosition();

	LocalError err1 = tsos.localError().positionError();
	LocalError err2 = hit->localPositionError();
	
	float errX = std::sqrt( err1.xx() + err2.xx() );
	float errY = std::sqrt( err1.yy() + err2.yy() );
	
	LogDebug("TrackerValidationVariables") << "Residual x/y " << res.x() << '/' << res.y() 
					       << ", Error x/y " << errX << '/' << errY;		

	// begin partly copied from Tifanalyser 

	const GeomDetUnit* detUnit = hit->detUnit();
	double dPhi = -999, dR = -999, dZ = -999, phiorientation = -999, zorientation = -999;
	double R = 0.;
	double origintointersect = 0.;	

	hitStruct.resX = res.x();
	hitStruct.resErrX = errX;
	hitStruct.resErrY = errY;

	if(detUnit) {
	  const Surface& surface = hit->detUnit()->surface();
	  LocalPoint lPModule(0.,0.,0.), lPhiDirection(1.,0.,0.), lROrZDirection(0.,1.,0.);
	  GlobalPoint gPModule       = surface.toGlobal(lPModule),
	    gPhiDirection  = surface.toGlobal(lPhiDirection),
	    gROrZDirection = surface.toGlobal(lROrZDirection);
	  phiorientation = deltaPhi(gPhiDirection.phi(),gPModule.phi()) >= 0 ? +1. : -1.;
	  zorientation = gROrZDirection.z() - gPModule.z() >= 0 ? +1. : -1.;
	  dPhi = tsos.globalPosition().phi() - hit->globalPosition().phi();
	  
	  if(IntSubDetID == PixelSubdetector::PixelBarrel || IntSubDetID == PixelSubdetector::PixelEndcap || 
	     IntSubDetID == StripSubdetector::TIB || 
	     IntSubDetID == StripSubdetector::TOB) {
	    hitStruct.resXprime = (res.x())*(phiorientation );
	    hitStruct.resXprimeErr = errX;
	    hitStruct.resYprime = (res.y())*(zorientation );
	    hitStruct.resYprimeErr = errY;
	    dZ = gROrZDirection.z() - gPModule.z();
	  } else if (IntSubDetID == StripSubdetector::TID || IntSubDetID == StripSubdetector::TEC) {
	    const RadialStripTopology* theTopol = dynamic_cast<const RadialStripTopology*>(&(detUnit->topology()));
	    origintointersect =  static_cast<float>(theTopol->originToIntersection());
	    

	    MeasurementPoint theMeasHitPos = theTopol->measurementPosition(hit->localPosition());
	    MeasurementPoint theMeasStatePos = theTopol->measurementPosition(tsos.localPosition());
	    Measurement2DVector residual =  theMeasStatePos - theMeasHitPos;
	    
	    MeasurementError theMeasHitErr = theTopol->measurementError(hit->localPosition(),err2);
	    MeasurementError theMeasStateErr = theTopol->measurementError(tsos.localPosition(),err1);

	    double localPitch = theTopol->localPitch(hit->localPosition());
	    double localStripLength =  theTopol->localStripLength(hit->localPosition());
	    float xPrime = residual.x()*localPitch ;
	    float measErr = std::sqrt( theMeasHitErr.uu()+theMeasStateErr.uu()) *localPitch;
	    float measErrY = std::sqrt( theMeasHitErr.vv()+theMeasStateErr.vv()) *localStripLength;

	    R = origintointersect;
	    dR = theTopol->yDistanceToIntersection( tsos.localPosition().y()) - 
	      theTopol->yDistanceToIntersection( hit->localPosition().y());
	    
	    hitStruct.resXprime = xPrime;
	    hitStruct.resXprimeErr = measErr;
	    hitStruct.resYprime = residual.y()*localStripLength;
	    hitStruct.resYprimeErr =  measErrY * localStripLength;
	    
	  } else {
	    edm::LogWarning("TrackerValidationVariables") << "@SUB=TrackerValidationVariables::fillHitQuantities" 
							  << "No valid tracker subdetector " << IntSubDetID;
	    hitStruct.resXprime = -999;
	  }	 
	  
	}
	
	
	if(dR != -999) hitStruct.resY = dR;
	else if(dZ != -999) hitStruct.resY = res.y() * (dZ >=0.? +1 : -1) ;
	else hitStruct.resY = res.y();
	
	hitStruct.rawDetId = IntRawDetID;
	hitStruct.phi = tsos.globalDirection().phi();

	// first try for overlapp residuals
	// based on Code from Keith and Wolfgang
	if(itTraj+1 != itTrajEnd) {
	  TransientTrackingRecHit::ConstRecHitPointer hit2 = (itTraj+1)->recHit();
	  TrackerAlignableId ali1, ali2;
	  if(hit2->isValid() && 
	     ali1.typeAndLayerFromDetId(hit->geographicalId()) == ali2.typeAndLayerFromDetId(hit2->geographicalId())  &&
	     hit2->geographicalId().rawId() !=  SiStripDetId::SiStripDetId(IntRawDetID).partnerDetId()  
	     ) {	    
	    
	    float overlapPath_;
	    TrajectoryStateCombiner combiner_;
	    AnalyticalPropagator propagator(&(*magneticField_));
	    // forward and backward predicted states at module 1
	    TrajectoryStateOnSurface fwdPred1 = (itTraj)->forwardPredictedState();
	    TrajectoryStateOnSurface bwdPred1 = (itTraj)->backwardPredictedState();
	    if ( !fwdPred1.isValid() || !bwdPred1.isValid() )  continue;
	    // backward predicted state at module 2
	    TrajectoryStateOnSurface bwdPred2 = (itTraj+1)->backwardPredictedState();
	    TrajectoryStateOnSurface fwdPred2 = (itTraj+1)->forwardPredictedState();
	    if ( !bwdPred2.isValid() )  continue;
	    // extrapolation bwdPred2 to module 1
	    TrajectoryStateOnSurface bwdPred2At1 = propagator.propagate(bwdPred2,fwdPred1.surface());
	    if ( !bwdPred2At1.isValid() )  continue;
	    // combination with fwdPred1 (ref. state, best estimate without hits 1 and 2)
	    TrajectoryStateOnSurface comb1 = combiner_.combine(fwdPred1,bwdPred2At1);
	    if ( !comb1.isValid() )  continue;
	    
	    //
	    // propagation of reference parameters to module 2
	    //
	    std::pair<TrajectoryStateOnSurface,double> tsosWithS =
	      propagator.propagateWithPath(comb1,bwdPred2.surface());
	    TrajectoryStateOnSurface comb1At2 = tsosWithS.first;
	    
	    // Alternative possibility, not used at present
	    //TrajectoryStateOnSurface comb1At2 = propagator.propagate(comb1,bwdPred2.surface());
	    
	    if ( !comb1At2.isValid() )  continue;
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
  }  

}

void 
TrackerValidationVariables::fillTrackQuantities(const edm::Event& iEvent,
					std::vector<AVTrackStruct>& v_avtrackout)
{
  // get track collection from the event
  edm::InputTag TkTag = conf_.getParameter<edm::InputTag>("Tracks");
  edm::Handle<reco::TrackCollection> RecoTracks;
  iEvent.getByLabel(TkTag,RecoTracks);
  LogDebug("TrackerValidationVariables")<<"track collection size "<< RecoTracks->size();
  
  // Put here all track based quantities such as eta, phi, pt,.... 
  int i=0;
  for( reco::TrackCollection::const_iterator RecoTrack = RecoTracks->begin(), RecoTrackEnd = RecoTracks->end();
       RecoTrack !=RecoTrackEnd ; ++i, ++RecoTrack) {
    AVTrackStruct trackStruct;
    trackStruct.pt = RecoTrack->pt();
    trackStruct.ptError = RecoTrack->ptError();
    trackStruct.px = RecoTrack->px();
    trackStruct.py = RecoTrack->py();
    trackStruct.pz = RecoTrack->pz();
    trackStruct.eta = RecoTrack->eta();
    trackStruct.phi = RecoTrack->phi();
    trackStruct.chi2 = RecoTrack->chi2();
    trackStruct.normchi2 = RecoTrack->normalizedChi2();
    GlobalPoint gPoint(RecoTrack->vx(), RecoTrack->vy(), RecoTrack->vz());
    double theLocalMagField = magneticField_->inTesla(gPoint).z();
    trackStruct.kappa = -RecoTrack->charge()*0.002998*theLocalMagField/RecoTrack->pt();
    trackStruct.charge = RecoTrack->charge();
    trackStruct.d0 = RecoTrack->d0();
    trackStruct.dz = RecoTrack->dz();
    v_avtrackout.push_back(trackStruct);
  }

}

