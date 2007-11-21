/**
 *  Class: GlobalMuonTrackMatcher
 *
 * 
 *  $Date: 2007/09/28 20:21:10 $
 *  $Revision: 1.2 $
 *
 *  Authors :
 *  \author Chang Liu  - Purdue University
 *  \author Norbert Neumeister - Purdue University
 *  \author Adam Everett - Purdue University
 *
 */

#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonTrackMatcher.h"

//---------------
// C++ Headers --
//---------------


//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Utilities/Timing/interface/TimingReport.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
//#include "TrackingTools/GeomPropagators/interface/StateOnMuonBound.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"

#include "DataFormats/GeometrySurface/interface/TangentPlane.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

using namespace std;
using namespace reco;

//
// constructor
//
GlobalMuonTrackMatcher::GlobalMuonTrackMatcher(const edm::ParameterSet& par, 
                                               const MuonServiceProxy* service) : 
   theService(service) {
  
  theMaxChi2 =  par.getParameter<double>("Chi2Cut");
  theDeltaEta = par.getParameter<double>("DeltaEtaCut");
  theDeltaPhi = par.getParameter<double>("DeltaPhiCut");
  theMinP = par.getParameter<double>("MinP");
  theMinPt = par.getParameter<double>("MinPt");
  
  theOutPropagatorName = par.getParameter<string>("StateOnTrackerBoundOutPropagator");

}


//
// destructor
//
GlobalMuonTrackMatcher::~GlobalMuonTrackMatcher() {

}


/*!
 * Choose the Track from a TrackCollection which has smallest chi2 with
 * a given standalone muon Track.
 */
pair<bool, GlobalMuonTrackMatcher::TrackCand> 
GlobalMuonTrackMatcher::matchOne(const TrackCand& staCand,
				 const vector<TrackCand>& tkTs) const {

  return pair<bool, TrackCand>(false, staCand);
  
}


/*!
 * Choose a vector of Tracks from a TrackCollection that are compatible
 * with a given standalone Track.  The order of checks for compatability are
 * \li matching-chi2 less than MaxChi2
 * \li gloabl position of TSOS on tracker bound
 * \li global momentum direction
 * \see matchChi()
 * \see matchPos()
 * \see matchMomAtIP()
 */
vector<GlobalMuonTrackMatcher::TrackCand>
GlobalMuonTrackMatcher::match(const TrackCand& staCand, 
                              const vector<TrackCand>& tkTs) const {

  const string category = "GlobalMuonTrackMatcher";  
  vector<TrackCand> result; 
  
  if ( tkTs.empty() ) return result;
  
  for (vector<TrackCand>::const_iterator is = tkTs.begin(); is != tkTs.end(); ++is) {

    /*    
    std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairTk 
      = convertToTSOSTk(staCand,*is);
    bool sameSurfaceTk = samePlane(tsosPairTk.first,tsosPairTk.second);

    std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairMu 
      = convertToTSOSMu(staCand,*is);
    bool sameSurfaceMu = samePlane(tsosPairMu.first,tsosPairMu.second);
    */

    std::pair<TrajectoryStateOnSurface, TrajectoryStateOnSurface> tsosPairMuHit
      = convertToTSOSMuHit(staCand,*is);
    bool sameSurfaceMuHit = samePlane(tsosPairMuHit.first,tsosPairMuHit.second);

    if(sameSurfaceMuHit && match_Rpos(tsosPairMuHit.first,tsosPairMuHit.second) < theDeltaEta ) result.push_back(*is); 

  }

  return result;

}



std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
GlobalMuonTrackMatcher::convertToTSOSTk(const TrackCand& staCand,
			    const TrackCand& tkCand) const {
  
  const string category = "GlobalMuonTrackMatcher";
  
  TransientTrack muTT(*staCand.second,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();

  TrajectoryStateOnSurface outerTkTsos;
  if(tkCand.first == 0) {
    //make sure the trackerTrack has enough momentum to reach the muon chambers
    if ( !(tkCand.second->p() < theMinP || tkCand.second->pt() < theMinPt )) {
      TrajectoryStateTransform tsTransform;
      outerTkTsos = tsTransform.outerStateOnSurface(*tkCand.second,*theService->trackingGeometry(),&*theService->magneticField());
    }
  } else {
    const GlobalVector& mom = tkCand.first->firstMeasurement().updatedState().globalMomentum();
    if(!(mom.mag() < theMinP || mom.perp() < theMinPt)) {
      outerTkTsos = (tkCand.first->direction() == alongMomentum) ? tkCand.first->lastMeasurement().updatedState() : tkCand.first->firstMeasurement().updatedState();
    }
  }
  
  if ( !impactMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(impactMuTSOS,outerTkTsos);
  
  // define StateOnTrackerBound objects  
  StateOnTrackerBound fromInside(&*theService->propagator(theOutPropagatorName));
  
  // extrapolate to outer tracker surface
  TrajectoryStateOnSurface tkTsosFromMu = fromInside(impactMuTSOS);
  TrajectoryStateOnSurface tkTsosFromTk = fromInside(outerTkTsos);

    
  if( !samePlane(tkTsosFromMu,tkTsosFromTk)) {
    bool same1, same2;
    //propagate tk to same surface as muon
    TrajectoryStateOnSurface newTkTsosFromTk, newTkTsosFromMu;
    if( tkTsosFromMu.isValid() ) newTkTsosFromTk = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,tkTsosFromMu.surface());
    same1 =  samePlane(newTkTsosFromTk,tkTsosFromMu);
    LogDebug(category) << "Propagating to same tracker surface (Mu):" << same1;
    if( !same1 ) {
      if( tkTsosFromTk.isValid() ) newTkTsosFromMu = theService->propagator(theOutPropagatorName)->propagate(impactMuTSOS,tkTsosFromTk.surface());
      same2 =  samePlane(newTkTsosFromMu,tkTsosFromTk);
      LogDebug(category) << "Propagating to same tracker surface (Tk):" << same2;
    }
    if(same1) tkTsosFromTk = newTkTsosFromTk;
    else if(same2) tkTsosFromMu = newTkTsosFromMu;
    else  {
      LogDebug(category) << "Could not propagate Muon and Tracker track to the same tracker bound!";
      return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(impactMuTSOS, outerTkTsos);
    }
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(tkTsosFromMu, tkTsosFromTk);
}

/*
std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
GlobalMuonTrackMatcher::convertToTSOSMu(const TrackCand& staCand,
			    const TrackCand& tkCand) const {
  
  const string category = "GlobalMuonTrackMatcher";
  
  TransientTrack muTT(*staCand.second,&*theService->magneticField(),theService->trackingGeometry());
  //TrajectoryStateOnSurface impactMuTSOS = muTT.impactPointState();
  TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();

  TrajectoryStateOnSurface outerTkTsos;
  if(tkCand.first == 0) {
    //make sure the trackerTrack has enough momentum to reach the muon chambers
    if ( !(tkCand.second->p() < theMinP || tkCand.second->pt() < theMinPt )) {
      TrajectoryStateTransform tsTransform;
      outerTkTsos = tsTransform.outerStateOnSurface(*tkCand.second,*theService->trackingGeometry(),&*theService->magneticField());
    }
  } else {
    const GlobalVector& mom = tkCand.first->firstMeasurement().updatedState().globalMomentum();
    if(!(mom.mag() < theMinP || mom.perp() < theMinPt)) {
      outerTkTsos = (tkCand.first->direction() == alongMomentum) ? tkCand.first->lastMeasurement().updatedState() : tkCand.first->firstMeasurement().updatedState();
    }
  }


  if ( !innerMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS,outerTkTsos);
  
  // define StateOnTrackerBound objects  
  StateOnMuonBound fromInside(&*theService->propagator(theOutPropagatorName));
  StateOnMuonBound fromOutside(&*theService->propagator(theOutPropagatorName));
  
  // extrapolate to outer tracker surface
  TrajectoryStateOnSurface muTsosFromMu = fromOutside(innerMuTSOS);
  TrajectoryStateOnSurface muTsosFromTk = fromInside(outerTkTsos);

    
  if( !samePlane(muTsosFromMu,muTsosFromTk)) {
    bool same1, same2;
    //propagate tk to same surface as muon
    TrajectoryStateOnSurface newMuTsosFromTk, newMuTsosFromMu;
    if( muTsosFromMu.isValid() ) newMuTsosFromTk = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,muTsosFromMu.surface());
    same1 =  samePlane(newMuTsosFromTk,muTsosFromMu);
    LogDebug(category) << "Propagating to same muon surface (Mu):" << same1;
    if( !same1 ) {
      if( muTsosFromTk.isValid() ) newMuTsosFromMu = theService->propagator(theOutPropagatorName)->propagate(innerMuTSOS,muTsosFromTk.surface());
      same2 =  samePlane(newMuTsosFromMu,muTsosFromTk);
      LogDebug(category) << "Propagating to same muon surface (Tk):" << same2;
    }
    if(same1) muTsosFromTk = newMuTsosFromTk;
    else if(same2) muTsosFromMu = newMuTsosFromMu;
    else {
      LogDebug(category) << "Could not propagate Muon and Tracker track to the same muon bound!";
      return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, outerTkTsos);
  }
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(muTsosFromMu, muTsosFromTk);
}
*/


std::pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>
GlobalMuonTrackMatcher::convertToTSOSMuHit(const TrackCand& staCand,
			       const TrackCand& tkCand) const {
  
  const string category = "GlobalMuonTrackMatcher";
  
  TransientTrack muTT(*staCand.second,&*theService->magneticField(),theService->trackingGeometry());
  TrajectoryStateOnSurface innerMuTSOS = muTT.innermostMeasurementState();


  TrajectoryStateOnSurface outerTkTsos;
  if(tkCand.first == 0) {
    //make sure the trackerTrack has enough momentum to reach the muon chambers
    if ( !(tkCand.second->p() < theMinP || tkCand.second->pt() < theMinPt )) {
      TrajectoryStateTransform tsTransform;
      outerTkTsos = tsTransform.outerStateOnSurface(*tkCand.second,*theService->trackingGeometry(),&*theService->magneticField());
    }
  } else {
    const GlobalVector& mom = tkCand.first->firstMeasurement().updatedState().globalMomentum();
    if(!(mom.mag() < theMinP || mom.perp() < theMinPt)) {
      outerTkTsos = (tkCand.first->direction() == alongMomentum) ? tkCand.first->lastMeasurement().updatedState() : tkCand.first->firstMeasurement().updatedState();
    }
  }


  if ( !innerMuTSOS.isValid() || !outerTkTsos.isValid() ) return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS,outerTkTsos);
  
  TrajectoryStateOnSurface tkAtMu = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,theService->trackingGeometry()->idToDet( DetId(staCand.second->innerDetId()) )->surface());
  //TrajectoryStateOnSurface tkAtMu = theService->propagator(theOutPropagatorName)->propagate(outerTkTsos,innerMuTSOS.surface());
  
  
  if( !samePlane(innerMuTSOS,tkAtMu)) {
    LogDebug(category) << "Could not propagate Muon and Tracker track to the same muon hit surface!";
    return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, outerTkTsos);    
  }
  
  return pair<TrajectoryStateOnSurface,TrajectoryStateOnSurface>(innerMuTSOS, tkAtMu);
}


bool
GlobalMuonTrackMatcher::samePlane(const TrajectoryStateOnSurface& tsos1,
				  const TrajectoryStateOnSurface& tsos2) const
{
  if( !tsos1.isValid() || !tsos2.isValid()) return false;
  const string category = "GlobalMuonTrackMatcher";

  const float maxtilt = 0.999;
  const float maxdist = 0.01; // in cm

  ReferenceCountingPointer<TangentPlane> p1(tsos1.surface().tangentPlane(tsos1.localPosition()));
  ReferenceCountingPointer<TangentPlane> p2(tsos2.surface().tangentPlane(tsos2.localPosition()));

  bool returnValue =  ( (fabs(p1->normalVector().dot(p2->normalVector())) > maxtilt) || (fabs((p1->toLocal(p2->position())).z()) < maxdist) ) ? true : false;

  return returnValue; 
  
}

double 
GlobalMuonTrackMatcher::matchChiAtSurface(const TrajectoryStateOnSurface& tsos1, 
			      const TrajectoryStateOnSurface& tsos2) const {
  
  const string category = "GlobalMuonTrackMatcher";
  
  if ( !tsos1.isValid() || !tsos2.isValid() ) return -1.;

  AlgebraicVector5 v(tsos1.localParameters().vector() - tsos2.localParameters().vector());
  AlgebraicSymMatrix55 m(tsos1.localError().matrix() + tsos2.localError().matrix());
  //LogDebug(category) << "vector v " << v;

  int ierr = ! m.Invert();

  if (ierr != 0) edm::LogInfo(category) << "Error inversing covariance matrix";

  double est = ROOT::Math::Similarity(v,m);

  //LogDebug(category) << "Chi2 " << est;

/*
  GlobalVector x = tsos1.globalParameters().position() - tsos2.globalParameters().position();
  AlgebraicVector v1(3); v1[0] = x.x(); v1[1] = x.y(); v1[2] = x.z();
  AlgebraicSymMatrix m1(tsos1.cartesianError().position().matrix() + tsos2.cartesianError().position().matrix());
  m1.invert(ierr);
  double est1 = m1.similarity(v1);
*/

  return est;

}

double
GlobalMuonTrackMatcher::match_R_IP(const reco::Track& staTrack, const reco::Track& tkTrack) const {
  return (deltaR<double>(staTrack.eta(),staTrack.phi(),
			 tkTrack.eta(),tkTrack.phi()));
}

double
GlobalMuonTrackMatcher::match_Rmom(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (deltaR<double>(sta.globalMomentum().eta(),sta.globalMomentum().phi(),
			 tk.globalMomentum().eta(),tk.globalMomentum().phi()));
}

double
GlobalMuonTrackMatcher::match_Rpos(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (deltaR<double>(sta.globalPosition().eta(),sta.globalPosition().phi(),
			 tk.globalPosition().eta(),tk.globalPosition().phi()));
}

double
GlobalMuonTrackMatcher::match_D(const TrajectoryStateOnSurface& sta, const TrajectoryStateOnSurface& tk) const {
  return (sta.globalPosition() - tk.globalPosition()).mag();
}
