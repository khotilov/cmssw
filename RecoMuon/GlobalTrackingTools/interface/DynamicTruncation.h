#ifndef RecoMuon_GlobalTrackingTools_DynamicTruncation_h
#define RecoMuon_GlobalTrackingTools_DynamicTruncation_h

/**
 *  Class: DynamicTruncation
 *
 *  Description:
 *  class for the dynamical stop of the KF according to the
 *  compatibility degree between the extrapolated track
 *  state and the reconstructed segment in the muon chambers
 *
 *  $Date: 2010/06/27 17:32:56 $
 *  $Revision: 1.2 $
 *
 *  Authors :
 *  D. Pagano & G. Bruno - UCL Louvain
 *
 **/

#include <memory>
#include "RecoMuon/GlobalTrackingTools/interface/DirectTrackerNavigation.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "RecoMuon/GlobalTrackingTools/interface/StateSegmentMatcher.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "RecoMuon/Navigation/interface/MuonNavigableLayer.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"
#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"


class DynamicTruncation {
  
 public:

  typedef TransientTrackingRecHit::ConstRecHitPointer ConstRecHitPointer;
  typedef TransientTrackingRecHit::ConstRecHitContainer ConstRecHitContainer;

  DynamicTruncation(const edm::Event&, const MuonServiceProxy&);

  ~DynamicTruncation();

  // Set thresholds for DYT
  void setThr(const std::vector<int>&);
  
  // Return the vector with the tracker plus the selected muon hits
  TransientTrackingRecHit::ConstRecHitContainer filter(const Trajectory&);
 
 private:

  void                 setThrs(float);
  void                 compatibleDets(TrajectoryStateOnSurface&, std::map<int, std::vector<DetId> >&);
  void                 filteringAlgo(std::map<int, std::vector<DetId> >&);
  double               getBest(std::vector<CSCSegment>&, TrajectoryStateOnSurface&, CSCSegment&); 
  double               getBest(std::vector<DTRecSegment4D>&, TrajectoryStateOnSurface&, DTRecSegment4D&); 
  void                 update(TrajectoryStateOnSurface&, ConstRecHitPointer);
  void                 updateWithDThits(ConstRecHitContainer&);
  void                 updateWithCSChits(ConstRecHitContainer&);
  ConstRecHitContainer returnGlobal(const Trajectory&);
  ConstRecHitContainer sort(ConstRecHitContainer&);
  
  ConstRecHitContainer result;
  
  int DTThr;
  int CSCThr;

  std::vector<int> DYTthrs;

  edm::ESHandle<Propagator> propagator;
  edm::ESHandle<Propagator> propagatorCompatibleDet;
  edm::ESHandle<GlobalTrackingGeometry> theG;
  edm::ESHandle<CSCGeometry> cscGeom;
  edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;
  edm::ESHandle<TrajectoryStateUpdator> updatorHandle;
  edm::ESHandle<MuonDetLayerGeometry> navMuon;
  DirectMuonNavigation *navigation;
  edm::ESHandle<MagneticField> magfield;
  const edm::Event* theEvent;
  const edm::EventSetup* theSetup;

  TrajectoryStateOnSurface currentState;
};

#endif


