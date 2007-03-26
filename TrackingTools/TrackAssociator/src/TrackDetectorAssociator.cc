// -*- C++ -*-
//
// Package:    TrackAssociator
// Class:      TrackDetectorAssociator
// 
/*

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Fri Apr 21 10:59:41 PDT 2006
// $Id: TrackDetectorAssociator.cc,v 1.7 2007/03/09 14:08:15 dmytro Exp $
//
//

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

// calorimeter info
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "Utilities/Timing/interface/TimingReport.h"
#include <stack>
#include <set>


#include "TrackingTools/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TimerStack.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include "HepPDT/ParticleID.hh"
//
// class declaration
//

using namespace reco;

TrackDetectorAssociator::TrackDetectorAssociator() 
{
   ivProp_ = 0;
   defProp_ = 0;
   useDefaultPropagator_ = false;
}

TrackDetectorAssociator::~TrackDetectorAssociator()
{
   if (defProp_) delete defProp_;
}

void TrackDetectorAssociator::setPropagator( Propagator* ptr)
{
   ivProp_ = ptr;
   cachedTrajectory_.setPropagator(ivProp_);
}

void TrackDetectorAssociator::useDefaultPropagator()
{
   useDefaultPropagator_ = true;
}


void TrackDetectorAssociator::init( const edm::EventSetup& iSetup )
{
   // access the calorimeter geometry
   iSetup.get<IdealGeometryRecord>().get(theCaloGeometry_);
   if (!theCaloGeometry_.isValid()) 
     throw cms::Exception("FatalError") << "Unable to find IdealGeometryRecord in event!\n";
   
   // get the tracking Geometry
   iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry_);
   if (!theTrackingGeometry_.isValid()) 
     throw cms::Exception("FatalError") << "Unable to find GlobalTrackingGeometryRecord in event!\n";
   
   if (useDefaultPropagator_ && ! defProp_ ) {
      // setup propagator
      edm::ESHandle<MagneticField> bField;
      iSetup.get<IdealMagneticFieldRecord>().get(bField);
      
      SteppingHelixPropagator* prop  = new SteppingHelixPropagator(&*bField,anyDirection);
      prop->setMaterialMode(false);
      prop->applyRadX0Correction(true);
      // prop->setDebug(true); // tmp
      defProp_ = prop;
      setPropagator(defProp_);
   }
   
	
}

TrackDetMatchInfo TrackDetectorAssociator::associate( const edm::Event& iEvent,
					      const edm::EventSetup& iSetup,
					      const FreeTrajectoryState& fts,
					      const AssociatorParameters& parameters )
{
   TrackDetMatchInfo info;
   TimerStack timers;
   SteppingHelixStateInfo trackOrigin(fts);
   
   init( iSetup );
   
   // get track trajectory
   // timers.push("TrackDetectorAssociator::fillEcal::propagation");
   // ECAL points (EB+EE)
   // If the phi angle between a track entrance and exit points is more
   // than 2 crystals, it is possible that the track will cross 3 crystals
   // and therefore one has to check at least 3 points along the track
   // trajectory inside ECAL. In order to have a chance to cross 4 crystalls
   // in the barrel, a track should have P_t as low as 3 GeV or smaller
   // If it's necessary, number of points along trajectory can be increased
   cachedTrajectory_.reset_trajectory();
   cachedTrajectory_.propagateAll(trackOrigin);
   cachedTrajectory_.getEcalTrajectory();
   cachedTrajectory_.getHcalTrajectory();
   cachedTrajectory_.getHOTrajectory();

   info.trkGlobPosAtEcal = getPoint( cachedTrajectory_.getStateAtEcal().position() );
   info.trkGlobPosAtHcal = getPoint( cachedTrajectory_.getStateAtHcal().position() );
   info.trkGlobPosAtHO = getPoint( cachedTrajectory_.getStateAtHO().position() );

   if (parameters.useEcal) fillEcal( iEvent, info, parameters);
   if (parameters.useCalo) fillCaloTowers( iEvent, info, parameters);
   if (parameters.useHcal) fillHcal( iEvent, info, parameters);
   if (parameters.useHO)   fillHO( iEvent, info, parameters);
   if (parameters.useMuon) fillMuon( iEvent, info, parameters);

   return info;
}

void TrackDetectorAssociator::fillEcal( const edm::Event& iEvent,
				TrackDetMatchInfo& info,
				const AssociatorParameters& parameters)
{
   TimerStack timers;
   timers.push("TrackDetectorAssociator::fillEcal");
   
   const std::vector<SteppingHelixStateInfo>& trajectoryStates = cachedTrajectory_.getEcalTrajectory();
   std::vector<GlobalPoint> trajectory;
   for(std::vector<SteppingHelixStateInfo>::const_iterator itr = trajectoryStates.begin();
       itr != trajectoryStates.end(); itr++) trajectory.push_back(itr->position());
   
   ecalDetIdAssociator_.setGeometry(&*theCaloGeometry_);
   
   if(trajectory.empty()) {
      LogTrace("TrackAssociator") << "ECAL track trajectory is empty; moving on\n";
      info.isGoodEcal = 0;
      return;
   }
   info.isGoodEcal = 1;

   // Find ECAL crystals
   timers.pop_and_push("TrackDetectorAssociator::fillEcal::access::EcalBarrel");
   edm::Handle<EBRecHitCollection> EBRecHits;
   iEvent.getByLabel (theEBRecHitCollectionLabel, EBRecHits);
   if (!EBRecHits.isValid()) throw cms::Exception("FatalError") << "Unable to find EBRecHitCollection in the event!\n";

   timers.pop_and_push("TrackDetectorAssociator::fillEcal::access::EcalEndcaps");
   edm::Handle<EERecHitCollection> EERecHits;
   iEvent.getByLabel (theEERecHitCollectionLabel, EERecHits);
   if (!EERecHits.isValid()) throw cms::Exception("FatalError") << "Unable to find EERecHitCollection in event!\n";

   timers.pop_and_push("TrackDetectorAssociator::fillEcal::matching");
   std::set<DetId> ecalIdsInRegion = ecalDetIdAssociator_.getDetIdsCloseToAPoint(trajectory[0],parameters.dREcalPreselection);
   LogTrace("TrackAssociator") << "ECAL hits in the region: " << ecalIdsInRegion.size();
   std::set<DetId> ecalIdsInACone =  ecalDetIdAssociator_.getDetIdsInACone(ecalIdsInRegion, trajectory, parameters.dREcal);
   LogTrace("TrackAssociator") << "ECAL hits in the cone: " << ecalIdsInACone.size();
   std::vector<DetId> crossedEcalIds =  ecalDetIdAssociator_.getCrossedDetIdsOrdered(ecalIdsInRegion, trajectory);
   LogTrace("TrackAssociator") << "ECAL crossed hits " << crossedEcalIds.size();
   
   info.crossedEcalIds = crossedEcalIds;
   
   // add EcalRecHits
   timers.pop_and_push("TrackDetectorAssociator::fillEcal::addEcalRecHits");
   for(std::vector<DetId>::const_iterator itr=crossedEcalIds.begin(); itr!=crossedEcalIds.end();itr++)
   {
      std::vector<EcalRecHit>::const_iterator ebHit = (*EBRecHits).find(*itr);
      std::vector<EcalRecHit>::const_iterator eeHit = (*EERecHits).find(*itr);
      if(ebHit != (*EBRecHits).end()) 
         info.crossedEcalRecHits.push_back(*ebHit);
      else if(eeHit != (*EERecHits).end()) 
         info.crossedEcalRecHits.push_back(*eeHit);
      else  
         LogTrace("TrackAssociator") << "Crossed EcalRecHit is not found for DetId: " << itr->rawId();
   }
   for(std::set<DetId>::const_iterator itr=ecalIdsInACone.begin(); itr!=ecalIdsInACone.end();itr++)
   {
      std::vector<EcalRecHit>::const_iterator ebHit = (*EBRecHits).find(*itr);
      std::vector<EcalRecHit>::const_iterator eeHit = (*EERecHits).find(*itr);
      if(ebHit != (*EBRecHits).end()) 
         info.ecalRecHits.push_back(*ebHit);
      else if(eeHit != (*EERecHits).end()) 
         info.ecalRecHits.push_back(*eeHit);
      else 
         LogTrace("TrackAssociator") << "EcalRecHit from the cone is not found for DetId: " << itr->rawId();
   }
}

void TrackDetectorAssociator::fillCaloTowers( const edm::Event& iEvent,
				      TrackDetMatchInfo& info,
				      const AssociatorParameters& parameters)
{
   TimerStack timers;
   timers.push("TrackDetectorAssociator::fillCaloTowers");

   caloDetIdAssociator_.setGeometry(&*theCaloGeometry_);
   
   // use ECAL and HCAL trajectories to match a tower. (HO isn't used for matching).
   std::vector<GlobalPoint> trajectory;
   const std::vector<SteppingHelixStateInfo>& ecalTrajectoryStates = cachedTrajectory_.getEcalTrajectory();
   const std::vector<SteppingHelixStateInfo>& hcalTrajectoryStates = cachedTrajectory_.getHcalTrajectory();
   for(std::vector<SteppingHelixStateInfo>::const_iterator itr = ecalTrajectoryStates.begin();
       itr != ecalTrajectoryStates.end(); itr++) trajectory.push_back(itr->position());
   for(std::vector<SteppingHelixStateInfo>::const_iterator itr = hcalTrajectoryStates.begin();
       itr != hcalTrajectoryStates.end(); itr++) trajectory.push_back(itr->position());
   
   if(trajectory.empty()) {
      LogTrace("TrackAssociator") << "HCAL trajectory is empty; moving on\n";
      info.isGoodCalo = 0;
      return;
   }
   info.isGoodCalo = 1;
   
   // find crossed CaloTowers
   timers.pop_and_push("TrackDetectorAssociator::fillCaloTowers::access::CaloTowers");
   edm::Handle<CaloTowerCollection> caloTowers;

   iEvent.getByLabel (theCaloTowerCollectionLabel, caloTowers);
   if (!caloTowers.isValid())  throw cms::Exception("FatalError") << "Unable to find CaloTowers in event!\n";
   
   timers.push("TrackDetectorAssociator::fillCaloTowers::matching");
   std::set<DetId> caloTowerIdsInRegion = caloDetIdAssociator_.getDetIdsCloseToAPoint(trajectory[0],parameters.dRHcalPreselection);
   LogTrace("TrackAssociator") << "Towers in the region: " << caloTowerIdsInRegion.size();
   std::set<DetId> caloTowerIdsInACone = caloDetIdAssociator_.getDetIdsInACone(caloTowerIdsInRegion, trajectory, parameters.dRHcal);
   LogTrace("TrackAssociator") << "Towers in the cone: " << caloTowerIdsInACone.size();
   std::vector<DetId> crossedCaloTowerIds = caloDetIdAssociator_.getCrossedDetIdsOrdered(caloTowerIdsInRegion, trajectory);
   LogTrace("TrackAssociator") << "Towers crossed: " << crossedCaloTowerIds.size();
   
   info.crossedTowerIds = crossedCaloTowerIds;
   
   // add CaloTowers
   timers.push("TrackDetectorAssociator::fillCaloTowers::addCaloTowers");
   for(std::vector<DetId>::const_iterator itr=crossedCaloTowerIds.begin(); itr!=crossedCaloTowerIds.end();itr++)
     {
	CaloTowerCollection::const_iterator tower = (*caloTowers).find(*itr);
	if(tower != (*caloTowers).end()) 
	  info.crossedTowers.push_back(*tower);
	else
	  LogTrace("TrackAssociator") << "Crossed CaloTower is not found for DetId: " << (*itr).rawId();
     }

   for(std::set<DetId>::const_iterator itr=caloTowerIdsInACone.begin(); itr!=caloTowerIdsInACone.end();itr++)
     {
	CaloTowerCollection::const_iterator tower = (*caloTowers).find(*itr);
	if(tower != (*caloTowers).end()) 
	  info.towers.push_back(*tower);
	else 
	  LogTrace("TrackAssociator") << "CaloTower from the cone is not found for DetId: " << (*itr).rawId();
     }
   
}

void TrackDetectorAssociator::fillHcal( const edm::Event& iEvent,
				TrackDetMatchInfo& info,
				const AssociatorParameters& parameters)
{
   TimerStack timers;
   timers.push("TrackDetectorAssociator::fillHcals");

   hcalDetIdAssociator_.setGeometry(&*theCaloGeometry_);
   
   const std::vector<SteppingHelixStateInfo>& trajectoryStates = cachedTrajectory_.getHcalTrajectory();
   std::vector<GlobalPoint> trajectory;
   for(std::vector<SteppingHelixStateInfo>::const_iterator itr = trajectoryStates.begin();
       itr != trajectoryStates.end(); itr++) trajectory.push_back(itr->position());

   if(trajectory.empty()) {
      LogTrace("TrackAssociator") << "HCAL trajectory is empty; moving on\n";
      info.isGoodHcal = 0;
      return;
   }
   info.isGoodHcal = 1;
   
   // find crossed Hcals
   timers.pop_and_push("TrackDetectorAssociator::fillHcal::access::Hcal");
   edm::Handle<HBHERecHitCollection> collection;

   iEvent.getByLabel (theHBHERecHitCollectionLabel, collection);
   if ( ! collection.isValid() ) throw cms::Exception("FatalError") << "Unable to find HBHERecHits in event!\n";
   
   timers.push("TrackDetectorAssociator::fillHcal::matching");
   std::set<DetId> idsInRegion = hcalDetIdAssociator_.getDetIdsCloseToAPoint(trajectory[0],parameters.dRHcalPreselection);
   LogTrace("TrackAssociator") << "HCAL hits in the region: " << idsInRegion.size() << "\n" << DetIdInfo::info(idsInRegion);
   std::set<DetId> idsInACone = hcalDetIdAssociator_.getDetIdsInACone(idsInRegion, trajectory, parameters.dRHcal);
   LogTrace("TrackAssociator") << "HCAL hits in the cone: " << idsInACone.size() << "\n" << DetIdInfo::info(idsInACone);
   std::vector<DetId> crossedIds = hcalDetIdAssociator_.getCrossedDetIdsOrdered(idsInRegion, trajectory);
   LogTrace("TrackAssociator") << "HCAL hits crossed: " << crossedIds.size() << "\n" << DetIdInfo::info(crossedIds);
   
   info.crossedHcalIds = crossedIds;
   
   // add Hcal
   timers.push("TrackDetectorAssociator::fillHcal::addHcal");
   for(std::vector<DetId>::const_iterator itr=crossedIds.begin(); itr!=crossedIds.end();itr++)
     {
	HBHERecHitCollection::const_iterator hit = (*collection).find(*itr);
	if( hit != (*collection).end() ) 
	  info.crossedHcalRecHits.push_back(*hit);
	else
	  LogTrace("TrackAssociator") << "Crossed HBHERecHit is not found for DetId: " << itr->rawId();
     }

   for(std::set<DetId>::const_iterator itr=idsInACone.begin(); itr!=idsInACone.end();itr++)
     {
	HBHERecHitCollection::const_iterator hit = (*collection).find(*itr);
	if( hit != (*collection).end() ) 
	  info.hcalRecHits.push_back(*hit);
	else 
	  LogTrace("TrackAssociator") << "HBHERecHit from the cone is not found for DetId: " << itr->rawId();
     }
}

void TrackDetectorAssociator::fillHO( const edm::Event& iEvent,
			      TrackDetMatchInfo& info,
			      const AssociatorParameters& parameters)
{
   TimerStack timers;
   timers.push("TrackDetectorAssociator::fillHO");

   hoDetIdAssociator_.setGeometry(&*theCaloGeometry_);
   
   const std::vector<SteppingHelixStateInfo>& trajectoryStates = cachedTrajectory_.getHOTrajectory();
   std::vector<GlobalPoint> trajectory;
   for(std::vector<SteppingHelixStateInfo>::const_iterator itr = trajectoryStates.begin();
       itr != trajectoryStates.end(); itr++) trajectory.push_back(itr->position());

   if(trajectory.empty()) {
      LogTrace("TrackAssociator") << "HO trajectory is empty; moving on\n";
      info.isGoodHO = 0;
      return;
   }
   info.isGoodHO = 1;
   
   // find crossed HOs
   timers.pop_and_push("TrackDetectorAssociator::fillHO::access::HO");
   edm::Handle<HORecHitCollection> collection;

   iEvent.getByLabel (theHORecHitCollectionLabel, collection);
   if ( ! collection.isValid() ) throw cms::Exception("FatalError") << "Unable to find HORecHits in event!\n";
   
   timers.push("TrackDetectorAssociator::fillHO::matching");
   std::set<DetId> idsInRegion = hoDetIdAssociator_.getDetIdsCloseToAPoint(trajectory[0],parameters.dRHcalPreselection);
   LogTrace("TrackAssociator") << "idsInRegion.size(): " << idsInRegion.size();
   std::set<DetId> idsInACone = hoDetIdAssociator_.getDetIdsInACone(idsInRegion, trajectory, parameters.dRHcal);
   LogTrace("TrackAssociator") << "idsInACone.size(): " << idsInACone.size();
   std::vector<DetId> crossedIds = hoDetIdAssociator_.getCrossedDetIdsOrdered(idsInRegion, trajectory);
   LogTrace("TrackAssociator") << "crossedIds.size(): " << crossedIds.size();
   
   info.crossedHOIds = crossedIds;
   
   // add HO
   timers.push("TrackDetectorAssociator::fillHO::addHO");
   for(std::vector<DetId>::const_iterator itr=crossedIds.begin(); itr!=crossedIds.end();itr++)
     {
	HORecHitCollection::const_iterator hit = (*collection).find(*itr);
	if( hit != (*collection).end() ) 
	  info.crossedHORecHits.push_back(*hit);
	else
	  LogTrace("TrackAssociator") << "Crossed HORecHit is not found for DetId: " << itr->rawId();
     }

   for(std::set<DetId>::const_iterator itr=idsInACone.begin(); itr!=idsInACone.end();itr++)
     {
	HORecHitCollection::const_iterator hit = (*collection).find(*itr);
	if( hit != (*collection).end() ) 
	  info.hoRecHits.push_back(*hit);
	else 
	  LogTrace("TrackAssociator") << "HORecHit from the cone is not found for DetId: " << itr->rawId();
     }
}

FreeTrajectoryState TrackDetectorAssociator::getFreeTrajectoryState( const edm::EventSetup& iSetup, 
							     const SimTrack& track, 
							     const SimVertex& vertex )
{
   edm::ESHandle<MagneticField> bField;
   iSetup.get<IdealMagneticFieldRecord>().get(bField);
   
   GlobalVector vector( track.momentum().x(), track.momentum().y(), track.momentum().z() );
   GlobalPoint point( vertex.position().x(), vertex.position().y(), vertex.position().z() );

   HepPDT::ParticleID id(track.type());
   int charge = id.threeCharge() < 0 ? -1 : 1;

   GlobalTrajectoryParameters tPars(point, vector, charge, &*bField);
   
   HepSymMatrix covT(6,1); covT *= 1e-6; // initialize to sigma=1e-3
   CartesianTrajectoryError tCov(covT);
   
   return FreeTrajectoryState(tPars, tCov);
}


FreeTrajectoryState TrackDetectorAssociator::getFreeTrajectoryState( const edm::EventSetup& iSetup,
							     const reco::Track& track )
{
   edm::ESHandle<MagneticField> bField;
   iSetup.get<IdealMagneticFieldRecord>().get(bField);
   
   GlobalVector vector( track.momentum().x(), track.momentum().y(), track.momentum().z() );

   GlobalPoint point( track.vertex().x(), track.vertex().y(),  track.vertex().z() );

   GlobalTrajectoryParameters tPars(point, vector, track.charge(), &*bField);
   
   // FIX THIS !!!
   // need to convert from perigee to global or helix (curvilinear) frame
   // for now just an arbitrary matrix.
   HepSymMatrix covT(6,1); covT *= 1e-6; // initialize to sigma=1e-3
   CartesianTrajectoryError tCov(covT);
   
   return FreeTrajectoryState(tPars, tCov);
}

void TrackDetectorAssociator::getMuonChamberMatches(std::vector<MuonChamberMatch>& matches,
					    const float dRMuonPreselection,
					    const float maxDistanceX,
					    const float maxDistanceY)
{
   // Strategy:
   //    Propagate through the whole detector, estimate change in eta and phi 
   //    along the trajectory, add this to dRMuon and find DetIds around this 
   //    direction using the map. Then propagate fast to each surface and apply 
   //    final matching criteria.

   // TimerStack timers(TimerStack::Disableable);
   // timers.push("MuonDetIdAssociator::getTrajectoryInMuonDetector");
   // timers.push("MuonDetIdAssociator::getTrajectoryInMuonDetector::propagation",TimerStack::FastMonitoring);
   // timers.pop();
   // get the direction first
   SteppingHelixStateInfo trajectoryPoint = cachedTrajectory_.getStateAtHcal();
   if (! trajectoryPoint.isValid() ) {
      LogTrace("TrackAssociator") << 
	"trajectory position at HCAL is not valid. Assume the track cannot reach muon detectors and skip it";
      return;
   }

   GlobalVector direction = trajectoryPoint.momentum().unit();
   LogTrace("TrackAssociator") << "muon direction: " << direction << "\n\t and corresponding point: " <<
     trajectoryPoint.position() <<"\n";
   
   float dEta = cachedTrajectory_.trajectoryDeltaEta();
   float dPhi = cachedTrajectory_.trajectoryDeltaPhi();
   float lookUpCone = ( dEta > dPhi ? dEta : dPhi ) + dRMuonPreselection;
   LogTrace("TrackAssociator") << "dEta, dPhi, lookUpCone" << dEta << ", " << dPhi << ", " << lookUpCone;
   
   // and find chamber DetIds

   // timers.push("MuonDetIdAssociator::getTrajectoryInMuonDetector::getDetIdsCloseToAPoint",TimerStack::FastMonitoring);
   std::set<DetId> muonIdsInRegion = muonDetIdAssociator_.getDetIdsCloseToAPoint(trajectoryPoint.position(), lookUpCone);
   // timers.pop_and_push("MuonDetIdAssociator::getTrajectoryInMuonDetector::matching",TimerStack::FastMonitoring);
   LogTrace("TrackAssociator") << "Number of chambers to check: " << muonIdsInRegion.size();
	
   for(std::set<DetId>::const_iterator detId = muonIdsInRegion.begin(); detId != muonIdsInRegion.end(); detId++)
     {
	const GeomDet* geomDet = muonDetIdAssociator_.getGeomDet(*detId);
	// timers.push("MuonDetIdAssociator::getTrajectoryInMuonDetector::matching::localPropagation",TimerStack::FastMonitoring);
	TrajectoryStateOnSurface stateOnSurface = cachedTrajectory_.propagate( &geomDet->surface() );
	if (! stateOnSurface.isValid()) {
	   LogTrace("TrackAssociator") << "Failed to propagate the track; moving on\n\t"<<
	     detId->rawId() << " not crossed\n"; ;
	   continue;
	}
	// timers.pop_and_push("MuonDetIdAssociator::getTrajectoryInMuonDetector::matching::geometryAccess",TimerStack::FastMonitoring);
	LocalPoint localPoint = geomDet->surface().toLocal(stateOnSurface.freeState()->position());
	float distanceX = fabs(localPoint.x()) - geomDet->surface().bounds().width()/2;
	float distanceY = fabs(localPoint.y()) - geomDet->surface().bounds().length()/2;
	// timers.pop_and_push("MuonDetIdAssociator::getTrajectoryInMuonDetector::matching::checking",TimerStack::FastMonitoring);
	if (distanceX < maxDistanceX && distanceY < maxDistanceY) {
	   LogTrace("TrackAssociator") << "found a match, DetId: " << detId->rawId();
	   MuonChamberMatch match;
	   match.tState = stateOnSurface;
	   match.localDistanceX = distanceX;
	   match.localDistanceY = distanceY;
	   match.id = *detId;
	   matches.push_back(match);
	}
	//timers.pop();
     }
   //timers.pop();
   
}


void TrackDetectorAssociator::fillMuon( const edm::Event& iEvent,
					TrackDetMatchInfo& info,
					const AssociatorParameters& parameters)
{
   TimerStack timers;
   timers.push("TrackDetectorAssociator::fillMuon");

   muonDetIdAssociator_.setGeometry(&*theTrackingGeometry_);

   // Get the segments from the event
   timers.push("TrackDetectorAssociator::fillMuon::access");
   edm::Handle<DTRecSegment4DCollection> dtSegments;
   iEvent.getByLabel (theDTRecSegment4DCollectionLabel, dtSegments);
   if (! dtSegments.isValid()) 
     throw cms::Exception("FatalError") << "Unable to find DTRecSegment4DCollection in event!\n";
   
   edm::Handle<CSCSegmentCollection> cscSegments;
   iEvent.getByLabel (theCSCSegmentCollectionLabel, cscSegments);
   if (! cscSegments.isValid()) 
     throw cms::Exception("FatalError") << "Unable to find CSCSegmentCollection in event!\n";

   ///// get a set of DetId's in a given direction
   
   // check the map of available segments
   // if there is no segments in a given direction at all,
   // then there is no point to fly there.
   // 
   // MISSING
   // Possible solution: quick search for presence of segments 
   // for the set of DetIds

   timers.pop_and_push("TrackDetectorAssociator::fillMuon::matchChembers");
   
   // get a set of matches corresponding to muon chambers
   std::vector<MuonChamberMatch> matchedChambers;
   getMuonChamberMatches(matchedChambers, parameters.dRMuonPreselection, parameters.muonMaxDistanceX, parameters.muonMaxDistanceY);
   LogTrace("TrackAssociator") << "Chambers matched: " << matchedChambers.size() << "\n";
   
   // Iterate over all chamber matches and fill segment matching 
   // info if it's available
   timers.pop_and_push("TrackDetectorAssociator::fillMuon::findSemgents");
   for(std::vector<MuonChamberMatch>::iterator matchedChamber = matchedChambers.begin(); 
       matchedChamber != matchedChambers.end(); matchedChamber++)
     {
	const GeomDet* geomDet = muonDetIdAssociator_.getGeomDet((*matchedChamber).id);
	// DT chamber
	if(const DTChamber* chamber = dynamic_cast<const DTChamber*>(geomDet) ) {
	   // Get the range for the corresponding segments
	   DTRecSegment4DCollection::range  range = dtSegments->get(chamber->id());
	   // Loop over the segments of this chamber
	   for (DTRecSegment4DCollection::const_iterator segment = range.first; segment!=range.second; segment++)
	     addMuonSegmentMatch(*matchedChamber, &(*segment), parameters);
	}else{
	   // CSC Chamber
	   if(const CSCChamber* chamber = dynamic_cast<const CSCChamber*>(geomDet) ) {
	      // Get the range for the corresponding segments
	      CSCSegmentCollection::range  range = cscSegments->get(chamber->id());
	      // Loop over the segments
	      for (CSCSegmentCollection::const_iterator segment = range.first; segment!=range.second; segment++)
		 addMuonSegmentMatch(*matchedChamber, &(*segment), parameters);
	   }else{
	      throw cms::Exception("FatalError") << "Failed to cast GeomDet object to either DTChamber or CSCChamber. Who is this guy anyway?\n";
	   }
	}
	info.chambers.push_back(*matchedChamber);
     }
}


void TrackDetectorAssociator::addMuonSegmentMatch(MuonChamberMatch& matchedChamber,
					  const RecSegment* segment,
					  const AssociatorParameters& parameters)
{
   LogTrace("TrackAssociator")
     << "Segment local position: " << segment->localPosition() << "\n"
     << std::hex << segment->geographicalId().rawId() << "\n";
   
   const GeomDet* chamber = muonDetIdAssociator_.getGeomDet(matchedChamber.id);
   TrajectoryStateOnSurface trajectoryStateOnSurface = matchedChamber.tState;
   GlobalPoint segmentGlobalPosition = chamber->toGlobal(segment->localPosition());

   LogTrace("TrackAssociator")
     << "Segment global position: " << segmentGlobalPosition << " \t (R_xy,eta,phi): "
     << segmentGlobalPosition.perp() << "," << segmentGlobalPosition.eta() << "," << segmentGlobalPosition.phi() << "\n";

   LogTrace("TrackAssociator")
     << "\teta hit: " << segmentGlobalPosition.eta() << " \tpropagator: " << trajectoryStateOnSurface.freeState()->position().eta() << "\n"
     << "\tphi hit: " << segmentGlobalPosition.phi() << " \tpropagator: " << trajectoryStateOnSurface.freeState()->position().phi() << std::endl;

   bool isGood = false;
   bool isDTOuterStation = false;
   if( const DTChamber* chamberDT = dynamic_cast<const DTChamber*>(chamber))
     if (chamberDT->id().station()==4)
       isDTOuterStation = true;
   if( isDTOuterStation )
     {
	isGood = fabs(segmentGlobalPosition.phi()-trajectoryStateOnSurface.freeState()->position().phi()) < parameters.dRMuon;
	// Be in chamber
	isGood &= fabs(segmentGlobalPosition.eta()-trajectoryStateOnSurface.freeState()->position().eta()) < .3;
     } else isGood = sqrt( pow(segmentGlobalPosition.eta()-trajectoryStateOnSurface.freeState()->position().eta(),2) + 
			   pow(segmentGlobalPosition.phi()-trajectoryStateOnSurface.freeState()->position().phi(),2)) < parameters.dRMuon;

   if(isGood) {
      MuonSegmentMatch muonSegment;
      muonSegment.segmentGlobalPosition = getPoint(segmentGlobalPosition);
      muonSegment.segmentLocalPosition = getPoint( segment->localPosition() );
      muonSegment.segmentLocalDirection = getVector( segment->localDirection() );
      muonSegment.segmentLocalErrorXX = segment->localPositionError().xx();
      muonSegment.segmentLocalErrorYY = segment->localPositionError().yy();
      muonSegment.segmentLocalErrorXY = segment->localPositionError().xy();
      muonSegment.segmentLocalErrorDxDz = segment->localDirectionError().xx();
      muonSegment.segmentLocalErrorDyDz = segment->localDirectionError().yy();
      
      // DANGEROUS - compiler cannot guaranty parameters ordering
      AlgebraicSymMatrix segmentCovMatrix = segment->parametersError();
      muonSegment.segmentLocalErrorXDxDz = segmentCovMatrix[2][0];
      muonSegment.segmentLocalErrorYDyDz = segmentCovMatrix[3][1];

      matchedChamber.segments.push_back(muonSegment);
   }
}

//********************** NON-CORE CODE ******************************//

std::vector<EcalRecHit> TrackDetectorAssociator::associateEcal( const edm::Event& iEvent,
							const edm::EventSetup& iSetup,
							const FreeTrajectoryState& trackOrigin,
							const double dR )
{
   AssociatorParameters parameters;
   parameters.useHcal = false;
   parameters.useMuon = false;
   parameters.dREcal = dR;
   TrackDetMatchInfo info( associate(iEvent, iSetup, trackOrigin, parameters ));
   if (dR>0) 
     return info.ecalRecHits;
   else
     return info.crossedEcalRecHits;
}

double TrackDetectorAssociator::getEcalEnergy( const edm::Event& iEvent,
				       const edm::EventSetup& iSetup,
				       const FreeTrajectoryState& trackOrigin,
				       const double dR )
{
   AssociatorParameters parameters;
   parameters.useHcal = false;
   parameters.useMuon = false;
   parameters.dREcal = dR;
   TrackDetMatchInfo info = associate(iEvent, iSetup, trackOrigin, parameters );
   if(dR>0) 
     return info.ecalConeEnergy();
   else
     return info.ecalEnergy();
}

std::vector<CaloTower> TrackDetectorAssociator::associateHcal( const edm::Event& iEvent,
						       const edm::EventSetup& iSetup,
						       const FreeTrajectoryState& trackOrigin,
						       const double dR )
{
   AssociatorParameters parameters;
   parameters.useEcal = false;
   parameters.useMuon = false;
   parameters.dRHcal = dR;
   TrackDetMatchInfo info( associate(iEvent, iSetup, trackOrigin, parameters ));
   if (dR>0) 
     return info.towers;
   else
     return info.crossedTowers;
   
}

double TrackDetectorAssociator::getHcalEnergy( const edm::Event& iEvent,
				       const edm::EventSetup& iSetup,
				       const FreeTrajectoryState& trackOrigin,
				       const double dR )
{
   AssociatorParameters parameters;
   parameters.useEcal = false;
   parameters.useMuon = false;
   parameters.dRHcal = dR;
   TrackDetMatchInfo info( associate(iEvent, iSetup, trackOrigin, parameters ));
   if (dR>0) 
     return info.hcalConeEnergy();
   else
     return info.hcalEnergy();
}
