//
// Package:         RecoTracker/SingleTrackPattern
// Class:           CosmicTrajectoryBuilder
// Original Author:  Michele Pioppi-INFN perugia
#include <vector>
#include <iostream>
#include <cmath>

#include "RecoTracker/SingleTrackPattern/interface/CosmicTrajectoryBuilder.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Vector/interface/GlobalPoint.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h" 
#include "TrackingTools/Records/interface/TransientRecHitRecord.h" 

CosmicTrajectoryBuilder::CosmicTrajectoryBuilder(const edm::ParameterSet& conf) : conf_(conf) { 
  //minimum number of hits per tracks

  theMinHits=conf_.getParameter<int>("MinHits");
  //cut on chi2
  chi2cut=conf_.getParameter<double>("Chi2Cut");
  //cut on final chi2
  finalchi2cut=conf_.getParameter<double>("FinalChi2Cut");
  edm::LogInfo("CosmicTrackFinder")<<"Minimum number of hits "<<theMinHits<<" Cut on Chi2= "<<chi2cut;
}


CosmicTrajectoryBuilder::~CosmicTrajectoryBuilder() {
}


void CosmicTrajectoryBuilder::init(const edm::EventSetup& es, bool seedplus){


  //services
  es.get<IdealMagneticFieldRecord>().get(magfield);
   es.get<TrackerDigiGeometryRecord>().get(tracker);

 
 
   if (seedplus) { 	 
     seed_plus=true; 	 
     thePropagator=      new PropagatorWithMaterial(alongMomentum,0.1057,&(*magfield) ); 	 
     thePropagatorOp=    new PropagatorWithMaterial(oppositeToMomentum,0.1057,&(*magfield) );} 	 
   else {
     seed_plus=false;
     thePropagator=      new PropagatorWithMaterial(oppositeToMomentum,0.1057,&(*magfield) ); 	
     thePropagatorOp=    new PropagatorWithMaterial(alongMomentum,0.1057,&(*magfield) );
   }

   theUpdator=       new KFUpdator();
   theEstimator=     new Chi2MeasurementEstimator(chi2cut);


   edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
   std::string builderName = conf_.getParameter<std::string>("TTRHBuilder");   
   es.get<TransientRecHitRecord>().get(builderName,theBuilder);
   

   RHBuilder=   theBuilder.product();




   theFitter=        new KFTrajectoryFitter(*thePropagator,
					    *theUpdator,	
					    *theEstimator) ;


    theSmoother=      new KFTrajectorySmoother(*thePropagatorOp,
					       *theUpdator,	
					       *theEstimator);
 
}

void CosmicTrajectoryBuilder::run(const TrajectorySeedCollection &collseed,
				  const SiStripRecHit2DLocalPosCollection &collstereo,
				  const SiStripRecHit2DLocalPosCollection &collrphi ,
				  const SiStripRecHit2DMatchedLocalPosCollection &collmatched,
				  const SiPixelRecHitCollection &collpixel,
				  const edm::EventSetup& es,
				  edm::Event& e,
				  vector<AlgoProduct> &algooutput)
{

  hits.clear();
  trajFit.clear();

  //order all the hits
  vector<const TrackingRecHit*> allHits= SortHits(collstereo,collrphi,collmatched,collpixel,collseed);
  //vector<const TrackingRecHit*> allHits= SortHits(collstereo,collrphi,collmatched,collseed);  
  std::vector<Trajectory> trajSmooth;
  std::vector<Trajectory>::iterator trajIter;
  
 
  TrajectorySeedCollection::const_iterator iseed;
  for(iseed=collseed.begin();iseed!=collseed.end();iseed++){
    Trajectory startingTraj = createStartingTrajectory(*iseed);
    AddHit(startingTraj,allHits);
 
    for (trajIter=trajFit.begin(); trajIter!=trajFit.end();trajIter++){
      trajSmooth=theSmoother->trajectories((*trajIter));
    }

    for (trajIter= trajSmooth.begin(); trajIter!=trajSmooth.end();trajIter++){
      if((*trajIter).isValid()){
	if ((*trajIter).chiSquared()<500){
	  AlgoProduct tk=makeTrack((*trajIter));
	  algooutput.push_back(tk);
	}else edm::LogError("CosmicTrackFinder")<<"Chi square of the track = "
						<<(*trajIter).chiSquared()
						<<" too high";
      }
    }
  }
};

Trajectory CosmicTrajectoryBuilder::createStartingTrajectory( const TrajectorySeed& seed) const
{
 
  Trajectory result( seed, seed.direction());

  std::vector<TM> seedMeas = seedMeasurements(seed);
  if ( !seedMeas.empty()) {
    for (std::vector<TM>::const_iterator i=seedMeas.begin(); i!=seedMeas.end(); i++){
      result.push(*i);
    }
  }
 
  return result;
}


std::vector<TrajectoryMeasurement> 
CosmicTrajectoryBuilder::seedMeasurements(const TrajectorySeed& seed) const
{
  std::vector<TrajectoryMeasurement> result;

  TrajectorySeed::range hitRange = seed.recHits();
  for (TrajectorySeed::const_iterator ihit = hitRange.first; 
       ihit != hitRange.second; ihit++) {
 
    TransientTrackingRecHit* recHit = RHBuilder->build(&(*ihit));
    
    const GeomDet* hitGeomDet = (&(*tracker))->idToDet( ihit->geographicalId());
    TSOS invalidState( new BasicSingleTrajectoryState( hitGeomDet->surface()));

    if (ihit == hitRange.second - 1) {
      TSOS  updatedState=startingTSOS(seed);
      result.push_back(TM( invalidState, updatedState, recHit));
      
    } 
    else {
      result.push_back(TM( invalidState, recHit));
    }
    
  }

  return result;
};





vector<const TrackingRecHit*> 
CosmicTrajectoryBuilder::SortHits(const SiStripRecHit2DLocalPosCollection &collstereo,
				  const SiStripRecHit2DLocalPosCollection &collrphi ,
				  const SiStripRecHit2DMatchedLocalPosCollection &collmatched,
				  const SiPixelRecHitCollection &collpixel,
				  const TrajectorySeedCollection &collseed){


  //The Hits with global y more than the seed are discarded
  //The Hits correspondign to the seed are discarded
  //At the end all the hits are sorted in y
  vector<const TrackingRecHit*> allHits;

  SiStripRecHit2DLocalPosCollection::const_iterator istrip;
  TrajectorySeedCollection::const_iterator seedbegin=collseed.begin();
  TrajectorySeed::range hRange= (*seedbegin).recHits();
  TrajectorySeed::const_iterator ihit;
  float yref=0.;
  for (ihit = hRange.first; 
       ihit != hRange.second; ihit++) {
    yref=RHBuilder->build(&(*ihit))->globalPosition().y();
    hits.push_back((RHBuilder->build(&(*ihit)))); 
    LogDebug("CosmicTrackFinder")<<"SEED HITS"<<RHBuilder->build(&(*ihit))->globalPosition();
  }

  
  if ((&collpixel)!=0){
    SiPixelRecHitCollection::const_iterator ipix;
    for(ipix=collpixel.begin();ipix!=collpixel.end();ipix++){
      float ych= RHBuilder->build(&(*ipix))->globalPosition().y();
      if ((seed_plus && (ych<yref)) || (!(seed_plus) && (ych>yref)))
	allHits.push_back(&(*ipix));
    }
  } 
  
  

  if ((&collrphi)!=0){
    for(istrip=collrphi.begin();istrip!=collrphi.end();istrip++){
      float ych= RHBuilder->build(&(*istrip))->globalPosition().y();
      if ((seed_plus && (ych<yref)) || (!(seed_plus) && (ych>yref)))
	allHits.push_back(&(*istrip));   
    }
  }




  if ((&collstereo)!=0){
    for(istrip=collstereo.begin();istrip!=collstereo.end();istrip++){
      float ych= RHBuilder->build(&(*istrip))->globalPosition().y();
      if ((seed_plus && (ych<yref)) || (!(seed_plus) && (ych>yref)))
	allHits.push_back(&(*istrip));
    }
  }

  SiStripRecHit2DMatchedLocalPosCollection::const_iterator istripm;
  if ((&collmatched)!=0){
    for(istripm=collmatched.begin();istripm!=collmatched.end();istripm++){
      float ych= RHBuilder->build(&(*istripm))->globalPosition().y();
      if ((seed_plus && (ych<yref)) || (!(seed_plus) && (ych>yref)))
	allHits.push_back(&(*istripm));
    }
  }

  if (seed_plus){
    stable_sort(allHits.begin(),allHits.end(),CompareHitY_plus(*tracker));
  }
  else {
    stable_sort(allHits.begin(),allHits.end(),CompareHitY(*tracker));
  }

  return allHits;
};

TrajectoryStateOnSurface
CosmicTrajectoryBuilder::startingTSOS(const TrajectorySeed& seed)const
{
  PTrajectoryStateOnDet pState( seed.startingState());
  const GeomDet* gdet  = (&(*tracker))->idToDet(DetId(pState.detId()));
  TSOS  State= tsTransform.transientState( pState, &(gdet->surface()), 
					   &(*magfield));
  return State;

}

void CosmicTrajectoryBuilder::AddHit(Trajectory &traj,
				     vector<const TrackingRecHit*>Hits){

  for (unsigned int icosmhit=0;icosmhit<Hits.size();icosmhit++){
    GlobalPoint gphit=RHBuilder->build(Hits[icosmhit])->globalPosition();
    unsigned int iraw= Hits[icosmhit]->geographicalId().rawId();
    LogDebug("CosmicTrackFinder")<<" HIT POSITION "<< gphit;
    TransientTrackingRecHit* tmphit=RHBuilder->build(Hits[icosmhit]);

     TSOS prSt= thePropagator->propagate(traj.lastMeasurement().updatedState(),
 					tracker->idToDet(Hits[icosmhit]->geographicalId())->surface());
     
     if (prSt.isValid()){
       LogDebug("CosmicTrackFinder") <<"STATE PROPAGATED AT DET "<<iraw<<prSt;
       float contr= theEstimator->estimate(prSt, *tmphit).second;
       if (contr<chi2cut)
	 {	 
	   TSOS UpdatedState= theUpdator->update( prSt, *tmphit);
	   if (UpdatedState.isValid()){
	     LogDebug("CosmicTrackFinder") <<"STATE UPDATED WITH HIT AT POSITION "<<gphit<<UpdatedState;
	     traj.push(TM(prSt,UpdatedState,RHBuilder->build(Hits[icosmhit])
			  , contr));
	     hits.push_back(&(*tmphit));
	     LogDebug("CosmicTrackFinder") <<"HIT SELECTED  position" <<gphit
					   <<" traj Chi2= "<<traj.chiSquared();
	     
	   }
	 }else edm::LogError("CosmicTrackFinder")<<" State can not be updated with hit at position "<<gphit;
     }else edm::LogError("CosmicTrackFinder")<<" State can not be propagated at det "<< iraw;
     
     
  }
  
  
  

  if ( qualityFilter( traj)){
    const TrajectorySeed& tmpseed=traj.seed();
    TSOS startingState=startingTSOS(tmpseed);     
    trajFit = theFitter->fit(tmpseed,hits, startingState );
  }
  
  
}


bool 
CosmicTrajectoryBuilder::qualityFilter(Trajectory traj){
  if ( traj.foundHits() >= theMinHits) {
    return true;
  }
  else {
    return false;
  }
}

std::pair<Trajectory, reco::Track>  CosmicTrajectoryBuilder::makeTrack(const Trajectory &traj){
 
  TSOS outertsos = traj.lastMeasurement().updatedState();
  TSOS Fitsos = traj.firstMeasurement().updatedState();

 
  int ndof =traj.foundHits()-5;
  if (ndof<0) ndof=0;

  TSCPBuilderNoMaterial tscpBuilder;
  TrajectoryStateClosestToPoint tscp;

  if (seed_plus) 	 
    tscp= tscpBuilder(*(outertsos.freeState()), 	 
		      outertsos.globalPosition());
  else
    tscp= tscpBuilder(*(Fitsos.freeState()),
		      Fitsos.globalPosition());

  reco::perigee::Parameters param = tscp.perigeeParameters();
  
  reco::perigee::Covariance covar = tscp.perigeeError();
  
 
  reco::Track theTrack(traj.chiSquared(),
		       int(ndof),
		       traj.foundHits(),
		       0,
		       traj.lostHits(),
		       param,
		       covar);
  
  AlgoProduct aProduct(traj,theTrack);
  return aProduct;

}
