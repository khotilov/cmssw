
#include "RecoTracker/CkfPattern/interface/GroupedCkfTrajectoryBuilder.h"
#include "RecoTracker/CkfPattern/interface/TrajectorySegmentBuilder.h"


#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"

#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "RecoTracker/CkfPattern/interface/GroupedTrajCandLess.h"
#include "TrackingTools/TrajectoryFiltering/interface/RegionalTrajectoryFilter.h"
#include "TrackingTools/PatternTools/interface/TempTrajectory.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "TrackingTools/MeasurementDet/interface/LayerMeasurements.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/DetLayers/interface/DetGroup.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/TrajectoryState/interface/BasicSingleTrajectoryState.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"

// only included for RecHit comparison operator:
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"

#include <algorithm> 

using namespace std;

//#define DBG2_GCTB

//#define STANDARD_INTERMEDIARYCLEAN

#ifdef STANDARD_INTERMEDIARYCLEAN
#include "RecoTracker/CkfPattern/interface/IntermediateTrajectoryCleaner.h"
#endif

/* ====== B.M. to be ported layer ===========
#ifdef DBG_GCTB
#include "RecoTracker/CkfPattern/src/ShowCand.h"
#endif
// #define DBG2_GCTB
#ifdef DBG2_GCTB
#include "RecoTracker/CkfPattern/src/SimIdPrinter.h"
#include "Tracker/TkDebugTools/interface/LayerFinderByDet.h"
#include "Tracker/TkLayout/interface/TkLayerName.h"
#endif
=================================== */


GroupedCkfTrajectoryBuilder::
GroupedCkfTrajectoryBuilder(const edm::ParameterSet&              conf,
			    const TrajectoryStateUpdator*         updator,
			    const Propagator*                     propagatorAlong,
			    const Propagator*                     propagatorOpposite,
			    const Chi2MeasurementEstimatorBase*   estimator,
			    const TransientTrackingRecHitBuilder* recHitBuilder,
			    const MeasurementTracker*             measurementTracker,
			    const TrajectoryFilter*               filter,
			    const TrajectoryFilter*               inOutFilter):


  BaseCkfTrajectoryBuilder(conf,
			   updator, propagatorAlong,propagatorOpposite,
			   estimator, recHitBuilder, measurementTracker, filter, inOutFilter)
{
  // fill data members from parameters (eventually data members could be dropped)
  //
  theMaxCand                  = conf.getParameter<int>("maxCand");

  theLostHitPenalty           = conf.getParameter<double>("lostHitPenalty");
  theFoundHitBonus            = conf.getParameter<double>("foundHitBonus");
  theIntermediateCleaning     = conf.getParameter<bool>("intermediateCleaning");
  theAlwaysUseInvalid         = conf.getParameter<bool>("alwaysUseInvalidHits");
  theLockHits                 = conf.getParameter<bool>("lockHits");
  theBestHitOnly              = conf.getParameter<bool>("bestHitOnly");
  theMinNrOf2dHitsForRebuild  = 2;
  theRequireSeedHitsInRebuild = conf.getParameter<bool>("requireSeedHitsInRebuild");
  theMinNrOfHitsForRebuild    = max(0,conf.getParameter<int>("minNrOfHitsForRebuild"));

  /* ======= B.M. to be ported layer ===========
  bool setOK = thePropagator->setMaxDirectionChange(1.6);
  if (!setOK) 
    cout  << "GroupedCkfTrajectoryBuilder WARNING: "
	  << "propagator does not support setMaxDirectionChange" 
	  << endl;
  //   addStopCondition(theMinPtStopCondition);

  theConfigurableCondition = createAlgo<TrajectoryFilter>(componentConfig("StopCondition"));
  ===================================== */

}

/*
  void GroupedCkfTrajectoryBuilder::setEvent(const edm::Event& event) const
  {
  theMeasurementTracker->update(event);
}
*/

GroupedCkfTrajectoryBuilder::TrajectoryContainer 
GroupedCkfTrajectoryBuilder::trajectories (const TrajectorySeed& seed) const 
{
  TrajectoryContainer ret; 
  ret.reserve(10);
  buildTrajectories(seed, ret, 0);
  return ret; 
}

GroupedCkfTrajectoryBuilder::TrajectoryContainer 
GroupedCkfTrajectoryBuilder::trajectories (const TrajectorySeed& seed, 
					   const TrackingRegion& region) const
{
  TrajectoryContainer ret; 
  ret.reserve(10);
  RegionalTrajectoryFilter regionalCondition(region);
  buildTrajectories(seed, ret, &regionalCondition);
  return ret; 
}

void 
GroupedCkfTrajectoryBuilder::trajectories (const TrajectorySeed& seed, GroupedCkfTrajectoryBuilder::TrajectoryContainer &ret) const 
{
  buildTrajectories(seed,ret,0);
}

void
GroupedCkfTrajectoryBuilder::trajectories (const TrajectorySeed& seed, 
                                            GroupedCkfTrajectoryBuilder::TrajectoryContainer &ret,
					    const TrackingRegion& region) const
{
  RegionalTrajectoryFilter regionalCondition(region);
  buildTrajectories(seed,ret,&regionalCondition);
}

void  
GroupedCkfTrajectoryBuilder::rebuildSeedingRegion(const TrajectorySeed& seed,
						  TrajectoryContainer& result) const {    
  TempTrajectory startingTraj = createStartingTrajectory( seed);
  TempTrajectoryContainer work;

  TrajectoryContainer final;

  work.reserve(result.size());
  for (TrajectoryContainer::iterator traj=result.begin();
       traj!=result.end(); ++traj) {
    if(traj->isValid()) work.push_back(TempTrajectory(*traj));
  }

  rebuildSeedingRegion(startingTraj,work);
  final.reserve(work.size());

  for (TempTrajectoryContainer::iterator traj=work.begin();
       traj!=work.end(); ++traj) {
    final.push_back(traj->toTrajectory());
  }
  
  result.swap(final);
  
}

void
GroupedCkfTrajectoryBuilder::buildTrajectories (const TrajectorySeed& seed,
                                                GroupedCkfTrajectoryBuilder::TrajectoryContainer &result,
						const TrajectoryFilter* regionalCondition) const
{
  //
  // Build trajectory outwards from seed
  //

  analyseSeed( seed);

  TempTrajectory startingTraj = createStartingTrajectory( seed);

  work_.clear();
  const bool inOut = true;
  groupedLimitedCandidates( startingTraj, regionalCondition, theForwardPropagator, inOut, work_);
  if ( work_.empty() )  return ;



  /*  rebuilding is de-coupled from standard building
  //
  // try to additional hits in the seeding region
  //
  if ( theMinNrOfHitsForRebuild>0 ) {
    // reverse direction
    //thePropagator->setPropagationDirection(oppositeDirection(seed.direction()));
    // rebuild part of the trajectory
    rebuildSeedingRegion(startingTraj,work);
  }

  */

  result.reserve(work_.size());
  for (TempTrajectoryContainer::const_iterator it = work_.begin(), ed = work_.end(); it != ed; ++it) {
      result.push_back( it->toTrajectory() );
  }

  work_.clear(); 
  if (work_.capacity() > work_MaxSize_) {  TempTrajectoryContainer().swap(work_); work_.reserve(work_MaxSize_/2); }

  analyseResult(result);

  LogDebug("CkfPattern")<< "GroupedCkfTrajectoryBuilder: returning result of size " << result.size();

}


void 
GroupedCkfTrajectoryBuilder::groupedLimitedCandidates (TempTrajectory& startingTraj, 
						       const TrajectoryFilter* regionalCondition,
						       const Propagator* propagator, 
                                                       bool inOut,
						       TempTrajectoryContainer& result) const
{
  unsigned int nIter=1;
  TempTrajectoryContainer candidates;
  TempTrajectoryContainer newCand;
  candidates.push_back( startingTraj);

  while ( !candidates.empty()) {

    newCand.clear();
    for (TempTrajectoryContainer::iterator traj=candidates.begin();
	 traj!=candidates.end(); traj++) {
      if ( !advanceOneLayer(*traj, regionalCondition, propagator, inOut, newCand, result) ) {
	LogDebug("CkfPattern")<< "GCTB: terminating after advanceOneLayer==false";
 	continue;
      }

      LogDebug("CkfPattern")<<"newCand(1): after advanced one layer:\n"<<PrintoutHelper::dumpCandidates(newCand);

      if ((int)newCand.size() > theMaxCand) {
	//ShowCand()(newCand);

 	sort( newCand.begin(), newCand.end(), GroupedTrajCandLess(theLostHitPenalty,theFoundHitBonus));
 	newCand.erase( newCand.begin()+theMaxCand, newCand.end());
      }
      LogDebug("CkfPattern")<<"newCand(2): after removing extra candidates.\n"<<PrintoutHelper::dumpCandidates(newCand);
    }

    LogDebug("CkfPattern") << "newCand.size() at end = " << newCand.size();
/*
    if (theIntermediateCleaning) {
      candidates.clear();
      candidates = groupedIntermediaryClean(newCand);
    } else {
      candidates.swap(newCand);
    }
*/
    if (theIntermediateCleaning) {
#ifdef STANDARD_INTERMEDIARYCLEAN
	IntermediateTrajectoryCleaner::clean(newCand);	
#else 
	groupedIntermediaryClean(newCand);
#endif	

    }	
    candidates.swap(newCand);

    LogDebug("CkfPattern") <<"candidates(3): "<<result.size()<<" candidates after "<<nIter++<<" groupedCKF iteration: \n"
      			   <<PrintoutHelper::dumpCandidates(result)
			   <<"\n "<<candidates.size()<<" running candidates are: \n"
			   <<PrintoutHelper::dumpCandidates(candidates);
  }
}

std::string whatIsTheNextStep(TempTrajectory& traj , std::pair<TrajectoryStateOnSurface,std::vector<const DetLayer*> >& stateAndLayers){
  std::stringstream buffer;
  vector<const DetLayer*> & nl = stateAndLayers.second;
#include "TrackingTools/DetLayers/interface/BarrelDetLayer.h"
#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
  //B.M. TkLayerName layerName;
  //B.M. buffer << "Started from " << layerName(traj.lastLayer()) 
  const BarrelDetLayer* sbdl = dynamic_cast<const BarrelDetLayer*>(traj.lastLayer());
  const ForwardDetLayer* sfdl = dynamic_cast<const ForwardDetLayer*>(traj.lastLayer());
  if (sbdl) {
    buffer << "Started from " << traj.lastLayer() << " r=" << sbdl->specificSurface().radius() 
	   << " phi=" << sbdl->specificSurface().phi() << endl;
  } else if (sfdl) {
    buffer << "Started from " << traj.lastLayer() << " z " << sfdl->specificSurface().position().z()
	   << " phi " << sfdl->specificSurface().phi() << endl;
  }
  buffer << "Trying to go to";
  for ( vector<const DetLayer*>::iterator il=nl.begin();
	il!=nl.end(); il++){ 
    //B.M. buffer << " " << layerName(*il)  << " " << *il << endl;
    const BarrelDetLayer* bdl = dynamic_cast<const BarrelDetLayer*>(*il);
    const ForwardDetLayer* fdl = dynamic_cast<const ForwardDetLayer*>(*il);
    
    if (bdl) buffer << " r " << bdl->specificSurface().radius() << endl;
    if (fdl) buffer << " z " << fdl->specificSurface().position().z() << endl;
    //buffer << " " << *il << endl;   
  }
  return buffer.str();
}

std::string whatIsTheStateToUse(TrajectoryStateOnSurface &initial, TrajectoryStateOnSurface & stateToUse, const DetLayer * l){
  std::stringstream buffer;
  buffer << "GCTB: starting from " 
         << " r / phi / z = " << stateToUse.globalPosition().perp()
	 << " / " << stateToUse.globalPosition().phi()
	 << " / " << stateToUse.globalPosition().z() 
         << " , pt / phi / pz /charge = " 
	 << stateToUse.globalMomentum().perp() << " / "  
	 << stateToUse.globalMomentum().phi() << " / " 
	 << stateToUse.globalMomentum().z() << " / " 
	 << stateToUse.charge()
	 << " for layer at "<< l << endl;
  buffer << "     errors:";
  for ( int i=0; i<5; i++ )  buffer << " " << sqrt(stateToUse.curvilinearError().matrix()(i,i));
  buffer << endl;
  
  //buffer << "GCTB: starting from r / phi / z = " << initial.globalPosition().perp()
  //<< " / " << initial.globalPosition().phi()
  //<< " / " << initial.globalPosition().z() << " , pt / pz = " 
  //<< initial.globalMomentum().perp() << " / " 
  //<< initial.globalMomentum().z() << " for layer at "
  //<< l << endl;
  //buffer << "     errors:";
  //for ( int i=0; i<5; i++ )  buffer << " " << sqrt(initial.curvilinearError().matrix()(i,i));
  //buffer << endl;
  return buffer.str();
}


bool 
GroupedCkfTrajectoryBuilder::advanceOneLayer (TempTrajectory& traj, 
					      const TrajectoryFilter* regionalCondition, 
					      const Propagator* propagator,
                                              bool inOut,
					      TempTrajectoryContainer& newCand, 
					      TempTrajectoryContainer& result) const
{
  std::pair<TSOS,std::vector<const DetLayer*> > stateAndLayers = findStateAndLayers(traj);
  vector<const DetLayer*>::iterator layerBegin = stateAndLayers.second.begin();
  vector<const DetLayer*>::iterator layerEnd   = stateAndLayers.second.end();

  //   if (nl.empty()) {
  //     addToResult(traj,result,inOut);
  //     return false;
  //   }
  
  LogDebug("CkfPattern")<<whatIsTheNextStep(traj, stateAndLayers);
  
  bool foundSegments(false);
  bool foundNewCandidates(false);
  for ( vector<const DetLayer*>::iterator il=layerBegin; 
	il!=layerEnd; il++) {

    TSOS stateToUse = stateAndLayers.first;
    if ((*il)==traj.lastLayer())
      {
	LogDebug("CkfPattern")<<" self propagating in advanceOneLayer.\n from: \n"<<stateToUse;
	//self navigation case
	// go to a middle point first
	TransverseImpactPointExtrapolator middle;
	GlobalPoint center(0,0,0);
	stateToUse = middle.extrapolate(stateToUse, center, *theForwardPropagator);
	
	if (!stateToUse.isValid()) continue;
	LogDebug("CkfPattern")<<"to: "<<stateToUse;
      }
    
    TrajectorySegmentBuilder layerBuilder(theMeasurementTracker,
					  theLayerMeasurements,
					  **il,*propagator,
					  *theUpdator,*theEstimator,
					  theLockHits,theBestHitOnly);

    LogDebug("CkfPattern")<<whatIsTheStateToUse(stateAndLayers.first,stateToUse,*il);
    
    TempTrajectoryContainer segments=
      layerBuilder.segments(stateToUse);

    LogDebug("CkfPattern")<< "GCTB: number of segments = " << segments.size();

    if ( !segments.empty() )  foundSegments = true;

    for ( TempTrajectoryContainer::const_iterator is=segments.begin();
	  is!=segments.end(); is++ ) {
      //
      // assume "invalid hit only" segment is last in list
      //
      const TempTrajectory::DataContainer & measurements = is->measurements();
      if ( !theAlwaysUseInvalid && is!=segments.begin() && measurements.size()==1 && 
	   (measurements.front().recHit()->getType() == TrackingRecHit::missing) )  break;
      //
      // create new candidate
      //
      TempTrajectory newTraj(traj);
      
      newTraj.push(*is);
      //GIO// for ( vector<TM>::const_iterator im=measurements.begin();
      //GIO//        im!=measurements.end(); im++ )  newTraj.push(*im);
      //if ( toBeContinued(newTraj,regionalCondition) ) { TOBE FIXED
      if ( toBeContinued(newTraj, inOut) ) {
	// Have added one more hit to track candidate
	
	LogDebug("CkfPattern")<<"GCTB: adding updated trajectory to candidates: inOut="<<inOut<<" hits="<<newTraj.foundHits();

	newCand.push_back(newTraj);
	foundNewCandidates = true;
      }
      else {
	// Have finished building this track. Check if it passes cuts.

	LogDebug("CkfPattern")<< "GCTB: adding completed trajectory to results if passes cuts: inOut="<<inOut<<" hits="<<newTraj.foundHits();

	addToResult(newTraj, result, inOut);
      }
    }
  }

  if ( !foundSegments ){
    LogDebug("CkfPattern")<< "GCTB: adding input trajectory to result";
    addToResult(traj, result, inOut);
  }
  return foundNewCandidates;
}

//TempTrajectoryContainer
void
GroupedCkfTrajectoryBuilder::groupedIntermediaryClean (TempTrajectoryContainer& theTrajectories) const 
{
  //if (theTrajectories.empty()) return TrajectoryContainer();
  //TrajectoryContainer result;
  if (theTrajectories.empty()) return;  
  //RecHitEqualByChannels recHitEqualByChannels(false, false);
  int firstLayerSize, secondLayerSize;
  vector<const DetLayer*> firstLayers, secondLayers;

  for (TempTrajectoryContainer::iterator firstTraj=theTrajectories.begin();
       firstTraj!=(theTrajectories.end()-1); firstTraj++) {

    if ( (!firstTraj->isValid()) ||
         (!firstTraj->lastMeasurement().recHit()->isValid()) ) continue;
    const TempTrajectory::DataContainer & firstMeasurements = firstTraj->measurements();
    layers(firstMeasurements, firstLayers);
    firstLayerSize = firstLayers.size();
    if ( firstLayerSize<4 )  continue;

    for (TempTrajectoryContainer::iterator secondTraj=(firstTraj+1);
       secondTraj!=theTrajectories.end(); secondTraj++) {

      if ( (!secondTraj->isValid()) ||
           (!secondTraj->lastMeasurement().recHit()->isValid()) ) continue;
      const TempTrajectory::DataContainer & secondMeasurements = secondTraj->measurements();
      layers(secondMeasurements, secondLayers);
      secondLayerSize = secondLayers.size();
      //
      // only candidates using the same last 3 layers are compared
      //
      if ( firstLayerSize!=secondLayerSize )  continue;
      if ( firstLayers[0]!=secondLayers[0] ||
	   firstLayers[1]!=secondLayers[1] ||
	   firstLayers[2]!=secondLayers[2] )  continue;

      TempTrajectory::DataContainer::const_iterator im1 = firstMeasurements.rbegin();
      TempTrajectory::DataContainer::const_iterator im2 = secondMeasurements.rbegin();
      //
      // check for identical hits in the last layer
      //
      bool unequal(false);
      const DetLayer* layerPtr = firstLayers[0];
      while ( im1!=firstMeasurements.rend()&&im2!=secondMeasurements.rend() ) {
	if ( im1->layer()!=layerPtr || im2->layer()!=layerPtr )  break;
	if ( !(im1->recHit()->isValid()) || !(im2->recHit()->isValid()) ||
	     //!recHitEqualByChannels(im1->recHit(),im2->recHit()) ) {
	     !im1->recHit()->hit()->sharesInput(im2->recHit()->hit(), TrackingRecHit::some) ) {
	  unequal = true;
	  break;
	}
	--im1;
	--im2;
      }
      if ( im1==firstMeasurements.rend() || im2==secondMeasurements.rend() ||
	   im1->layer()==layerPtr || im2->layer()==layerPtr || unequal )  continue;
      //
      // check for invalid hits in the layer -2
      // compare only candidates with invalid / valid combination
      //
      layerPtr = firstLayers[1];
      bool firstValid(true);
      while ( im1!=firstMeasurements.rend()&&im1->layer()==layerPtr ) {
	if ( !im1->recHit()->isValid() )  firstValid = false;
	--im1;
      }
      bool secondValid(true);
      while ( im2!=secondMeasurements.rend()&&im2->layer()==layerPtr ) {
	if ( !im2->recHit()->isValid() )  secondValid = false;
	--im2;
      }
      if ( !tkxor(firstValid,secondValid) )  continue;
      //
      // ask for identical hits in layer -3
      //
      unequal = false;
      layerPtr = firstLayers[2];
      while ( im1!=firstMeasurements.rend()&&im2!=secondMeasurements.rend() ) {
	if ( im1->layer()!=layerPtr || im2->layer()!=layerPtr )  break;
	if ( !(im1->recHit()->isValid()) || !(im2->recHit()->isValid()) ||
	     //!recHitEqualByChannels(im1->recHit(),im2->recHit()) ) {
	     !im1->recHit()->hit()->sharesInput(im2->recHit()->hit(), TrackingRecHit::some) ) {
	  unequal = true;
	  break;
	}
	--im1;
	--im2;
      }
      if ( im1==firstMeasurements.rend() || im2==secondMeasurements.rend() ||
	   im1->layer()==layerPtr || im2->layer()==layerPtr || unequal )  continue;

      if ( !firstValid ) {
	firstTraj->invalidate();
	break;
      }
      else {
	secondTraj->invalidate();
	break;
      }
    }
  }
/*
  for (TempTrajectoryContainer::const_iterator it = theTrajectories.begin();
       it != theTrajectories.end(); it++) {
    if(it->isValid()) result.push_back( *it);
  }

  return result;
*/
  theTrajectories.erase(std::remove_if( theTrajectories.begin(),theTrajectories.end(),
                                        std::not1(std::mem_fun_ref(&TempTrajectory::isValid))),
 //                                     boost::bind(&TempTrajectory::isValid,_1)), 
                        theTrajectories.end());
}

void
GroupedCkfTrajectoryBuilder::layers (const TempTrajectory::DataContainer& measurements,
                                     vector<const DetLayer*> &result) const 
{
  result.clear();

  if ( measurements.empty() )  return ;

  result.push_back(measurements.back().layer());
  TempTrajectory::DataContainer::const_iterator ifirst = measurements.rbegin();
  --ifirst;	 
  for ( TempTrajectory::DataContainer::const_iterator im=ifirst;
	im!=measurements.rend(); --im ) {
    if ( im->layer()!=result.back() )  result.push_back(im->layer());
  }

  for (vector<const DetLayer*>::const_iterator iter = result.begin(); iter != result.end(); iter++){
    if (!*iter) edm::LogWarning("CkfPattern")<< "Warning: null det layer!! ";
  }
}

void
GroupedCkfTrajectoryBuilder::rebuildSeedingRegion(TempTrajectory& startingTraj,
						  TempTrajectoryContainer& result) const
{
  //
  // Rebuilding of trajectories. Candidates are taken from result,
  // which will be replaced with the solutions after rebuild
  // (assume vector::swap is more efficient than building new container)
  //
  LogDebug("CkfPattern")<< "Starting to rebuild " << result.size() << " tracks";
  //
  // Fitter (need to create it here since the propagation direction
  // might change between different starting trajectories)
  //
  KFTrajectoryFitter fitter(&(*theBackwardPropagator),&updator(),&estimator());
  //
  TempTrajectoryContainer reFitted;
  TrajectorySeed::range rseedHits = startingTraj.seed().recHits();
  std::vector<const TrackingRecHit*> seedHits;
  //seedHits.insert(seedHits.end(), rseedHits.first, rseedHits.second);
  //for (TrajectorySeed::recHitContainer::const_iterator iter = rseedHits.first; iter != rseedHits.second; iter++){
  //	seedHits.push_back(&*iter);
  //}

  //unsigned int nSeed(seedHits.size());
  unsigned int nSeed(rseedHits.second-rseedHits.first);
  //seedHits.reserve(nSeed);
  TempTrajectoryContainer rebuiltTrajectories;
  for ( TempTrajectoryContainer::iterator it=result.begin();
	it!=result.end(); it++ ) {
    //
    // skip candidates which are not exceeding the seed size 
    // (e.g. because no Tracker layers outside seeding region) 
    //

    if ( it->measurements().size()<=startingTraj.measurements().size() ) {
      rebuiltTrajectories.push_back(*it);
      LogDebug("CkfPattern")<< "RebuildSeedingRegion skipped as in-out trajectory does not exceed seed size.";
      continue;
    }
    //
    // Refit - keep existing trajectory in case fit is not possible
    // or fails
    //
    backwardFit(*it,nSeed,fitter,reFitted,seedHits);
    if ( reFitted.size()!=1 ) {
      rebuiltTrajectories.push_back(*it);
      LogDebug("CkfPattern")<< "RebuildSeedingRegion skipped as backward fit failed";
      //			    << "after reFitted.size() " << reFitted.size();
      continue;
    }
    //LogDebug("CkfPattern")<<"after reFitted.size() " << reFitted.size();
    //
    // Rebuild seeding part. In case it fails: keep initial trajectory
    // (better to drop it??)
    //
    int nRebuilt =
      rebuildSeedingRegion (seedHits,reFitted.front(),rebuiltTrajectories);

    if ( nRebuilt==0 ) it->invalidate();  // won't use original in-out track

    if ( nRebuilt<=0 ) rebuiltTrajectories.push_back(*it);

  }
  //
  // Replace input trajectories with new ones
  //
  result.swap(rebuiltTrajectories);
  result.erase(std::remove_if( result.begin(),result.end(),
			       std::not1(std::mem_fun_ref(&TempTrajectory::isValid))),
	       result.end());
}

int
GroupedCkfTrajectoryBuilder::rebuildSeedingRegion(const std::vector<const TrackingRecHit*>& seedHits, 
						  TempTrajectory& candidate,
						  TempTrajectoryContainer& result) const 
{
  //
  // Starting from track found by in-out tracking phase, extrapolate it inwards through
  // the seeding region if possible in towards smaller Tracker radii, searching for additional
  // hits.
  // The resulting trajectories are returned in result,
  // the count is the return value.
  //
  TempTrajectoryContainer rebuiltTrajectories;
#ifdef DBG2_GCTB
/*  const LayerFinderByDet layerFinder;
  if ( !seedHits.empty() && seedHits.front().isValid() ) {
    DetLayer* seedLayer = layerFinder(seedHits.front().det());
    cout << "Seed hit at " << seedHits.front().globalPosition()
	 << " " << seedLayer << endl;
    cout << "Started from " 
	 << candidate.lastMeasurement().updatedState().globalPosition().perp() << " "
	 << candidate.lastMeasurement().updatedState().globalPosition().z() << endl;
    pair<bool,TrajectoryStateOnSurface> layerComp(false,TrajectoryStateOnSurface());
    if ( seedLayer ) layerComp =
      seedLayer->compatible(candidate.lastMeasurement().updatedState(),
			      propagator(),estimator());
    pair<bool,TrajectoryStateOnSurface> detComp =
      seedHits.front().det().compatible(candidate.lastMeasurement().updatedState(),
					propagator(),estimator());
    cout << "  layer compatibility = " << layerComp.first;
    cout << "  det compatibility = " << detComp.first;
    if ( detComp.first ) {
      cout << "  estimate = " 
	   << estimator().estimate(detComp.second,seedHits.front()).second ;
    }
    cout << endl;
  }*/
  cout << "Before backward building: #measurements = " 
       << candidate.measurements().size() ; //<< endl;;
#endif
  //
  // Use standard building with standard cuts. Maybe better to use different
  // cuts from "forward" building (e.g. no check on nr. of invalid hits)?
  //
  const bool inOut = false;
  groupedLimitedCandidates(candidate, (const TrajectoryFilter*)0, theBackwardPropagator, inOut, rebuiltTrajectories);

  LogDebug("CkfPattern")<<" After backward building: "<<PrintoutHelper::dumpCandidates(rebuiltTrajectories);

  //
  // Check & count resulting candidates
  //
  int nrOfTrajectories(0);
  bool orig_ok = false;
  //const RecHitEqualByChannels recHitEqual(false,false);
  //vector<TM> oldMeasurements(candidate.measurements());
  for ( TempTrajectoryContainer::iterator it=rebuiltTrajectories.begin();
	it!=rebuiltTrajectories.end(); it++ ) {

    TempTrajectory::DataContainer newMeasurements(it->measurements());
    //
    // Verify presence of seeding hits?
    //
    if ( theRequireSeedHitsInRebuild ) {
      orig_ok = true;
      // no hits found (and possibly some invalid hits discarded): drop track
      if ( newMeasurements.size()<=candidate.measurements().size() ){  
	LogDebug("CkfPattern") << "newMeasurements.size()<=candidate.measurements().size()";
	continue;
      }	
      // verify presence of hits
      //GIO//if ( !verifyHits(newMeasurements.begin()+candidate.measurements().size(),
      //GIO//		       newMeasurements.end(),seedHits) ){
      if ( !verifyHits(newMeasurements.rbegin(), 
                       newMeasurements.size() - candidate.measurements().size(),
		       seedHits) ){
	LogDebug("CkfPattern")<< "seed hits not found in rebuild";
	continue; 
      }
    }
    //
    // construct final trajectory in the right order
    //
    TempTrajectory reversedTrajectory(it->seed(),it->seed().direction());
    for (TempTrajectory::DataContainer::const_iterator im=newMeasurements.rbegin(), ed = newMeasurements.rend();
	  im != ed; --im ) {
      reversedTrajectory.push(*im);
    }
    // save & count result
    result.push_back(reversedTrajectory);
    nrOfTrajectories++;

    LogDebug("CkgPattern")<<"New traj direction = " << reversedTrajectory.direction()<<"\n"
			  <<PrintoutHelper::dumpMeasurements(reversedTrajectory.measurements());
  }
  // If nrOfTrajectories = 0 and orig_ok = true, this means that a track was actually found on the
  // out-in step (meeting those requirements) but did not have the seed hits in it.
  // In this case when we return we will go ahead and use the original in-out track.

  // If nrOfTrajectories = 0 and orig_ok = false, this means that the out-in step failed to
  // find any track.  Two cases are a technical failure in fitting the original seed hits or
  // because the track did not meet the out-in criteria (which may be stronger than the out-in
  // criteria).  In this case we will NOT allow the original in-out track to be used.

  if ( (nrOfTrajectories == 0) && orig_ok ) {
    nrOfTrajectories = -1;
  }
  return nrOfTrajectories;
}

void
GroupedCkfTrajectoryBuilder::backwardFit (TempTrajectory& candidate, unsigned int nSeed,
						    const TrajectoryFitter& fitter,
						    TempTrajectoryContainer& fittedTracks,
						    std::vector<const TrackingRecHit*>& remainingHits) const
{
  //
  // clear array of non-fitted hits
  //
  remainingHits.clear();
  fittedTracks.clear();
  //
  // skip candidates which are not exceeding the seed size
  // (e.g. Because no Tracker layers exist outside seeding region)
  //
  if ( candidate.measurements().size()<=nSeed ) {
    fittedTracks.clear();
    return;
  }

  LogDebug("CkfPattern")<<"nSeed " << nSeed << endl
			<< "Old traj direction = " << candidate.direction() << endl
			<<PrintoutHelper::dumpMeasurements(candidate.measurements());

  //
  // backward fit trajectory.
  // (Will try to fit only hits outside the seeding region. However,
  // if there are not enough of these, it will also use the seeding hits).
  //
  TempTrajectory::DataContainer oldMeasurements(candidate.measurements());
//   int nOld(oldMeasurements.size());
//   const unsigned int nHitAllMin(5);
//   const unsigned int nHit2dMin(2);
  unsigned int nHit(0);    // number of valid hits after seeding region
  //unsigned int nHit2d(0);  // number of valid hits after seeding region with 2D info
  // use all hits except the first n (from seed), but require minimum
  // specified in configuration.
  //  Swapped over next two lines.
  unsigned int nHitMin = max(candidate.foundHits()-nSeed,theMinNrOfHitsForRebuild);
  //  unsigned int nHitMin = oldMeasurements.size()-nSeed;
  // we want to rebuild only if the number of VALID measurements excluding the seed measurements is higher than the cut
  if (nHitMin<theMinNrOfHitsForRebuild){
	fittedTracks.clear();
    	return;
  }

  LogDebug("CkfPattern")/* << "nHitMin " << nHitMin*/ <<"Sizes: " << oldMeasurements.size() << " / ";
  //
  // create input trajectory for backward fit
  //
  Trajectory fwdTraj(candidate.seed(),oppositeDirection(candidate.direction()));
  //const TrajectorySeed seed = TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), oppositeDirection(candidate.direction()));
  //Trajectory fwdTraj(seed, oppositeDirection(candidate.direction()));
  std::vector<const DetLayer*> bwdDetLayer; 
  for ( TempTrajectory::DataContainer::const_iterator im=oldMeasurements.rbegin();
	im!=oldMeasurements.rend(); --im) {
    const TrackingRecHit* hit = im->recHit()->hit();
    //
    // add hits until required number is reached
    //
    if ( nHit<nHitMin ){//|| nHit2d<theMinNrOf2dHitsForRebuild ) {
      fwdTraj.push(*im);
      bwdDetLayer.push_back(im->layer());
      //
      // count valid / 2D hits
      //
      if ( hit->isValid() ) {
	nHit++;
	//if ( hit.isMatched() ||
	//     hit.det().detUnits().front()->type().module()==pixel )
        //nHit2d++;
      }
    }
    //if (nHit==nHitMin) lastBwdDetLayer=im->layer();	
    //
    // keep remaining (valid) hits for verification
    //
    else if ( hit->isValid() ) {
      //std::cout << "Adding a remaining hit" << std::endl;
      remainingHits.push_back(hit);
    }
  }
  //
  // Fit only if required number of valid hits can be used
  //
  if ( nHit<nHitMin ){  //|| nHit2d<theMinNrOf2dHitsForRebuild ) {
    fittedTracks.clear();
    return;
  }
  //
  // Do the backward fit (important: start from scaled, not random cov. matrix!)
  //
  TrajectoryStateOnSurface firstTsos(fwdTraj.firstMeasurement().updatedState());
  //cout << "firstTsos "<< firstTsos << endl;
  firstTsos.rescaleError(10.);
  //TrajectoryContainer bwdFitted(fitter.fit(fwdTraj.seed(),fwdTraj.recHits(),firstTsos));
  TrajectoryContainer bwdFitted(fitter.fit(
  		TrajectorySeed(PTrajectoryStateOnDet(), TrajectorySeed::recHitContainer(), oppositeDirection(candidate.direction())),
  		fwdTraj.recHits(),firstTsos));
  if (bwdFitted.size()){
    LogDebug("CkfPattern")<<"Obtained " << bwdFitted.size() << " bwdFitted trajectories with measurement size " << bwdFitted.front().measurements().size();
	TempTrajectory fitted(fwdTraj.seed(), fwdTraj.direction());
        vector<TM> tmsbf = bwdFitted.front().measurements();
	int iDetLayer=0;
	//this is ugly but the TM in the fitted track do not contain the DetLayer.
	//So we have to cache the detLayer pointers and replug them in.
	//For the backward building it would be enaugh to cache the last DetLayer, 
	//but for the intermediary cleaning we need all
 	for ( vector<TM>::const_iterator im=tmsbf.begin();im!=tmsbf.end(); im++ ) {
		fitted.push(TM( (*im).forwardPredictedState(),
				(*im).backwardPredictedState(),
				(*im).updatedState(),
				(*im).recHit(),
				(*im).estimate(),
				bwdDetLayer[iDetLayer]));

		LogDebug("CkfPattern")<<PrintoutHelper::dumpMeasurement(*im);
		iDetLayer++;
	}
/*
	TM lastMeas = bwdFitted.front().lastMeasurement();
	fitted.pop();
	fitted.push(TM(lastMeas.forwardPredictedState(), 
			       lastMeas.backwardPredictedState(), 
			       lastMeas.updatedState(),
			       lastMeas.recHit(),
			       lastMeas.estimate(),
                               lastBwdDetLayer));*/
	fittedTracks.push_back(fitted);
  }
  //
  // save result
  //
  //fittedTracks.swap(bwdFitted);
  //cout << "Obtained " << fittedTracks.size() << " fittedTracks trajectories with measurement size " << fittedTracks.front().measurements().size() << endl;
}

bool
GroupedCkfTrajectoryBuilder::verifyHits (TempTrajectory::DataContainer::const_iterator rbegin,
                                         size_t maxDepth,
					 const std::vector<const TrackingRecHit*>& hits) const
{
  //
  // verify presence of the seeding hits
  //
  LogDebug("CkfPattern")<<"Checking for " << hits.size() << " hits in "
			<< maxDepth << " measurements" << endl;

  TempTrajectory::DataContainer::const_iterator rend = rbegin; 
  while (maxDepth > 0) { --maxDepth; --rend; }
  for ( vector<const TrackingRecHit*>::const_iterator ir=hits.begin();
	ir!=hits.end(); ir++ ) {
    // assume that all seeding hits are valid!
    bool foundHit(false);
    for ( TempTrajectory::DataContainer::const_iterator im=rbegin; im!=rend; --im ) {
      if ( im->recHit()->isValid() && (*ir)->sharesInput(im->recHit()->hit(), TrackingRecHit::some) ) {
	foundHit = true;
	break;
      }
    }
    if ( !foundHit )  return false;
  }
  return true;
}




