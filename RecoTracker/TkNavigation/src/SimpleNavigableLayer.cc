#include "RecoTracker/TkNavigation/interface/SimpleNavigableLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include <set>
#include <TrackingTools/TrackAssociator/interface/DetIdInfo.h>

using namespace std;

TrajectoryStateOnSurface SimpleNavigableLayer::crossingState(const FreeTrajectoryState& fts,
							     PropagationDirection dir) const{
  TSOS propState;
  //self propagating. step one: go close to the center
  GlobalPoint initialPoint = fts.position();
  TransverseImpactPointExtrapolator middle;
  GlobalPoint center(0,0,0);
  propState = middle.extrapolate(fts, center, propagator(dir));
  if ( !propState.isValid()) return TrajectoryStateOnSurface();
  
  FreeTrajectoryState & dest = *propState.freeState();
  GlobalPoint middlePoint = dest.position();
  const double toCloseToEachOther=1e-4;
  if ( (middlePoint-initialPoint).mag() < toCloseToEachOther){
    LogDebug("SimpleNavigableLayer")<<"initial state and PCA are identical. Things are bound to fail. Do not add the link.";
    return TrajectoryStateOnSurface();
  }
  
  std::string dirS;
  if (dir==alongMomentum) dirS = "alongMomentum";
  else if (dir==oppositeToMomentum) dirS = "oppositeToMomentum";
  else dirS = "anyDirection";
  
  LogDebug("SimpleNavigableLayer")<<"self propagating("<< dir <<") from:\n"
				  <<fts<<"\n"
				  <<dest<<"\n"
				  <<" and the direction is: "<<dir<<" = "<<dirS;
  
  //second propagation to go on the other side of the barrel
  //propState = propagator(dir).propagate( dest, detLayer()->specificSurface());
  propState = propagator(dir).propagate( dest, detLayer()->surface());
  if ( !propState.isValid()) return TrajectoryStateOnSurface();
  
  FreeTrajectoryState & dest2 = *propState.freeState();
  GlobalPoint finalPoint = dest2.position();
  LogDebug("SimpleNavigableLayer")<<"second propagation("<< dir <<") to: \n"
				  <<dest2;
  double finalDot = (middlePoint - initialPoint).basicVector().dot((finalPoint-middlePoint).basicVector());
  if (finalDot<0){ // check that before and after are in different side.
    LogDebug("SimpleNavigableLayer")<<"switch side back: ABORT.";
    return TrajectoryStateOnSurface();
  }
  return propState;
}

bool SimpleNavigableLayer::wellInside( const FreeTrajectoryState& fts,
				       PropagationDirection dir,
				       const BarrelDetLayer* bl,
				       DLC& result) const
{
  TSOS propState ;

  if (bl==detLayer()){
    propState = crossingState(fts,dir);
    if (!propState.isValid()) return false;
  }else{
    propState= propagator(dir).propagate( fts, bl->specificSurface());
  }

  if ( !propState.isValid()) return false;
 
  //if requested check that the layer is crossed on the right side
  if (theCheckCrossingSide){
    bool backTobackTransverse = (fts.position().x()*propState.globalPosition().x() + fts.position().y()*propState.globalPosition().y())<0;
    bool backToback = propState.globalPosition().basicVector().dot(fts.position().basicVector())<0;

    if (backTobackTransverse || backToback ){
      LogTrace("TkNavigation") << "Crossing over prevented!\nStaring from (x,y,z,r) (" 
			       << fts.position().x()<<","<< fts.position().y()<<","<< fts.position().z()<<","<<fts.position().perp()
			       << ") going to TSOS (x,y,z,r)" 
			       << propState.globalPosition().x()<<","<< propState.globalPosition().y()<<","<< propState.globalPosition().z()<<","<<propState.globalPosition().perp()<<")";
      return false;
    
    /*
    //we have to check the crossing side only if we are going to something smaller
    if (fts.position().perp()>bl->specificSurface().radius() || 
    fabs(fts.position().z())>bl->surface().bounds().length()/2. ){
    if (propState.globalPosition().basicVector().dot(fts.position().basicVector())<0){
    LogTrace("TkNavigation") << "Crossing over prevented!\nStaring from (x,y,z,r) (" 
    << fts.position().x()<<","<< fts.position().y()<<","<< fts.position().z()<<","<<fts.position().perp()
    << ") going to TSOS (x,y,z,r)" 
    << propState.globalPosition().x()<<","<< propState.globalPosition().y()<<","<< propState.globalPosition().z()<<","<<propState.globalPosition().perp()<<")";; 
    return false;
    }
    } 	
    */
    }}

  const Bounds& bounds( bl->specificSurface().bounds());
  float length = bounds.length() / 2.f;

  // take into account the thickness of the layer
  float deltaZ = bounds.thickness()/2. / 
    fabs( tan( propState.globalDirection().theta()));

  // take into account the error on the predicted state
  const float nSigma = theEpsilon;  // temporary reuse of epsilon
  if (propState.hasError()) {
    deltaZ += nSigma * sqrt( fts.cartesianError().position().czz());
  }

  // cout << "SimpleNavigableLayer BarrelDetLayer deltaZ = " << deltaZ << endl;

  float zpos = propState.globalPosition().z();
  if ( fabs( zpos) < length + deltaZ) result.push_back( bl);

  if ( fabs( zpos) < length - deltaZ) return true;
  else return false;
}

bool SimpleNavigableLayer::wellInside( const FreeTrajectoryState& fts,
				       PropagationDirection dir,
				       const ForwardDetLayer* fl,
				       DLC& result) const
{
  TSOS propState = propagator(dir).propagate( fts, fl->specificSurface());
  if ( !propState.isValid()) return false;

  if (fl==detLayer()){
    LogDebug("SimpleNavigableLayer")<<"self propagating from:\n"
				    <<fts<<"\n to \n"
				    <<*propState.freeState();
  }

  //if requested avoids crossing over the tracker 
  if (theCheckCrossingSide){
    bool backTobackTransverse = (fts.position().x()*propState.globalPosition().x() + fts.position().y()*propState.globalPosition().y())<0;
    bool backToback = propState.globalPosition().basicVector().dot(fts.position().basicVector())<0;

    if (backTobackTransverse || backToback ){
      LogTrace("TkNavigation") << "Crossing over prevented!\nStaring from (x,y,z,r) (" 
			       << fts.position().x()<<","<< fts.position().y()<<","<< fts.position().z()<<","<<fts.position().perp()
			       << ") going to TSOS (x,y,z,r)" 
			       << propState.globalPosition().x()<<","<< propState.globalPosition().y()<<","<< propState.globalPosition().z()<<","<<propState.globalPosition().perp()<<")";; 
      return false;
    
  //	if (fts.position().z()*propState.globalPosition().z() < 0) return false;
    }}


  float rpos = propState.globalPosition().perp();
  float innerR = fl->specificSurface().innerRadius();
  float outerR = fl->specificSurface().outerRadius();
 
  // take into account the thickness of the layer
  float deltaR = fl->surface().bounds().thickness()/2. *
    fabs( tan( propState.localDirection().theta()));

  // take into account the error on the predicted state
  const float nSigma = theEpsilon;
  if (propState.hasError()) {
    LocalError err = propState.localError().positionError();
    // ignore correlation for the moment...
    deltaR += nSigma * sqrt(err.xx() + err.yy());
  }

  // cout << "SimpleNavigableLayer BarrelDetLayer deltaR = " << deltaR << endl;

  if ( innerR-deltaR < rpos && rpos < outerR+deltaR) result.push_back( fl);
  
  if ( innerR+deltaR < rpos && rpos < outerR-deltaR) return true;
  else return false;
}

bool SimpleNavigableLayer::wellInside( const FreeTrajectoryState& fts,
				       PropagationDirection dir,
				       const DLC& layers,
				       DLC& result) const
{

  // cout << "Entering SimpleNavigableLayer::wellInside" << endl;

  for (DLC::const_iterator i = layers.begin(); i != layers.end(); i++) {
    const BarrelDetLayer* bl = dynamic_cast<const BarrelDetLayer*>(*i);
    if ( bl != 0) {
      if (wellInside( fts, dir, bl, result)) return true;
    }
    else {
      const ForwardDetLayer* fl = dynamic_cast<const ForwardDetLayer*>(*i);
      if ( fl == 0) edm::LogError("TkNavigation") << "dynamic_cast<const ForwardDetLayer*> failed" ;
      if (wellInside( fts, dir, fl, result)) return true;
    }
  }
  return false;
}

bool SimpleNavigableLayer::wellInside( const FreeTrajectoryState& fts,
				       PropagationDirection dir,
				       ConstBDLI begin, ConstBDLI end,
				       DLC& result) const
{
  for ( ConstBDLI i = begin; i < end; i++) {
    if (wellInside( fts, dir, *i, result)) return true;
  }
  return false;
}

bool SimpleNavigableLayer::wellInside( const FreeTrajectoryState& fts,
				       PropagationDirection dir,
				       ConstFDLI begin, ConstFDLI end,
				       DLC& result) const
{
  for ( ConstFDLI i = begin; i < end; i++) {
    if (wellInside( fts, dir, *i, result)) return true;
  }
  return false;
}

Propagator& SimpleNavigableLayer::propagator( PropagationDirection dir) const
{
#ifndef CMS_NO_MUTABLE
  thePropagator.setPropagationDirection(dir);
  return thePropagator;
#else
  SimpleNavigableLayer* mthis = const_cast<SimpleNavigableLayer*>(this);
  mthis->thePropagator.setPropagationDirection(dir);
  return mthis->thePropagator;
#endif
}

void SimpleNavigableLayer::pushResult( DLC& result, const BDLC& tmp) const 
{
  for ( ConstBDLI i = tmp.begin(); i != tmp.end(); i++) {
    result.push_back(*i);
  }
}

void SimpleNavigableLayer::pushResult( DLC& result, const FDLC& tmp) const 
{
  for ( ConstFDLI i = tmp.begin(); i != tmp.end(); i++) {
    result.push_back(*i);
  }
}

std::vector< const DetLayer * > SimpleNavigableLayer::compatibleLayers (const FreeTrajectoryState &fts, 
									PropagationDirection timeDirection,
									int& counter) const {
  typedef std::vector<const DetLayer*> Lvect;
  typedef std::set<const DetLayer *> Lset;  

  Lvect empty,result;
  result.reserve(15); //that's a max
  Lset collect; //a container of unique instances. to avoid duplicates
  
  Lset layerToTry,nextLayerToTry;//set used for iterations
  //initiate the first iteration
  Lvect someLayers(nextLayers(fts,timeDirection));
  layerToTry.insert(someLayers.begin(),someLayers.end());
  
  while (!layerToTry.empty() && (counter++)<=150){
    LogDebug("SimpleNavigableLayer")
      <<counter<<"] going to check on : "<<layerToTry.size()<<" next layers.";
    //clear this set first, it will be swaped with layerToTry
    nextLayerToTry.clear();
    for (Lset::iterator toTry=layerToTry.begin();toTry!=layerToTry.end();++toTry){
      //add the layer you tried.
      LogDebug("SimpleNavigableLayer")
	<<counter<<"] adding layer with pointer: "<<(*toTry)
	<<" first detid: "<<DetIdInfo::info((*toTry)->basicComponents().front()->geographicalId());
      collect.insert( (*toTry) );
      
      //find the next layers from it
      Lvect nextLayers( (*toTry)->nextLayers(fts,timeDirection) );
      LogDebug("SimpleNavigableLayer")
	<<counter<<"] this layer has : "<<nextLayers.size()<<" next layers.";
      for (Lvect::iterator next=nextLayers.begin();next!=nextLayers.end();++next){
	//if the next layer is not already in the collected layers, add it for next round.
	bool alreadyTested=false;
	for (Lset::iterator testMe=collect.begin();testMe!=collect.end();++testMe){ //find method does not work well!
	  if ((*testMe) == (*next))
	    {alreadyTested=true;break;}
	}
	//	if ( collect.find(*next)!=collect.end()){ //find method does not work it seems...
	if (!alreadyTested){
	  nextLayerToTry.insert( *next );
	  LogDebug("SimpleNavigableLayer")
	    <<counter<<"] to try next, layer with pointer: "<<(*next)
	    <<" first detid: "<<DetIdInfo::info((*next)->basicComponents().front()->geographicalId());
	}
	else{
	  LogDebug("SimpleNavigableLayer")
	    <<counter<<"] skipping for next tryout, layer with pointer: "<<(*next)
	    <<" first detid: "<<DetIdInfo::info((*next)->basicComponents().front()->geographicalId());
	  //for (Lset::iterator testMe=collect.begin();testMe!=collect.end();++testMe){  LogTrace("SimpleNavigableLayer") <<"pointer collected:" <<(*testMe); }
	}
      }
      LogDebug("SimpleNavigableLayer")
	<<counter<<"] "<<nextLayerToTry.size()<<" layers to try next so far.";
    }
    //swap now that you where to go next.
    layerToTry.swap(nextLayerToTry);
  }
  if(counter>=150) {
    edm::LogWarning("SimpleNavigableLayer") << "WARNING: compatibleLayers() more than 150 iterations!!! Bailing out..";
    counter = -1;return empty;
  }

  result.insert(result.end(),collect.begin(),collect.end()); //cannot do a swap it seems
  LogDebug("SimpleNavigableLayer")
      <<"Number of compatible layers: "<<result.size();
  return result;


}
