#include "RecoTracker/TkDetLayers/interface/PixelBlade.h"
#include "RecoTracker/TkDetLayers/interface/BladeShapeBuilderFromDet.h"

#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"

#include "Utilities/General/interface/CMSexception.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"


typedef GeometricSearchDet::DetWithState DetWithState;

PixelBlade::PixelBlade(vector<const GeomDet*>& frontDets,
		       vector<const GeomDet*>& backDets):		       
  theFrontDets(frontDets), theBackDets(backDets) 
{
  theDets.assign(theFrontDets.begin(),theFrontDets.end());
  theDets.insert(theDets.end(),theBackDets.begin(),theBackDets.end());

  theDiskSector      = BladeShapeBuilderFromDet()(theDets);  
  theFrontDiskSector = BladeShapeBuilderFromDet()(theFrontDets);
  theBackDiskSector  = BladeShapeBuilderFromDet()(theBackDets);   


  /*--------- DEBUG INFO --------------
  cout << "DEBUG INFO for PixelBlade" << endl;
  cout << "this: " << this << endl;
  cout << "PixelForwardLayer.surfcace.z(): " 
       << this->surface().position().z() << endl;
  cout << "PixelForwardLayer.surfcace.innerR(): " 
       << this->specificSurface().innerRadius() << endl;
  cout << "PixelForwardLayer.surfcace.outerR(): " 
       << this->specificSurface().outerRadius() << endl;
  -----------------------------------*/

}


vector<const GeometricSearchDet*> 
PixelBlade::components() const{
  return vector<const GeometricSearchDet*>();
}

pair<bool, TrajectoryStateOnSurface>
PixelBlade::compatible( const TrajectoryStateOnSurface& ts, const Propagator&, 
			const MeasurementEstimator&) const{
  cout << "temporary dummy implementation of PixelBlade::compatible()!!" << endl;
  return pair<bool,TrajectoryStateOnSurface>();
}


vector<DetWithState> 
PixelBlade::compatibleDets( const TrajectoryStateOnSurface& startingState,
			    const Propagator& prop, 
			    const MeasurementEstimator& est) const{

  // standard implementation of compatibleDets() for class which have 
  // groupedCompatibleDets implemented.
  // This code should be moved in a common place intead of being 
  // copied many times.
  
  vector<DetWithState> result;  
  vector<DetGroup> vectorGroups = groupedCompatibleDets(startingState,prop,est);
  for(vector<DetGroup>::const_iterator itDG=vectorGroups.begin();
      itDG!=vectorGroups.end();itDG++){
    for(vector<DetGroupElement>::const_iterator itDGE=itDG->begin();
	itDGE!=itDG->end();itDGE++){
      result.push_back(DetWithState(itDGE->det(),itDGE->trajectoryState()));
    }
  }
  return result;  
}



vector<DetGroup> 
PixelBlade::groupedCompatibleDets( const TrajectoryStateOnSurface& tsos,
				   const Propagator& prop,
				   const MeasurementEstimator& est) const
{
  vector<DetGroup> closestResult;
  SubLayerCrossings  crossings; 
  try{
    crossings = computeCrossings( tsos, prop.propagationDirection());  
  }
  catch(Genexception& err){ //In ORCA, it was a DetLogicError exception
    cout << "Aie, got an exception in PixelBlade::groupedCompatibleDets:" 
	 << err.what() << endl;
    return closestResult;
  }    
  addClosest( tsos, prop, est, crossings.closest(), closestResult);

  if (closestResult.empty()){
    vector<DetGroup> nextResult;
    addClosest( tsos, prop, est, crossings.other(), nextResult);
    if(nextResult.empty())    return nextResult;

    DetGroupElement nextGel( nextResult.front().front());  
    int crossingSide = LayerCrossingSide().barrelSide( nextGel.trajectoryState(), prop);
    DetGroupMerger merger;
    return  merger.orderAndMergeTwoLevels( closestResult, nextResult, 
					   crossings.closestIndex(), crossingSide);   
  }
  
  DetGroupElement closestGel( closestResult.front().front());
  float window = computeWindowSize( closestGel.det(), closestGel.trajectoryState(), est);

  searchNeighbors( tsos, prop, est, crossings.closest(), window,
		   closestResult, false);

  vector<DetGroup> nextResult;
  searchNeighbors( tsos, prop, est, crossings.other(), window,
		   nextResult, true);

  int crossingSide = LayerCrossingSide().barrelSide( closestGel.trajectoryState(), prop);
  DetGroupMerger merger;
  return merger.orderAndMergeTwoLevels( closestResult, nextResult, 
					crossings.closestIndex(), crossingSide);

}


SubLayerCrossings 
PixelBlade::computeCrossings( const TrajectoryStateOnSurface& startingState,
			      PropagationDirection propDir) const
{
  HelixPlaneCrossing::PositionType startPos( startingState.globalPosition());
  HelixPlaneCrossing::DirectionType startDir( startingState.globalMomentum());
  double rho( startingState.transverseCurvature());

  HelixArbitraryPlaneCrossing crossing( startPos, startDir, rho, propDir);

  pair<bool,double> innerPath = crossing.pathLength( *theFrontDiskSector);
  if (!innerPath.first) {
    cout << "ERROR in PixelBlade: inner subBlade not crossed by track" << endl;
    //throw DetLogicError("PixelBlade: inner subBlade not crossed by track");
    throw Genexception("PixelBlade: inner subBlade not crossed by track");
  }
  GlobalPoint gInnerPoint( crossing.position(innerPath.second));
  //Code for use of binfinder
  //int innerIndex = theInnerBinFinder.binIndex(gInnerPoint.perp());  
  //float innerDist = fabs( theInnerBinFinder.binPosition(innerIndex) - gInnerPoint.z());
  int innerIndex = findBin(gInnerPoint.perp(),0);
  float innerDist = fabs( findPosition(innerIndex,0).perp() - gInnerPoint.perp());
  SubLayerCrossing innerSLC( 0, innerIndex, gInnerPoint);

  pair<bool,double> outerPath = crossing.pathLength( *theBackDiskSector);
  if (!outerPath.first) {
    cout << "ERROR in PixelBlade: outer subBlade not crossed by track" << endl;
    //throw DetLogicError("PixelBlade: outer subBlade not crossed by track");
    throw Genexception("PixelBlade: outer subBlade not crossed by track");
  }
  GlobalPoint gOuterPoint( crossing.position(outerPath.second));
  //Code for use of binfinder
  //int outerIndex = theOuterBinFinder.binIndex(gOuterPoint.perp());
  //float outerDist = fabs( theOuterBinFinder.binPosition(outerIndex) - gOuterPoint.perp());
  int outerIndex  = findBin(gOuterPoint.perp(),1);
  float outerDist = fabs( findPosition(outerIndex,1).perp() - gOuterPoint.perp());
  SubLayerCrossing outerSLC( 1, outerIndex, gOuterPoint);

  if (innerDist < outerDist) {
    return SubLayerCrossings( innerSLC, outerSLC, 0);
  }
  else {
    return SubLayerCrossings( outerSLC, innerSLC, 1);
  } 
}




bool 
PixelBlade::addClosest( const TrajectoryStateOnSurface& tsos,
			const Propagator& prop,
			const MeasurementEstimator& est,
			const SubLayerCrossing& crossing,
			vector<DetGroup>& result) const
{

  const vector<const GeomDet*>& sBlade( subBlade( crossing.subLayerIndex()));
  return CompatibleDetToGroupAdder().add( *sBlade[crossing.closestDetIndex()], 
					  tsos, prop, est, result);
}


float PixelBlade::computeWindowSize( const GeomDet* det, 
				     const TrajectoryStateOnSurface& tsos, 
				     const MeasurementEstimator& est) const
{
  return
    est.maximalLocalDisplacement(tsos, dynamic_cast<const BoundPlane&>(det->surface())).x();
}




void PixelBlade::searchNeighbors( const TrajectoryStateOnSurface& tsos,
				  const Propagator& prop,
				  const MeasurementEstimator& est,
				  const SubLayerCrossing& crossing,
				  float window, 
				  vector<DetGroup>& result,
				  bool checkClosest) const
{
  GlobalPoint gCrossingPos = crossing.position();

  const vector<const GeomDet*>& sBlade( subBlade( crossing.subLayerIndex()));
 
  int closestIndex = crossing.closestDetIndex();
  int negStartIndex = closestIndex-1;
  int posStartIndex = closestIndex+1;

  if (checkClosest) { // must decide if the closest is on the neg or pos side
    if (gCrossingPos.perp() < sBlade[closestIndex]->surface().position().perp()) {
      posStartIndex = closestIndex;
    }
    else {
      negStartIndex = closestIndex;
    }
  }

  CompatibleDetToGroupAdder adder;
  for (int idet=negStartIndex; idet >= 0; idet--) {
    if (!overlap( gCrossingPos, *sBlade[idet], window)) break;
    if (!adder.add( *sBlade[idet], tsos, prop, est, result)) break;
  }
  for (int idet=posStartIndex; idet < static_cast<int>(sBlade.size()); idet++) {
    if (!overlap( gCrossingPos, *sBlade[idet], window)) break;
    if (!adder.add( *sBlade[idet], tsos, prop, est, result)) break;
  }
}



bool PixelBlade::overlap( const GlobalPoint& crossPoint, const GeomDet& det, float window) const
{
  // check if the z window around TSOS overlaps with the detector theDet (with a 1% margin added)
  
  //   const float tolerance = 0.1;
  const float relativeMargin = 1.01;

  LocalPoint localCrossPoint( det.surface().toLocal(crossPoint));
  //   if (fabs(localCrossPoint.z()) > tolerance) {
  //     cout << "PixelBlade::overlap calculation assumes point on surface, but it is off by "
  // 	 << localCrossPoint.z() << endl;
  //   }

  float localX = localCrossPoint.x();
  float detHalfLength = det.surface().bounds().length()/2.;

  //   cout << "PixelBlade::overlap: Det at " << det.position() << " hit at " << localY 
  //        << " Window " << window << " halflength "  << detHalfLength << endl;
  
  if ( ( fabs(localX)-window) < relativeMargin*detHalfLength ) { // FIXME: margin hard-wired!
    return true;
  } else {
    return false;
  }
}

int 
PixelBlade::findBin( float R,int diskSectorIndex) const 
{
  vector<const GeomDet*> localDets = diskSectorIndex==0 ? theFrontDets : theBackDets;
  
  int theBin = -1;
  float rDiff = 200.;
  for (vector<const GeomDet*>::const_iterator i=localDets.begin(); i !=localDets.end(); i++){
    float testDiff = fabs( R - (**i).surface().position().perp());
    if ( testDiff < rDiff) {
      rDiff = testDiff;
      theBin = i - localDets.begin();
    }
  }
  return theBin;
}



GlobalPoint 
PixelBlade::findPosition(int index,int diskSectorType) const 
{
  vector<const GeomDet*> diskSector = diskSectorType == 0 ? theFrontDets : theBackDets; 
  return (diskSector[index])->surface().position();
}

