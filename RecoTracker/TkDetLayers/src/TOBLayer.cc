#include "RecoTracker/TkDetLayers/interface/TOBLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"
#include "RecoTracker/TkDetLayers/interface/GlobalDetRodRangeZPhi.h"

#include "TrackingTools/DetLayers/interface/DetLayerException.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/HelixBarrelCylinderCrossing.h"
#include "TrackingTools/DetLayers/interface/CylinderBuilderFromDet.h"
#include "TrackingTools/DetLayers/interface/PhiLess.h"
#include "TrackingTools/DetLayers/interface/rangesIntersect.h"

using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;

TOBLayer::TOBLayer(vector<const TOBRod*>& innerRods,
		   vector<const TOBRod*>& outerRods) : 
  theInnerComps(innerRods.begin(),innerRods.end()), 
  theOuterComps(outerRods.begin(),outerRods.end())
{
  theComps.assign(theInnerComps.begin(),theInnerComps.end());
  theComps.insert(theComps.end(),theOuterComps.begin(),theOuterComps.end());
    
  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin();
      it!=theComps.end();it++){  
    theBasicComps.insert(theBasicComps.end(),	
			 (**it).basicComponents().begin(),
			 (**it).basicComponents().end());
  }

  
  theInnerCylinder = cylinder( theInnerComps);
  theOuterCylinder = cylinder( theOuterComps);

  theInnerBinFinder = BinFinderType(theInnerComps.front()->position().phi(),
				    theInnerComps.size());
  theOuterBinFinder = BinFinderType(theOuterComps.front()->position().phi(),
				    theOuterComps.size());

  
  BarrelDetLayer::initialize();

   //--------- DEBUG INFO --------------
  LogDebug("TkDetLayers") << "==== DEBUG TOBLayer =====" ; 
  LogDebug("TkDetLayers") << "innerCyl radius, thickness, lenght: " 
			  << theInnerCylinder->radius() << " , "
			  << theInnerCylinder->bounds().thickness() << " , "
			  << theInnerCylinder->bounds().length() ;
 
  LogDebug("TkDetLayers") << "outerCyl radius, thickness, lenght: " 
			  << theOuterCylinder->radius() << " , "
			  << theOuterCylinder->bounds().thickness() << " , "
			  << theOuterCylinder->bounds().length() ;

  LogDebug("TkDetLayers") << "Cyl radius, thickness, lenght: " 
			  << specificSurface().radius() << " , "
			  << specificSurface().bounds().thickness() << " , "
			  << specificSurface().bounds().length() ;

  for (vector<const GeometricSearchDet*>::const_iterator i=theInnerComps.begin();
       i != theInnerComps.end(); i++){
     LogDebug("TkDetLayers") << "inner TOBRod pos z,perp,eta,phi: " 
			     << (**i).position().z() << " , " 
			     << (**i).position().perp() << " , " 
			     << (**i).position().eta() << " , " 
			     << (**i).position().phi() ;
  }
  
  for (vector<const GeometricSearchDet*>::const_iterator i=theOuterComps.begin();
       i != theOuterComps.end(); i++){
    LogDebug("TkDetLayers") << "outer TOBRod pos z,perp,eta,phi: " 
			    << (**i).position().z() << " , " 
			    << (**i).position().perp() << " , " 
			    << (**i).position().eta() << " , " 
			    << (**i).position().phi() ;
  }
  LogDebug("TkDetLayers") << "==== end DEBUG TOBLayer =====" ;   
}


TOBLayer::~TOBLayer(){
  vector<const GeometricSearchDet*>::const_iterator i;
  for (i=theComps.begin(); i!=theComps.end(); i++) {
    delete *i;
  }
} 


  
vector<DetWithState> 
TOBLayer::compatibleDets( const TrajectoryStateOnSurface& startingState,
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
TOBLayer::groupedCompatibleDets( const TrajectoryStateOnSurface& tsos,
				 const Propagator& prop,
				 const MeasurementEstimator& est) const
{
  vector<DetGroup> closestResult;
  SubLayerCrossings crossings;
  crossings = computeCrossings( tsos, prop.propagationDirection());
  if(! crossings.isValid()) return closestResult;

  addClosest( tsos, prop, est, crossings.closest(), closestResult);
  if (closestResult.empty()){
    vector<DetGroup> nextResult;
    addClosest( tsos, prop, est, crossings.other(), nextResult);
    return nextResult;
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



// private methods for the implementation of groupedCompatibleDets()

SubLayerCrossings TOBLayer::computeCrossings( const TrajectoryStateOnSurface& startingState,
					      PropagationDirection propDir) const
{
  GlobalPoint startPos( startingState.globalPosition());
  GlobalVector startDir( startingState.globalMomentum());
  double rho( startingState.transverseCurvature());

  HelixBarrelCylinderCrossing innerCrossing( startPos, startDir, rho,
					     propDir,*theInnerCylinder);
  if (!innerCrossing.hasSolution()) return SubLayerCrossings();

  GlobalPoint gInnerPoint( innerCrossing.position());
  int innerIndex = theInnerBinFinder.binIndex(gInnerPoint.phi());
  float innerDist = theInnerBinFinder.binPosition(innerIndex) - gInnerPoint.phi();
  SubLayerCrossing innerSLC( 0, innerIndex, gInnerPoint);

  HelixBarrelCylinderCrossing outerCrossing( startPos, startDir, rho,
					     propDir,*theOuterCylinder);
  if (!outerCrossing.hasSolution()) return SubLayerCrossings();

  GlobalPoint gOuterPoint( outerCrossing.position());
  int outerIndex = theOuterBinFinder.binIndex(gOuterPoint.phi());
  float outerDist = theOuterBinFinder.binPosition(outerIndex) - gOuterPoint.phi() ;
  SubLayerCrossing outerSLC( 1, outerIndex, gOuterPoint);
  
  innerDist *= PhiLess()( theInnerBinFinder.binPosition(innerIndex),gInnerPoint.phi()) ? -1. : 1.; 
  outerDist *= PhiLess()( theOuterBinFinder.binPosition(outerIndex),gOuterPoint.phi()) ? -1. : 1.; 
  if (innerDist < 0.) { innerDist += 2.*Geom::pi();}
  if (outerDist < 0.) { outerDist += 2.*Geom::pi();}
  

  if (innerDist < outerDist) {
    return SubLayerCrossings( innerSLC, outerSLC, 0);
  }
  else {
    return SubLayerCrossings( outerSLC, innerSLC, 1);
  } 
}

bool TOBLayer::addClosest( const TrajectoryStateOnSurface& tsos,
			   const Propagator& prop,
			   const MeasurementEstimator& est,
			   const SubLayerCrossing& crossing,
			   vector<DetGroup>& result) const
{
  const vector<const GeometricSearchDet*>& sub( subLayer( crossing.subLayerIndex()));
  const GeometricSearchDet* det(sub[crossing.closestDetIndex()]);
  return CompatibleDetToGroupAdder().add( *det, tsos, prop, est, result);
}

float TOBLayer::computeWindowSize( const GeomDet* det, 
				   const TrajectoryStateOnSurface& tsos, 
				   const MeasurementEstimator& est) const
{
  double xmax = 
    est.maximalLocalDisplacement(tsos, dynamic_cast<const BoundPlane&>(det->surface())).x();
  return calculatePhiWindow( xmax, *det, tsos);
}


double TOBLayer::calculatePhiWindow( double Xmax, const GeomDet& det,
				     const TrajectoryStateOnSurface& state) const
{

  LocalPoint startPoint = state.localPosition();
  LocalVector shift( Xmax , 0. , 0.);
  LocalPoint shift1 = startPoint + shift;
  LocalPoint shift2 = startPoint + (-shift); 
  //LocalPoint shift2( startPoint); //original code;
  //shift2 -= shift;

  double phi1 = det.surface().toGlobal(shift1).phi();
  double phi2 = det.surface().toGlobal(shift2).phi();
  double phiStart = state.globalPosition().phi();
  double phiWin = min(fabs(phiStart-phi1),fabs(phiStart-phi2));

  return phiWin;
}


void TOBLayer::searchNeighbors( const TrajectoryStateOnSurface& tsos,
				const Propagator& prop,
				const MeasurementEstimator& est,
				const SubLayerCrossing& crossing,
				float window, 
				vector<DetGroup>& result,
				bool checkClosest) const
{
  GlobalPoint gCrossingPos = crossing.position();

  const vector<const GeometricSearchDet*>& sLayer( subLayer( crossing.subLayerIndex()));
 
  int closestIndex = crossing.closestDetIndex();
  int negStartIndex = closestIndex-1;
  int posStartIndex = closestIndex+1;

  if (checkClosest) { // must decide if the closest is on the neg or pos side
    if ( PhiLess()( gCrossingPos.phi(), sLayer[closestIndex]->position().phi())) {
      posStartIndex = closestIndex;
    }
    else {
      negStartIndex = closestIndex;
    }
  }

  const BinFinderType& binFinder = (crossing.subLayerIndex()==0 ? theInnerBinFinder : theOuterBinFinder);

  CompatibleDetToGroupAdder adder;
  int quarter = sLayer.size()/4;
  for (int idet=negStartIndex; idet >= negStartIndex - quarter; idet--) {
    const GeometricSearchDet* neighborRod = sLayer[binFinder.binIndex(idet)];
    if (!overlap( gCrossingPos, *neighborRod, window)) break;
    if (!adder.add( *neighborRod, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
  for (int idet=posStartIndex; idet < posStartIndex + quarter; idet++) {
    const GeometricSearchDet* neighborRod = sLayer[binFinder.binIndex(idet)];
    if (!overlap( gCrossingPos, *neighborRod, window)) break;
    if (!adder.add( *neighborRod, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
}

bool TOBLayer::overlap( const GlobalPoint& gpos, const GeometricSearchDet& gsdet, float phiWin) const
{
  GlobalPoint crossPoint(gpos);

  // introduce offset (extrapolated point and true propagated point differ by 0.0003 - 0.00033, 
  // due to thickness of Rod of 1 cm) 
  const float phiOffset = 0.00034;  //...TOBE CHECKED LATER...
  phiWin += phiOffset;

  // detector phi range
  const TOBRod& theRod = dynamic_cast<const TOBRod&>(gsdet);
  GlobalDetRodRangeZPhi rodRange( theRod.specificSurface());
  pair<float,float> phiRange(crossPoint.phi()-phiWin, crossPoint.phi()+phiWin);

  //   // debug
  //   edm::LogInfo(TkDetLayers) ;
  //   edm::LogInfo(TkDetLayers) << " overlapInPhi: position, det phi range " 
  //        << "("<< rod.position().perp() << ", " << rod.position().phi() << ")  "
  //        << rodRange.phiRange().first << " " << rodRange.phiRange().second ;
  //   edm::LogInfo(TkDetLayers) << " overlapInPhi: cross point phi, window " << crossPoint.phi() << " " << phiWin ;
  //   edm::LogInfo(TkDetLayers) << " overlapInPhi: search window: " << crossPoint.phi()-phiWin << "  " << crossPoint.phi()+phiWin ;

  if ( rangesIntersect(phiRange, rodRange.phiRange(), PhiLess())) {
    return true;
  } else {
    return false;
  }
} 


BoundCylinder* TOBLayer::cylinder( const vector<const GeometricSearchDet*>& rods) const 
{
  vector<const GeomDet*> tmp;
  for (vector<const GeometricSearchDet*>::const_iterator it=rods.begin(); it!=rods.end(); it++) {    
    tmp.insert(tmp.end(),(*it)->basicComponents().begin(),(*it)->basicComponents().end());
  }
  return CylinderBuilderFromDet()( tmp.begin(), tmp.end());
}
