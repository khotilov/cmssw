#include "RecoTracker/TkDetLayers/interface/TIBLayer.h"

#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"

#include "TrackingTools/DetLayers/interface/DetLayerException.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/HelixBarrelCylinderCrossing.h"
#include "Geometry/Surface/interface/SimpleCylinderBounds.h"
#include "TrackingTools/DetLayers/src/DetLessZ.h"


using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;

TIBLayer::TIBLayer(vector<const TIBRing*>& innerRings,
		   vector<const TIBRing*>& outerRings) : 
  theInnerComps(innerRings.begin(),innerRings.end()), 
  theOuterComps(outerRings.begin(),outerRings.end())
{
  theComps.assign(theInnerComps.begin(),theInnerComps.end());
  theComps.insert(theComps.end(),theOuterComps.begin(),theOuterComps.end());
  
  sort(theComps.begin(),theComps.end(),DetLessZ());
  sort(theInnerComps.begin(),theInnerComps.end(),DetLessZ());
  sort(theOuterComps.begin(),theOuterComps.end(),DetLessZ());
  
  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin();
      it!=theComps.end();it++){  
    theBasicComps.insert(theBasicComps.end(),	
			 (**it).basicComponents().begin(),
			 (**it).basicComponents().end());
  }

  // initialize the surface
  theInnerCylinder = cylinder( theInnerComps);
  theOuterCylinder = cylinder( theOuterComps);
  setSurface( cylinder(theComps) );

  /*
  cout << "==== DEBUG TIBLayer =====" << endl; 
  for (vector<const TIBRing*>::const_iterator i=theInnerRings.begin();
       i != theInnerRings.end(); i++){
    cout << "inner TIBRing pos z,radius,eta,phi: " 
	 << (**i).position().z() << " , " 
	 << (**i).specificSurface().radius() << " , " 
	 << (**i).position().eta() << " , " 
	 << (**i).position().phi() << endl;
  }

  for (vector<const TIBRing*>::const_iterator i=theOuterRings.begin();
       i != theOuterRings.end(); i++){
    cout << "outer TIBRing pos z,radius,eta,phi: " 
	 << (**i).position().z() << " , " 
	 << (**i).specificSurface().radius() << " , " 
	 << (**i).position().eta() << " , " 
	 << (**i).position().phi() << endl;
  }
  */


  // initialise the bin finders
  //  vector<const GeometricSearchDet*> tmpIn;
  //for (vector<const TIBRing*>::const_iterator i=theInnerRings.begin();
  //     i != theInnerRings.end(); i++) tmpIn.push_back(*i);
  theInnerBinFinder = GeneralBinFinderInZforGeometricSearchDet<float>(theInnerComps.begin(), 
								      theInnerComps.end());
 
  theOuterBinFinder = GeneralBinFinderInZforGeometricSearchDet<float>(theOuterComps.begin(),
								      theOuterComps.end());
}

TIBLayer::~TIBLayer(){
  vector<const GeometricSearchDet*>::const_iterator i;
  for (i=theComps.begin(); i!=theComps.end(); i++) {
    delete *i;
  }
} 

  
vector<DetWithState> 
TIBLayer::compatibleDets( const TrajectoryStateOnSurface& startingState,
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


BoundCylinder* 
TIBLayer::cylinder( const vector<const GeometricSearchDet*>& rings)
{
  float leftPos = rings.front()->surface().position().z();
  float rightPos = rings.back()->surface().position().z();

  const TIBRing* frontRing = dynamic_cast<const TIBRing*>(rings.front());
  const TIBRing* backRing  = dynamic_cast<const TIBRing*>(rings.back());
  float r = frontRing->specificSurface().radius(); 
  const Bounds& leftBounds  = frontRing->specificSurface().bounds();
  const Bounds& rightBounds = backRing->specificSurface().bounds();

  //float r = rings.front()->specificSurface().radius();
  //const Bounds& leftBounds = rings.front()->specificSurface().bounds();
  //const Bounds& rightBounds = rings.back()->specificSurface().bounds();

  float thick = leftBounds.thickness() / 2;
  float zmin = leftPos  - leftBounds.length() / 2;
  float zmax = rightPos + rightBounds.length() / 2;
  float zpos = 0.5*(leftPos+rightPos);

  return new BoundCylinder( Surface::PositionType( 0, 0, zpos), 
			    rings.front()->surface().rotation(),
			    SimpleCylinderBounds( r-thick, r+thick,
						  zmin-zpos, zmax-zpos));
}



vector<DetGroup> 
TIBLayer::groupedCompatibleDets( const TrajectoryStateOnSurface& tsos,
				 const Propagator& prop,
				 const MeasurementEstimator& est) const{
  
  vector<DetGroup> closestResult;
  SubLayerCrossings  crossings; 
  try{
    crossings = computeCrossings( tsos, prop.propagationDirection());  
  }
  catch(DetLayerException& err){
    //cout << "Aie, got a DetLayerException in TIBLayer::groupedCompatibleDets:" 
    //	 << err.what() << endl;
    return closestResult;
  }

  

  addClosest( tsos, prop, est, crossings.closest(), closestResult);
  // this differs from compatibleDets logic, which checks next in such cases!!!
  if (closestResult.empty())    return closestResult;
  
  

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

SubLayerCrossings TIBLayer::computeCrossings( const TrajectoryStateOnSurface& startingState,
						      PropagationDirection propDir) const
{
  GlobalPoint startPos( startingState.globalPosition());
  GlobalVector startDir( startingState.globalMomentum());
  double rho( startingState.transverseCurvature());

  HelixBarrelCylinderCrossing innerCrossing( startPos, startDir, rho,
					     propDir,*theInnerCylinder);
  if (!innerCrossing.hasSolution()) {
    //cout << "ERROR in TIBLayer: inner cylinder not crossed by track" << endl;
    throw DetLayerException("TIBLayer: inner cylinder not crossed by track");
  }

  GlobalPoint gInnerPoint( innerCrossing.position());
  int innerIndex = theInnerBinFinder.binIndex(gInnerPoint.z());
  const GeometricSearchDet* innerRing( theInnerComps[innerIndex]);
  float innerDist = fabs( innerRing->surface().position().z() - gInnerPoint.z());
  SubLayerCrossing innerSLC( 0, innerIndex, gInnerPoint);

  HelixBarrelCylinderCrossing outerCrossing( startPos, startDir, rho,
					     propDir,*theOuterCylinder);
  if (!outerCrossing.hasSolution()) {
    //cout << "ERROR in TIBLayer: outer cylinder not crossed by track" << endl;
    throw DetLayerException("TIBLayer: inner cylinder not crossed by track");
  }

  GlobalPoint gOuterPoint( outerCrossing.position());
  int outerIndex = theOuterBinFinder.binIndex(gOuterPoint.z());
  const GeometricSearchDet* outerRing( theOuterComps[outerIndex]);
  float outerDist = fabs( outerRing->surface().position().z() - gOuterPoint.z());
  SubLayerCrossing outerSLC( 1, outerIndex, gOuterPoint);

  if (innerDist < outerDist) {
    return SubLayerCrossings( innerSLC, outerSLC, 0);
  }
  else {
    return SubLayerCrossings( outerSLC, innerSLC, 1);
  } 
}

bool TIBLayer::addClosest( const TrajectoryStateOnSurface& tsos,
				      const Propagator& prop,
				      const MeasurementEstimator& est,
				      const SubLayerCrossing& crossing,
				      vector<DetGroup>& result) const
{
//   cout << "Entering TIBLayer::addClosest" << endl;

  const vector<const GeometricSearchDet*>& sub( subLayer( crossing.subLayerIndex()));
  const Det* det(sub[crossing.closestDetIndex()]);
  return CompatibleDetToGroupAdder().add( *det, tsos, prop, est, result);
}

void TIBLayer::searchNeighbors( const TrajectoryStateOnSurface& tsos,
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
    if (gCrossingPos.z() < sLayer[closestIndex]->surface().position().z()) {
      posStartIndex = closestIndex;
    }
    else {
      negStartIndex = closestIndex;
    }
  }

  CompatibleDetToGroupAdder adder;
  for (int idet=negStartIndex; idet >= 0; idet--) {
    const GeometricSearchDet* neighborRing = sLayer[idet];
    if (!overlap( gCrossingPos, *neighborRing, window)) break;
    if (!adder.add( *neighborRing, tsos, prop, est, result)) break;
  }
  for (int idet=posStartIndex; idet < static_cast<int>(sLayer.size()); idet++) {
    const GeometricSearchDet* neighborRing = sLayer[idet];
    if (!overlap( gCrossingPos, *neighborRing, window)) break;
    if (!adder.add( *neighborRing, tsos, prop, est, result)) break;
  }
}

bool TIBLayer::overlap( const GlobalPoint& crossPoint,
			const GeometricSearchDet& det, 
			float window) const
{
  float halfLength = det.surface().bounds().length()/2.;

//   cout << " TIBLayer: checking ring with z " << det.position().z();

  if ( fabs( crossPoint.z()-det.position().z()) < (halfLength + window)) {
//     cout << "    PASSED" << endl;
    return true;
  } else {
//     cout << "    FAILED " << endl;
    return false;
  }
}

float TIBLayer::computeWindowSize( const GeomDet* det, 
				   const TrajectoryStateOnSurface& tsos, 
				   const MeasurementEstimator& est) const
{
  // we assume the outer and inner rings have about same thickness...

  // improve interface to avoid dynamic_cast?
  const BoundPlane& plane( dynamic_cast<const BoundPlane&>(tsos.surface()));

//   cout << "TIBLayer::computeWindowSize: Y axis of tangent plane is"
//        << plane.toGlobal( LocalVector(0,1,0)) << endl;

  MeasurementEstimator::Local2DVector localError( est.maximalLocalDisplacement(tsos, plane));
  float yError = localError.y();

  float tanTheta = tan( tsos.globalMomentum().theta());
  float thickCorrection = det->surface().bounds().thickness() / (2.*fabs( tanTheta));

  // FIXME: correct this in case of wide phi window !  

  return yError + thickCorrection;
}



