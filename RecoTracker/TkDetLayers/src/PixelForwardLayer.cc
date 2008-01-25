#include "RecoTracker/TkDetLayers/interface/PixelForwardLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GeometrySurface/interface/BoundingBox.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

#include "TrackingTools/DetLayers/interface/simple_stat.h"
#include "TrackingTools/DetLayers/interface/PhiLess.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing2Order.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"


#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"

using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;

PixelForwardLayer::PixelForwardLayer(vector<const PixelBlade*>& blades):
  theComps(blades.begin(),blades.end())
{
  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin();
      it!=theComps.end();it++){  
    theBasicComps.insert(theBasicComps.end(),	
			 (**it).basicComponents().begin(),
			 (**it).basicComponents().end());
  }

  //They should be already phi-ordered. TO BE CHECKED!!
  //sort( theBlades.begin(), theBlades.end(), PhiLess());
  setSurface( computeSurface() );
  
  //Is a "periodic" binFinderInPhi enough?. TO BE CHECKED!!
  theBinFinder = BinFinderType( theComps.front()->surface().position().phi(),
				theComps.size());

  //--------- DEBUG INFO --------------
  LogDebug("TkDetLayers") << "DEBUG INFO for PixelForwardLayer" << "\n"
			  << "PixelForwardLayer.surfcace.phi(): " 
			  << this->surface().position().phi() << "\n"
			  << "PixelForwardLayer.surfcace.z(): " 
			  << this->surface().position().z() << "\n"
			  << "PixelForwardLayer.surfcace.innerR(): " 
			  << this->specificSurface().innerRadius() << "\n"
			  << "PixelForwardLayer.surfcace.outerR(): " 
			  << this->specificSurface().outerRadius() ;

  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin(); 
      it!=theComps.end(); it++){
    LogDebug("TkDetLayers") << "blades phi,z,r: " 
			    << (*it)->surface().position().phi() << " , "
			    << (*it)->surface().position().z() <<   " , "
			    << (*it)->surface().position().perp();
  }
  //-----------------------------------

    
}

PixelForwardLayer::~PixelForwardLayer(){
  vector<const GeometricSearchDet*>::const_iterator i;
  for (i=theComps.begin(); i!=theComps.end(); i++) {
    delete *i;
  }
} 

void
PixelForwardLayer::groupedCompatibleDetsV( const TrajectoryStateOnSurface& tsos,
					  const Propagator& prop,
					   const MeasurementEstimator& est,
					   std::vector<DetGroup> & result) const {
  vector<DetGroup> closestResult;
  SubTurbineCrossings  crossings; 

  crossings = computeCrossings( tsos, prop.propagationDirection());
  if (!crossings.isValid){
    //edm::LogInfo("TkDetLayers") << "computeCrossings returns invalid in PixelForwardLayer::groupedCompatibleDets:";
    return;
  }

  typedef CompatibleDetToGroupAdder Adder;
  Adder::add( *theComps[theBinFinder.binIndex(crossings.closestIndex)], 
	     tsos, prop, est, closestResult);

  if(closestResult.empty()){
    Adder::add( *theComps[theBinFinder.binIndex(crossings.nextIndex)], 
	       tsos, prop, est, result);
    return;
  }      

  DetGroupElement closestGel( closestResult.front().front());
  float window = computeWindowSize( closestGel.det(), closestGel.trajectoryState(), est);

  float detWidth = closestGel.det()->surface().bounds().width();
  if (crossings.nextDistance < detWidth + window) {
    vector<DetGroup> nextResult;
    if (Adder::add( *theComps[theBinFinder.binIndex(crossings.nextIndex)], 
		   tsos, prop, est, nextResult)) {
      int crossingSide = LayerCrossingSide().endcapSide( tsos, prop);
      int theHelicity = computeHelicity(theComps[theBinFinder.binIndex(crossings.closestIndex)],
					theComps[theBinFinder.binIndex(crossings.nextIndex)] );
      DetGroupMerger::orderAndMergeTwoLevels( closestResult, nextResult, result, 
					      theHelicity, crossingSide);
    }
    else {
      result.swap(closestResult);
    }
  }
  else {
    result.swap(closestResult);
  }

  // only loop over neighbors (other than closest and next) if window is BIG
  if (window > 0.5*detWidth) {
    searchNeighbors( tsos, prop, est, crossings, window, result);
  } 
}



void 
PixelForwardLayer::searchNeighbors( const TrajectoryStateOnSurface& tsos,
				    const Propagator& prop,
				    const MeasurementEstimator& est,
				    const SubTurbineCrossings& crossings,
				    float window, 
				    vector<DetGroup>& result) const
{
   typedef CompatibleDetToGroupAdder Adder;
  int crossingSide = LayerCrossingSide().endcapSide( tsos, prop);
  typedef DetGroupMerger Merger;

  int negStart = min( crossings.closestIndex, crossings.nextIndex) - 1;
  int posStart = max( crossings.closestIndex, crossings.nextIndex) + 1;

  int quarter = theComps.size()/4;
 
  vector<DetGroup> tmp;
  vector<DetGroup> newResult;
  for (int idet=negStart; idet >= negStart - quarter+1; idet--) {
    tmp.clear();
    newResult.clear();
    const GeometricSearchDet* neighbor = theComps[theBinFinder.binIndex(idet)];
    // if (!overlap( gCrossingPos, *neighbor, window)) break; // mybe not needed?
    // maybe also add shallow crossing angle test here???
    if (!Adder::add( *neighbor, tsos, prop, est, tmp)) break;
    int theHelicity = computeHelicity(theComps[theBinFinder.binIndex(idet)],
				      theComps[theBinFinder.binIndex(idet+1)] );
    Merger::orderAndMergeTwoLevels( tmp, result, newResult, theHelicity, crossingSide);
    result.swap(newResult);
  }
  for (int idet=posStart; idet < posStart + quarter-1; idet++) {
    tmp.clear();
    newResult.clear();
    const GeometricSearchDet* neighbor = theComps[theBinFinder.binIndex(idet)];
    // if (!overlap( gCrossingPos, *neighbor, window)) break; // mybe not needed?
    // maybe also add shallow crossing angle test here???
    if (!Adder::add( *neighbor, tsos, prop, est, tmp)) break;
    int theHelicity = computeHelicity(theComps[theBinFinder.binIndex(idet-1)],
				      theComps[theBinFinder.binIndex(idet)] );
    Merger::orderAndMergeTwoLevels( result, tmp, newResult, theHelicity, crossingSide);
    result.swap(newResult);
  }
}

int 
PixelForwardLayer::computeHelicity(const GeometricSearchDet* firstBlade,const GeometricSearchDet* secondBlade) const
{  
  if( fabs(firstBlade->position().z()) < fabs(secondBlade->position().z()) ) return 0;
  return 1;
}

PixelForwardLayer::SubTurbineCrossings 
PixelForwardLayer::computeCrossings( const TrajectoryStateOnSurface& startingState,
				     PropagationDirection propDir) const
{  
  typedef MeasurementEstimator::Local2DVector Local2DVector;

  HelixPlaneCrossing::PositionType startPos( startingState.globalPosition());
  HelixPlaneCrossing::DirectionType startDir( startingState.globalMomentum());
  float rho( startingState.transverseCurvature());

  HelixArbitraryPlaneCrossing turbineCrossing( startPos, startDir, rho,
					       propDir);

  pair<bool,double> thePath = turbineCrossing.pathLength( specificSurface() );
  
  if (!thePath.first) {
    //edm::LogInfo("TkDetLayers") << "ERROR in PixelForwardLayer: disk not crossed by track" ;
    return SubTurbineCrossings();
  }

  HelixPlaneCrossing::PositionType  turbinePoint( turbineCrossing.position(thePath.second));
  HelixPlaneCrossing::DirectionType turbineDir( turbineCrossing.direction(thePath.second));
  int closestIndex = theBinFinder.binIndex(turbinePoint.phi());

  const BoundPlane& closestPlane( dynamic_cast<const BoundPlane&>( 
    theComps[closestIndex]->surface()));


  HelixArbitraryPlaneCrossing2Order theBladeCrossing(turbinePoint, turbineDir, rho);

  pair<bool,double> theClosestBladePath = theBladeCrossing.pathLength( closestPlane );
  LocalPoint closestPos = closestPlane.toLocal(GlobalPoint(theBladeCrossing.position(theClosestBladePath.second)) );
    
  float closestDist = closestPos.x(); // use fact that local X perp to global Y

  //int next = turbinePoint.phi() - closestPlane.position().phi() > 0 ? closest+1 : closest-1;
  int nextIndex = PhiLess()( closestPlane.position().phi(), turbinePoint.phi()) ? 
    closestIndex+1 : closestIndex-1;

  const BoundPlane& nextPlane( dynamic_cast<const BoundPlane&>( 
    theComps[ theBinFinder.binIndex(nextIndex)]->surface()));

  pair<bool,double> theNextBladePath    = theBladeCrossing.pathLength( nextPlane );
  LocalPoint nextPos = nextPlane.toLocal(GlobalPoint(theBladeCrossing.position(theNextBladePath.second)) );

  float nextDist = nextPos.x();

  if (fabs(closestDist) < fabs(nextDist)) {
    return SubTurbineCrossings( closestIndex, nextIndex, nextDist);
  }
  else {
    return SubTurbineCrossings( nextIndex, closestIndex, closestDist);
  }
}

float 
PixelForwardLayer::computeWindowSize( const GeomDet* det, 
				      const TrajectoryStateOnSurface& tsos, 
				      const MeasurementEstimator& est) const
{
  return est.maximalLocalDisplacement(tsos, det->surface()).x();
}


