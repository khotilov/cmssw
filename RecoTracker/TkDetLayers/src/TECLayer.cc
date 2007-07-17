#include "RecoTracker/TkDetLayers/interface/TECLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"

#include "TrackingTools/DetLayers/interface/DetLayerException.h"
#include "TrackingTools/DetLayers/interface/rangesIntersect.h"
#include "TrackingTools/GeomPropagators/interface/HelixForwardPlaneCrossing.h"
#include "TrackingTools/DetLayers/interface/PhiLess.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "RecoTracker/TkDetLayers/interface/TkDetUtil.h"


using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;

TECLayer::TECLayer(vector<const TECPetal*>& innerPetals,
		   vector<const TECPetal*>& outerPetals) : 
  theFrontComps(innerPetals.begin(),innerPetals.end()), 
  theBackComps(outerPetals.begin(),outerPetals.end())
{
  theComps.assign(theFrontComps.begin(),theFrontComps.end());
  theComps.insert(theComps.end(),theBackComps.begin(),theBackComps.end());

  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin();
      it!=theComps.end();it++){  
    theBasicComps.insert(theBasicComps.end(),	
			 (**it).basicComponents().begin(),
			 (**it).basicComponents().end());
  }


  //This should be no necessary. TO BE CHECKED
  //sort(theFrontPetals.begin(), theFrontPetals.end(), PetalLessPhi());
  //sort(theBackPetals.begin(), theBackPetals.end(), PetalLessPhi());

  // building disk for front and back petals
  setSurface( computeDisk( theComps ) );
  theFrontDisk = computeDisk( theFrontComps);
  theBackDisk  = computeDisk( theBackComps);

  // set up the bin finders
  theFrontBinFinder = BinFinderPhi(theFrontComps.front()->position().phi(),
  				   theFrontComps.size());
  theBackBinFinder  = BinFinderPhi(theBackComps.front()->position().phi(),
				   theBackComps.size());  

  //--------- DEBUG INFO --------------
  LogDebug("TkDetLayers") << "DEBUG INFO for TECLayer" << "\n"
			  << "TECLayer z,perp, innerRadius, outerR: " 
			  << this->position().z()    << " , "
			  << this->position().perp() << " , "
			  << this->specificSurface().innerRadius() << " , "
			  << this->specificSurface().outerRadius() ;
  

  for(vector<const GeometricSearchDet*>::const_iterator it=theFrontComps.begin(); 
      it!=theFrontComps.end(); it++){
    LogDebug("TkDetLayers") << "frontPetal phi,z,r: " 
	 << (*it)->surface().position().phi() << " , "
	 << (*it)->surface().position().z() <<   " , "
	 << (*it)->surface().position().perp() ;
  }

  for(vector<const GeometricSearchDet*>::const_iterator it=theBackComps.begin(); 
      it!=theBackComps.end(); it++){
    LogDebug("TkDetLayers") << "backPetal phi,z,r: " 
	 << (*it)->surface().position().phi() << " , "
	 << (*it)->surface().position().z() <<   " , "
	 << (*it)->surface().position().perp() ;
  }
  //----------------------------------- 


}



TECLayer::~TECLayer(){
  vector<const GeometricSearchDet*>::const_iterator i;
  for (i=theComps.begin(); i!=theComps.end(); i++) {
    delete *i;
  }
} 
  

void
TECLayer::groupedCompatibleDetsV( const TrajectoryStateOnSurface& tsos,
					  const Propagator& prop,
					   const MeasurementEstimator& est,
					   std::vector<DetGroup> & result) const {  
  SubLayerCrossings  crossings; 
  crossings = computeCrossings( tsos, prop.propagationDirection());
  if(! crossings.isValid()) return;

  vector<DetGroup> closestResult;
  addClosest( tsos, prop, est, crossings.closest(), closestResult); 
  LogDebug("TkDetLayers") << "in TECLayer, closestResult.size(): " << closestResult.size();

  // this differs from other groupedCompatibleDets logic, which DON'T check next in such cases!!!
  if(closestResult.empty()){
    vector<DetGroup> nextResult;
    addClosest( tsos, prop, est, crossings.other(), nextResult);   
    LogDebug("TkDetLayers") << "in TECLayer, nextResult.size(): " << nextResult.size();
    if(nextResult.empty())       return;
    

    DetGroupElement nextGel( nextResult.front().front());  
    int crossingSide = LayerCrossingSide().endcapSide( nextGel.trajectoryState(), prop);
    DetGroupMerger::orderAndMergeTwoLevels( closestResult, nextResult, result, 
					    crossings.closestIndex(), crossingSide);   
  }  
  else {
    DetGroupElement closestGel( closestResult.front().front());  
    float phiWindow = tkDetUtil::computeWindowSize( closestGel.det(), closestGel.trajectoryState(), est); 
    searchNeighbors( tsos, prop, est, crossings.closest(), phiWindow,
		     closestResult, false); 
    vector<DetGroup> nextResult;  
    searchNeighbors( tsos, prop, est, crossings.other(), phiWindow,
		     nextResult, true); 
    
    int crossingSide = LayerCrossingSide().endcapSide( closestGel.trajectoryState(), prop);
    DetGroupMerger::orderAndMergeTwoLevels( closestResult, nextResult, result,
					    crossings.closestIndex(), crossingSide);
  }
}

SubLayerCrossings TECLayer::computeCrossings(const TrajectoryStateOnSurface& startingState,
					     PropagationDirection propDir) const
{
  double rho( startingState.transverseCurvature());
  
  HelixPlaneCrossing::PositionType startPos( startingState.globalPosition() );
  HelixPlaneCrossing::DirectionType startDir( startingState.globalMomentum() );
  HelixForwardPlaneCrossing crossing(startPos,startDir,rho,propDir);

  pair<bool,double> frontPath = crossing.pathLength( *theFrontDisk);
  if (!frontPath.first) SubLayerCrossings();


  GlobalPoint gFrontPoint(crossing.position(frontPath.second));
  LogDebug("TkDetLayers") 
    << "in TECLayer,front crossing point: r,z,phi: (" 
    << gFrontPoint.perp() << ","
    << gFrontPoint.z() << "," 
    << gFrontPoint.phi() << ")" << endl;
  

  int frontIndex = theFrontBinFinder.binIndex(gFrontPoint.barePhi()); 
  SubLayerCrossing frontSLC( 0, frontIndex, gFrontPoint);



  pair<bool,double> backPath = crossing.pathLength( *theBackDisk);
  if (!backPath.first) SubLayerCrossings();


  GlobalPoint gBackPoint( crossing.position(backPath.second));
  LogDebug("TkDetLayers") 
    << "in TECLayer,back crossing point: r,z,phi: (" 
    << gBackPoint.perp() << "," 
    << gFrontPoint.z() << "," 
    << gBackPoint.phi() << ")" << endl;


  int backIndex = theBackBinFinder.binIndex(gBackPoint.barePhi());
  SubLayerCrossing backSLC( 1, backIndex, gBackPoint);

  
  // 0ss: frontDisk has index=0, backDisk has index=1
  float frontDist = std::abs(Geom::deltaPhi( double(gFrontPoint.barePhi()), 
					     double(theFrontDets[frontIndex]->surface().phi())));
  float backDist = std::abs(Geom::deltaPhi( double(gBackPoint.barePhi()), 
					    double(theBackDets[backIndex]->surface().phi())));
  

  if (frontDist < backDist) {
    return SubLayerCrossings( frontSLC, backSLC, 0);
  }
  else {
    return SubLayerCrossings( backSLC, frontSLC, 1);
  } 
}

bool TECLayer::addClosest( const TrajectoryStateOnSurface& tsos,
			   const Propagator& prop,
			   const MeasurementEstimator& est,
			   const SubLayerCrossing& crossing,
			   vector<DetGroup>& result) const
{
  const vector<const GeometricSearchDet*>& sub( subLayer( crossing.subLayerIndex()));
  const GeometricSearchDet* det(sub[crossing.closestDetIndex()]);

  LogDebug("TkDetLayers")  
    << "in TECLayer, adding petal at r,z,phi: (" 
    << det->position().perp() << "," 
    << det->position().z() << "," 
    << det->position().phi() << ")" << endl;

  return CompatibleDetToGroupAdder().add( *det, tsos, prop, est, result); 
}

void TECLayer::searchNeighbors( const TrajectoryStateOnSurface& tsos,
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

  const BinFinderPhi& binFinder = (crossing.subLayerIndex()==0 ? theFrontBinFinder : theBackBinFinder);

  typedef CompatibleDetToGroupAdder Adder;
  int half = sLayer.size()/2;  // to check if dets are called twice....
  for (int idet=negStartIndex; idet >= negStartIndex - half; idet--) {
    const GeometricSearchDet & neighborPetal = *sLayer[binFinder.binIndex(idet)];
    if (!overlap( gCrossingPos, neighborPetal, window)) break;
    if (!Adder::add( neighborPetal, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
  for (int idet=posStartIndex; idet < posStartIndex + half; idet++) {
    const GeometricSearchDet & neighborPetal = *sLayer[binFinder.binIndex(idet)];
    if (!overlap( gCrossingPos, neighborPetal, window)) break;
    if (!Adder::add( neighborPetal, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
}

bool TECLayer::overlap( const GlobalPoint& gpos, const GeometricSearchDet& gsdet, float phiWin) const
{
  float phi = gpos.barePhi();
  const TECPetal& petal = dynamic_cast<const TECPetal&>(gsdet);
  pair<float,float> phiRange(phi-phiWin,phi+phiWin);
  pair<float,float> petalPhiRange(petal.surface().phi() - 0.5*petal.specificSurface().phiExtension(),
				  petal.surface().phi() + 0.5*petal.specificSurface().phiExtension());


  return rangesIntersect(phiRange, petalPhiRange, PhiLess());
}



BoundDisk*
TECLayer::computeDisk( vector<const GeometricSearchDet*>& petals) const
{
  // Attention: it is assumed that the petals do belong to one layer, and are all
  // of the same rmin/rmax extension !!  
  
  const TECPetal* frontPetal = dynamic_cast<const TECPetal*>(petals.front());

  float rmin = frontPetal->specificSurface().innerRadius();
  float rmax = frontPetal->specificSurface().outerRadius();
  
  float theZmax(petals.front()->position().z());
  float theZmin(theZmax);
  for ( vector<const GeometricSearchDet*>::const_iterator i = petals.begin(); i != petals.end(); i++ ) {
    float zmin = (**i).position().z() - (**i).surface().bounds().thickness()/2.;
    float zmax = (**i).position().z() + (**i).surface().bounds().thickness()/2.;
    theZmax = max( theZmax, zmax);
    theZmin = min( theZmin, zmin);
  }

  float zPos = (theZmax+theZmin)/2.;
  PositionType pos(0.,0.,zPos);
  RotationType rot;

  return new BoundDisk( pos, rot,SimpleDiskBounds(rmin, rmax,    
						  theZmin-zPos, theZmax-zPos));
}

