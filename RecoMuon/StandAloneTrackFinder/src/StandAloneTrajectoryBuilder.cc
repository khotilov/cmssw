/** \class StandAloneTrajectoryBuilder
 *  Concrete class for the STA Muon reco 
 *
 *  $Date: 2007/02/16 13:31:23 $
 *  $Revision: 1.38 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 *  \author Stefano Lacaprara - INFN Legnaro
 */

#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneTrajectoryBuilder.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonRefitter.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonBackwardFilter.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"

#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

using namespace edm;
using namespace std;

StandAloneMuonTrajectoryBuilder::StandAloneMuonTrajectoryBuilder(const ParameterSet& par, 
								 const MuonServiceProxy* service):theService(service){
  const std::string metname = "Muon|RecoMuon|StandAloneMuonTrajectoryBuilder";
  
  LogTrace(metname) << "constructor called" << endl;

  // The navigation type:
  // "Direct","Standard"
  theNavigationType = par.getParameter<string>("NavigationType");
  
  // The inward-outward fitter (starts from seed state)
  ParameterSet refitterPSet = par.getParameter<ParameterSet>("RefitterParameters");
  refitterPSet.addParameter<string>("NavigationType",theNavigationType);
  theRefitter = new StandAloneMuonRefitter(refitterPSet,theService);

  // Fit direction
  string seedPosition = par.getParameter<string>("SeedPosition");
  
  if (seedPosition == "in" ) theSeedPosition = recoMuon::in;
  else if (seedPosition == "out" ) theSeedPosition = recoMuon::out;
  else 
    throw cms::Exception("StandAloneMuonRefitter constructor") 
      <<"Wrong seed position chosen in StandAloneMuonRefitter::StandAloneMuonRefitter ParameterSet"
      << "\n"
      << "Possible choices are:"
      << "\n"
      << "SeedPosition = in or SeedPosition = out";
  
  // Propagator for the seed extrapolation
  theSeedPropagatorName = par.getParameter<string>("SeedPropagator");
  
  // Disable/Enable the backward filter
  doBackwardRefit = par.getParameter<bool>("DoBackwardRefit");
  
  // Disable/Enable the smoothing of the trajectory
  doSmoothing = par.getParameter<bool>("DoSmoothing");
   
  if(doBackwardRefit){
    // The outward-inward fitter (starts from theRefitter outermost state)
    ParameterSet bwFilterPSet = par.getParameter<ParameterSet>("BWFilterParameters");
    //  theBWFilter = new StandAloneMuonBackwardFilter(bwFilterPSet,theService); // FIXME
    bwFilterPSet.addParameter<string>("NavigationType",theNavigationType);
    theBWFilter = new StandAloneMuonRefitter(bwFilterPSet,theService);

    theBWSeedType = bwFilterPSet.getParameter<string>("BWSeedType");
  }

  if(doSmoothing){
    // The outward-inward fitter (starts from theBWFilter innermost state)
    ParameterSet smootherPSet = par.getParameter<ParameterSet>("SmootherParameters");
    theSmoother = new StandAloneMuonSmoother(smootherPSet,theService);
  }
} 

void StandAloneMuonTrajectoryBuilder::setEvent(const edm::Event& event){
  theRefitter->setEvent(event);
   if(doBackwardRefit) theBWFilter->setEvent(event);
}

StandAloneMuonTrajectoryBuilder::~StandAloneMuonTrajectoryBuilder(){

  LogTrace("Muon|RecoMuon|StandAloneMuonTrajectoryBuilder") 
    << "StandAloneMuonTrajectoryBuilder destructor called" << endl;
  
  if(theRefitter) delete theRefitter;
  if(doBackwardRefit && theBWFilter) delete theBWFilter;
  if(doSmoothing && theSmoother) delete theSmoother;
}


MuonTrajectoryBuilder::TrajectoryContainer 
StandAloneMuonTrajectoryBuilder::trajectories(const TrajectorySeed& seed){ 

  const std::string metname = "Muon|RecoMuon|StandAloneMuonTrajectoryBuilder";
  MuonPatternRecoDumper debug;

  // the trajectory container. In principle starting from one seed we can
  // obtain more than one trajectory. TODO: this feature is not yet implemented!
  TrajectoryContainer trajectoryContainer;

  PropagationDirection fwDirection = (theSeedPosition == recoMuon::in) ? alongMomentum : oppositeToMomentum;  
  Trajectory trajectoryFW(seed,fwDirection);

  DetLayerWithState inputFromSeed = propagateTheSeedTSOS(seed); // it returns DetLayer-TSOS pair
  
  // refine the FTS given by the seed

  // the trajectory is filled in the refitter::refit
  refitter()->refit(inputFromSeed.second,inputFromSeed.first,trajectoryFW);

  // Get the last TSOS
  TrajectoryStateOnSurface tsosAfterRefit = refitter()->lastUpdatedTSOS();

  LogTrace(metname) << "StandAloneMuonTrajectoryBuilder REFITTER OUTPUT " << endl ;
  LogTrace(metname) << debug.dumpTSOS(tsosAfterRefit);
  

  if( refitter()->layers().size() ) 
    LogTrace(metname) << debug.dumpLayer( refitter()->lastDetLayer());
  else return trajectoryContainer; 
  
  LogTrace(metname) << "Number of DT/CSC/RPC chamber used (fw): " 
       << refitter()->getDTChamberUsed() << "/"
       << refitter()->getCSCChamberUsed() << "/"
       << refitter()->getRPCChamberUsed() <<endl;
  LogTrace(metname) << "Momentum: " <<tsosAfterRefit.freeState()->momentum();
  

  if(!doBackwardRefit){
    LogTrace(metname) << "Only forward refit requested. Any backward refit will be performed!"<<endl;
    
    if (  refitter()->getTotalChamberUsed() >= 2 && 
	  ((refitter()->getDTChamberUsed() + refitter()->getCSCChamberUsed()) >0 ||
	   refitter()->onlyRPC()) ){

      // Smoothing
      if (doSmoothing && !trajectoryFW.empty()){
	pair<bool,Trajectory> smoothingResult = smoother()->smooth(trajectoryFW);
	if (smoothingResult.first){
	  trajectoryContainer.push_back(new Trajectory(smoothingResult.second));
	  LogTrace(metname) << "StandAloneMuonTrajectoryBuilder SMOOTHER OUTPUT " << endl ;
	  LogTrace(metname) << debug.dumpTSOS(smoothingResult.second.lastMeasurement().updatedState());
	}
	else
	  trajectoryContainer.push_back(new Trajectory(trajectoryFW));
      }
      else
	trajectoryContainer.push_back(new Trajectory(trajectoryFW));
      
      LogTrace(metname)<< "Trajectory saved" << endl;
    }
    else LogTrace(metname)<< "Trajectory NOT saved. No enough number of tracking chamber used!" << endl;
    
    return trajectoryContainer;
  }

  // FIXME put the possible choices: (factory???)
  // fw_low-granularity + bw_high-granularity
  // fw_high-granularity + smoother
  // fw_low-granularity + bw_high-granularity + smoother (not yet sure...)

  // BackwardFiltering

  TrajectorySeed seedForBW;

  if(theBWSeedType == "noSeed"){
    TrajectorySeed seedVoid;
    seedForBW = seedVoid;
  }
  else if (theBWSeedType == "fromFWFit"){
    
    TrajectoryStateTransform tsTransform;
    
    PTrajectoryStateOnDet *seedTSOS =
      tsTransform.persistentState( tsosAfterRefit, trajectoryFW.lastMeasurement().recHit()->geographicalId().rawId());
    
    edm::OwnVector<TrackingRecHit> recHitsContainer; // FIXME!!
    TrajectorySeed fwFit(*seedTSOS,recHitsContainer,alongMomentum);

    seedForBW = fwFit;
  }
  else if (theBWSeedType == "fromGenerator"){
    seedForBW = seed;
  }
  else
    LogWarning(metname) << "Wrong seed type for the backward filter!";

  PropagationDirection bwDirection = (theSeedPosition == recoMuon::in) ?  oppositeToMomentum : alongMomentum;
  Trajectory trajectoryBW(seedForBW,bwDirection);

  // FIXME! under check!
  bwfilter()->refit(tsosAfterRefit,refitter()->lastDetLayer(),trajectoryBW);

  // Get the last TSOS
  TrajectoryStateOnSurface tsosAfterBWRefit = bwfilter()->lastUpdatedTSOS();

  LogTrace(metname) << "StandAloneMuonTrajectoryBuilder BW FILTER OUTPUT " << endl ;
  LogTrace(metname) << debug.dumpTSOS(tsosAfterBWRefit);

  LogTrace(metname) 
    << "Number of RecHits: " << trajectoryBW.foundHits() << "\n"
    << "Number of DT/CSC/RPC chamber used (bw): " 
    << bwfilter()->getDTChamberUsed() << "/"
    << bwfilter()->getCSCChamberUsed() << "/" 
    << bwfilter()->getRPCChamberUsed();
  
  // The trajectory is good if there are at least 2 chamber used in total and at
  // least 1 "tracking" (DT or CSC)
  if (  bwfilter()->getTotalChamberUsed() >= 2 && 
	((bwfilter()->getDTChamberUsed() + bwfilter()->getCSCChamberUsed()) >0 ||
	 bwfilter()->onlyRPC()) ) {
    
    if (doSmoothing && !trajectoryBW.empty()){
      pair<bool,Trajectory> smoothingResult = smoother()->smooth(trajectoryBW);
      if (smoothingResult.first){
     	trajectoryContainer.push_back(new Trajectory(smoothingResult.second));
	LogTrace(metname) << "StandAloneMuonTrajectoryBuilder SMOOTHER OUTPUT " << endl ;
	LogTrace(metname) << debug.dumpTSOS(smoothingResult.second.lastMeasurement().updatedState());
      }
      else
	trajectoryContainer.push_back(new Trajectory(trajectoryBW));
    }
    else
      trajectoryContainer.push_back(new Trajectory(trajectoryBW));
    
    LogTrace(metname)<< "Trajectory saved" << endl;
    
  }
  //if the trajectory is not saved, but at least two chamber are used in the
  //forward filtering, try to build a new trajectory starting from the old
  //trajectory w/o the latest measurement and a looser chi2 cut
  else if ( refitter()->getTotalChamberUsed() >= 2 ) {
    LogTrace(metname)<< "Trajectory NOT saved. Second Attempt." << endl
		     << "FIRST MEASUREMENT KILLED" << endl; // FIXME: why???
    // FIXME:
    // a better choice could be: identify the poorest one, exclude it, redo
    // the fw and bw filtering. Or maybe redo only the bw without the excluded
    // measure. As first step I will port the ORCA algo, then I will move to the
    // above pattern.
    
  }
  else
    LogTrace(metname)<< "Trajectory NOT saved" << endl;
  return trajectoryContainer;
}


StandAloneMuonTrajectoryBuilder::DetLayerWithState
StandAloneMuonTrajectoryBuilder::propagateTheSeedTSOS(const TrajectorySeed& seed){

  const std::string metname = "Muon|RecoMuon|StandAloneMuonTrajectoryBuilder";
  MuonPatternRecoDumper debug;

  // Get the Trajectory State on Det (persistent version of a TSOS) from the seed
  PTrajectoryStateOnDet pTSOD = seed.startingState();
  
  // Transform it in a TrajectoryStateOnSurface
  LogTrace(metname)<<"Transform PTrajectoryStateOnDet in a TrajectoryStateOnSurface"<<endl;
  TrajectoryStateTransform tsTransform;

  DetId seedDetId(pTSOD.detId());

  const GeomDet* gdet = theService->trackingGeometry()->idToDet( seedDetId );

  TrajectoryStateOnSurface initialState = tsTransform.transientState(pTSOD, &(gdet->surface()), 
								     &*theService->magneticField());

  LogTrace(metname)<<"Seed's Pt: "<<initialState.freeState()->momentum().perp()<<endl;

  LogTrace(metname)<< "Seed's id: "<< endl ;
  LogTrace(metname) << debug.dumpMuonId(seedDetId);
  
  // Get the layer on which the seed relies
  const DetLayer *initialLayer = theService->detLayerGeometry()->idToLayer( seedDetId );

  LogTrace(metname)<< "Seed's detLayer: "<< endl ;
  LogTrace(metname) << debug.dumpLayer(initialLayer);

  LogTrace(metname)<< "TrajectoryStateOnSurface before propagation:" << endl;
  LogTrace(metname) << debug.dumpTSOS(initialState);


  PropagationDirection detLayerOrder = (theSeedPosition == recoMuon::in) ? oppositeToMomentum : alongMomentum;

  // ask for compatible layers
  vector<const DetLayer*> detLayers;

  if(theNavigationType == "Standard")
    detLayers = initialLayer->compatibleLayers( *initialState.freeState(),detLayerOrder); 
  else if (theNavigationType == "Direct"){
    DirectMuonNavigation navigation( &*theService->detLayerGeometry() );
    detLayers = navigation.compatibleLayers( *initialState.freeState(),detLayerOrder);
  }
  else
    edm::LogError(metname) << "No Properly Navigation Selected!!"<<endl;

 
  LogTrace(metname) << "There are "<< detLayers.size() <<" compatible layers"<<endl;
  
  DetLayerWithState result = DetLayerWithState(initialLayer,initialState);

  if(detLayers.size()){

    LogTrace(metname) << "Compatible layers:"<<endl;
    for( vector<const DetLayer*>::const_iterator layer = detLayers.begin(); 
	 layer != detLayers.end(); layer++){
      LogTrace(metname) << debug.dumpMuonId((*layer)->basicComponents().front()->geographicalId());
      LogTrace(metname) << debug.dumpLayer(*layer);
    }

    const DetLayer* finalLayer = detLayers.back();

    if(theSeedPosition == recoMuon::in) LogTrace(metname) << "Most internal one:"<<endl;
    else LogTrace(metname) << "Most external one:"<<endl;
    
    LogTrace(metname) << debug.dumpLayer(finalLayer);
    
    const TrajectoryStateOnSurface propagatedState = 
      theService->propagator(theSeedPropagatorName)->propagate(initialState,
							       finalLayer->surface());

    if(propagatedState.isValid()){
      result = DetLayerWithState(finalLayer,propagatedState);
      
      LogTrace(metname) << "Trajectory State on Surface after the extrapolation"<<endl;
      LogTrace(metname) << debug.dumpTSOS(propagatedState);
      LogTrace(metname) << debug.dumpLayer(finalLayer);
    }
    else 
      LogTrace(metname)<< "Extrapolation failed. Keep the original seed state" <<endl;
  }
  else
    LogTrace(metname)<< "No compatible layers. Keep the original seed state" <<endl;
  
  return result;
}
