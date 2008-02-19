/** \class StandAloneMuonRefitter
 *  The inward-outward fitter (starts from seed state).
 *
 *  $Date: 2008/02/19 08:46:56 $
 *  $Revision: 1.39 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 *  \author S. Lacaprara - INFN Legnaro
 */
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonRefitter.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FIXME: remove this
#include "FWCore/Framework/interface/Event.h"

#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonTrajectoryUpdator.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <vector>

using namespace edm;
using namespace std;

StandAloneMuonRefitter::StandAloneMuonRefitter(const ParameterSet& par,
					       const MuonServiceProxy* service)
:theService(service),
 theOverlappingChambersFlag(true)
{
  // Fit direction
  string fitDirectionName = par.getParameter<string>("FitDirection");

  if (fitDirectionName == "insideOut" ) theFitDirection = insideOut;
  else if (fitDirectionName == "outsideIn" ) theFitDirection = outsideIn;
  else 
    throw cms::Exception("StandAloneMuonRefitter constructor") 
      <<"Wrong fit direction chosen in StandAloneMuonRefitter::StandAloneMuonRefitter ParameterSet"
      << "\n"
      << "Possible choices are:"
      << "\n"
      << "FitDirection = insideOut or FitDirection = outsideIn";
  
  // The max allowed chi2 to accept a rechit in the fit
  theMaxChi2 = par.getParameter<double>("MaxChi2");

  // The errors of the trajectory state are multiplied by nSigma 
  // to define acceptance of BoundPlane and maximalLocalDisplacement
  theNSigma = par.getParameter<double>("NumberOfSigma"); // default = 3.

  // The navigation type:
  // "Direct","Standard"
  theNavigationType = par.getParameter<string>("NavigationType");
  
  // The estimator: makes the decision wheter a measure is good or not
  // it isn't used by the updator which does the real fit. In fact, in principle,
  // a looser request onto the measure set can be requested 
  // (w.r.t. the request on the accept/reject measure in the fit)
  theEstimator = new Chi2MeasurementEstimator(theMaxChi2,theNSigma);
  
  thePropagatorName = par.getParameter<string>("Propagator");

  theBestMeasurementFinder = new MuonBestMeasurementFinder();

  // Muon trajectory updator parameters
  ParameterSet muonUpdatorPSet = par.getParameter<ParameterSet>("MuonTrajectoryUpdatorParameters");
  
  // the updator needs the fit direction
  theMuonUpdator = new MuonTrajectoryUpdator(muonUpdatorPSet,
					     fitDirection() );

  // Measurement Extractor: enable the measure for each muon sub detector
  bool enableDTMeasurement = par.getParameter<bool>("EnableDTMeasurement");
  bool enableCSCMeasurement = par.getParameter<bool>("EnableCSCMeasurement");
  bool enableRPCMeasurement = par.getParameter<bool>("EnableRPCMeasurement");

  theMeasurementExtractor = new MuonDetLayerMeasurements(par.getParameter<InputTag>("DTRecSegmentLabel"),
							 par.getParameter<InputTag>("CSCRecSegmentLabel"),
							 par.getParameter<InputTag>("RPCRecSegmentLabel"),
							 enableDTMeasurement,
							 enableCSCMeasurement,
							 enableRPCMeasurement);
  
  theRPCLoneliness = (!(enableDTMeasurement && enableCSCMeasurement)) ? enableRPCMeasurement : false;
}

StandAloneMuonRefitter::~StandAloneMuonRefitter(){

  LogTrace("Muon|RecoMuon|StandAloneMuonRefitter")
    <<"StandAloneMuonRefitter destructor called"<<endl;
  
  delete theEstimator;
  delete theMuonUpdator;
  delete theMeasurementExtractor;
  delete theBestMeasurementFinder;
}

const Propagator* StandAloneMuonRefitter::propagator() const { 
  return &*theService->propagator(thePropagatorName); 
}

/// Return the propagation direction
PropagationDirection StandAloneMuonRefitter::propagationDirection() const{
  if( fitDirection() == 0 ) return alongMomentum;
  else if ( fitDirection() == 1 ) return oppositeToMomentum;
  else return anyDirection;
}


void StandAloneMuonRefitter::reset(){
  totalChambers = dtChambers = cscChambers = rpcChambers = 0;
  
  theLastUpdatedTSOS =  theLastButOneUpdatedTSOS = TrajectoryStateOnSurface();

  theMuonUpdator->makeFirstTime();

  theDetLayers.clear();
}

void StandAloneMuonRefitter::setEvent(const Event& event){
  theMeasurementExtractor->setEvent(event);
}


void StandAloneMuonRefitter::incrementChamberCounters(const DetLayer *layer){

  if(layer->subDetector()==GeomDetEnumerators::DT) dtChambers++; 
  else if(layer->subDetector()==GeomDetEnumerators::CSC) cscChambers++; 
  else if(layer->subDetector()==GeomDetEnumerators::RPCBarrel || layer->subDetector()==GeomDetEnumerators::RPCEndcap) rpcChambers++; 
  else 
    LogError("Muon|RecoMuon|StandAloneMuonRefitter")
      << "Unrecognized module type in incrementChamberCounters";
  // FIXME:
  //   << layer->module() << " " <<layer->Part() << endl;
  
  totalChambers++;
}


vector<const DetLayer*> StandAloneMuonRefitter::compatibleLayers(const DetLayer *initialLayer,
								 FreeTrajectoryState& fts,
								 PropagationDirection propDir){
  vector<const DetLayer*> detLayers;

  if(theNavigationType == "Standard"){
    // ask for compatible layers
    detLayers = initialLayer->compatibleLayers(fts,propDir);  
    // I have to fit by hand the first layer until the seedTSOS is defined on the first rechit layer
    // In fact the first layer is not returned by initialLayer->compatibleLayers.
    detLayers.insert(detLayers.begin(),initialLayer);
  }
  else if (theNavigationType == "Direct"){
    DirectMuonNavigation navigation(&*theService->detLayerGeometry());
    detLayers = navigation.compatibleLayers(fts,propDir);
  }
  else
    edm::LogError("Muon|RecoMuon|StandAloneMuonRefitter") << "No Properly Navigation Selected!!"<<endl;
  
  return detLayers;
}


void StandAloneMuonRefitter::refit(const TrajectoryStateOnSurface& initialTSOS,
				   const DetLayer* initialLayer, Trajectory &trajectory){
  
  const std::string metname = "Muon|RecoMuon|StandAloneMuonRefitter";

  // reset the refitter each seed refinement
  reset();
  
  MuonPatternRecoDumper debug;
  
  LogTrace(metname) << "Starting the refit"<<endl; 

  // this is the most outward TSOS (updated or predicted) onto a DetLayer
  TrajectoryStateOnSurface lastTSOS = theLastUpdatedTSOS = theLastButOneUpdatedTSOS = initialTSOS;
  
  double eta0 = initialTSOS.freeTrajectoryState()->momentum().eta();
  vector<const DetLayer*> detLayers = compatibleLayers(initialLayer,*initialTSOS.freeTrajectoryState(),
						       propagationDirection());  

  const DetLayer* lastLayer=initialLayer;
  
  while(!detLayers.empty()){
    LogTrace(metname)<<"compatible layers found: "<<detLayers.size()<<endl;

    // the layers are ordered in agreement with the fit/propagation direction 
    for (vector<const DetLayer*>::const_iterator layer = detLayers.begin(); 
	 layer!= detLayers.end(); ++layer ) {
      
      //    bool firstTime = true;
      LogTrace(metname) << debug.dumpLayer(*layer);
      LogTrace(metname) << "search Trajectory Measurement from: " << lastTSOS.globalPosition();
      
      // pick the best measurement from each group
      std::vector<TrajectoryMeasurement> bestMeasurements = findBestMeasurements(*layer, lastTSOS);
      
      // RB: Different ways can be choosen if no bestMeasurement is available:
      // 1- check on lastTSOS-initialTSOS eta difference
      // 2- check on lastTSOS-lastButOneUpdatedTSOS eta difference
      // After this choice:
      // A- extract the measurements compatible with the initialTSOS (seed)
      // B- extract the measurements compatible with the lastButOneUpdatedTSOS
      // In ORCA the choice was 1A. Here I will try 1B and if it fail I'll try 1A
      // another possibility could be 2B and then 1A.

      // if no measurement found and the current TSOS has an eta very different
      // wrt the initial one (i.e. seed), then try to find the measurements
      // according to the lastButOne FTS. (1B)
      double lastdEta = fabs(lastTSOS.freeTrajectoryState()->momentum().eta() - eta0);
      if( bestMeasurements.empty() && lastdEta > 0.1) {
	LogTrace(metname) << "No measurement and big eta variation wrt seed" << endl
			  << "trying with lastButOneUpdatedTSOS";
	bestMeasurements = findBestMeasurements(*layer, theLastButOneUpdatedTSOS);
      }
    
      //if no measurement found and the current FTS has an eta very different
      //wrt the initial one (i.e. seed), then try to find the measurements
      //according to the initial FTS. (1A)
      if( bestMeasurements.empty() && lastdEta > 0.1) {
	LogTrace(metname) << "No measurement and big eta variation wrt seed" << endl
			  << "tryng with seed TSOS";
	bestMeasurements = findBestMeasurements(*layer, initialTSOS);
      }
    
      // FIXME: uncomment this line!!
      // if(!bestMeasurement && firstTime) break;

      if(!bestMeasurements.empty()) {
	bool added = false;
	for(std::vector<TrajectoryMeasurement>::const_iterator tmItr = bestMeasurements.begin();
	    tmItr != bestMeasurements.end(); ++tmItr){
	  added |= update(*layer, &(*tmItr), trajectory);
	}
	if(added) {
	  lastTSOS = theLastUpdatedTSOS;
	  incrementChamberCounters(*layer);
	  theDetLayers.push_back(*layer);
	  lastLayer = (*layer);
	  //   detLayers = (*layer)->compatibleLayers(*theLastUpdatedTSOS.freeTrajectoryState(),
	  // 						 propagationDirection());  
	  break;
	}
      }
      else LogTrace(metname)<<"No best measurement found"<<endl;
      lastLayer = *layer;
    } // loop over layers
    detLayers = lastLayer->compatibleLayers(*theLastUpdatedTSOS.freeTrajectoryState(),
					    propagationDirection());  
    
  }
}

std::vector<TrajectoryMeasurement>
StandAloneMuonRefitter::findBestMeasurements(const DetLayer* layer,
                                             const TrajectoryStateOnSurface& tsos){

  const std::string metname = "Muon|RecoMuon|StandAloneMuonRefitter";

  std::vector<TrajectoryMeasurement> result;
  std::vector<TrajectoryMeasurement> measurements;

  if(theOverlappingChambersFlag && layer->hasGroups()){
    
    std::vector<TrajectoryMeasurementGroup> measurementGroups =
      theMeasurementExtractor->groupedMeasurements(layer, tsos, *propagator(), *estimator());

    for(std::vector<TrajectoryMeasurementGroup>::const_iterator tmGroupItr = measurementGroups.begin();
        tmGroupItr != measurementGroups.end(); ++tmGroupItr){
    
      measurements = tmGroupItr->measurements();
      LogTrace(metname) << "Number of Trajectory Measurement: " << measurements.size();
      
      const TrajectoryMeasurement* bestMeasurement 
	= bestMeasurementFinder()->findBestMeasurement(measurements,  propagator());
      
      if(bestMeasurement) result.push_back(*bestMeasurement);
    }
  } 
  else{
    measurements = theMeasurementExtractor->measurements(layer, tsos, *propagator(), *estimator());
    LogTrace(metname) << "Number of Trajectory Measurement: " << measurements.size();
    const TrajectoryMeasurement* bestMeasurement 
      = bestMeasurementFinder()->findBestMeasurement(measurements,  
						     propagator());
    if(bestMeasurement) result.push_back(*bestMeasurement);
  }
  return result;
}




bool StandAloneMuonRefitter::update(const DetLayer * layer, 
                                    const TrajectoryMeasurement * meas, 
                                    Trajectory & trajectory)
{
  const std::string metname = "Muon|RecoMuon|StandAloneMuonRefitter";
  MuonPatternRecoDumper debug;

  LogTrace(metname)<<"best measurement found" << "\n"
                   <<"updating the trajectory..."<<endl;
  pair<bool,TrajectoryStateOnSurface> result = updator()->update(meas,
                                                                 trajectory,
                                                                 propagator());
  LogTrace(metname)<<"trajectory updated: "<<result.first<<endl;
  LogTrace(metname) << debug.dumpTSOS(result.second);

  if(result.first){
    theLastButOneUpdatedTSOS = theLastUpdatedTSOS;
    theLastUpdatedTSOS = result.second;
  }
  return result.first;
}
