/**  \class GlobalMuonProducer
 * 
 *   Global muon reconstructor:
 *   reconstructs muons using DT, CSC, RPC and tracker
 *   information,<BR>
 *   starting from a standalone reonstructed muon.
 *
 *   $Date: 2007/02/16 13:33:11 $
 *   $Revision: 1.23 $
 *
 *   \author  R.Bellan - INFN TO
 */

// Framework
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoMuon/GlobalMuonProducer/src/GlobalMuonProducer.h"

// TrackFinder and specific GLB Trajectory Builder
#include "RecoMuon/GlobalTrackFinder/interface/GlobalMuonTrajectoryBuilder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackFinder.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

// Input and output collection
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

using namespace edm;
using namespace std;

//
// constructor with config
//
GlobalMuonProducer::GlobalMuonProducer(const ParameterSet& parameterSet) {

  LogTrace("Muon|RecoMuon|GlobalMuonProducer") << "constructor called" << endl;

  // Parameter set for the Builder
  ParameterSet trajectoryBuilderParameters = parameterSet.getParameter<ParameterSet>("GLBTrajBuilderParameters");
  InputTag trackCollectionTag = parameterSet.getParameter<InputTag>("TrackerCollectionLabel");
  trajectoryBuilderParameters.addParameter<InputTag>("TrackerCollectionLabel",trackCollectionTag);

  // STA Muon Collection Label
  theSTACollectionLabel = parameterSet.getParameter<InputTag>("MuonCollectionLabel");

  // STA semi-persistent flag
  theSTATrajectoryFlag = parameterSet.getParameter<bool>("MuonTrajectoryAvailable");

  // service parameters
  ParameterSet serviceParameters = parameterSet.getParameter<ParameterSet>("ServiceParameters");

  // TrackLoader parameters
  ParameterSet trackLoaderParameters = parameterSet.getParameter<ParameterSet>("TrackLoaderParameters");
  
  // the services
  theService = new MuonServiceProxy(serviceParameters);
  
  // instantiate the concrete trajectory builder in the Track Finder
  MuonTrackLoader* mtl = new MuonTrackLoader(trackLoaderParameters,theService);
  GlobalMuonTrajectoryBuilder* gmtb = new GlobalMuonTrajectoryBuilder(trajectoryBuilderParameters, theService);

  theTrackFinder = new MuonTrackFinder(gmtb, mtl);
  
  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  produces<vector<Trajectory> >() ;

  produces<reco::MuonCollection>();

}


//
// destructor
//
GlobalMuonProducer::~GlobalMuonProducer() {

  LogTrace("Muon|RecoMuon|GlobalMuonProducer") << "destructor called" << endl;
  if (theService) delete theService;
  if (theTrackFinder) delete theTrackFinder;

}


//
// reconstruct muons
//
void GlobalMuonProducer::produce(Event& event, const EventSetup& eventSetup) {
  const string metname = "Muon|RecoMuon|GlobalMuonProducer";  
  LogTrace(metname)<<endl<<endl<<endl;
  LogTrace(metname)<<"Global Muon Reconstruction started"<<endl;  

  typedef vector<Trajectory> TrajColl;

  // Update the services
  theService->update(eventSetup);

  // Take the STA muon container(s)
  LogTrace(metname)<<"Taking the Stand Alone Muons "<<theSTACollectionLabel.label()<<endl;

  Handle<reco::TrackCollection> staMuons;
  event.getByLabel(theSTACollectionLabel,staMuons);

  Handle<vector<Trajectory> > staMuonsTraj;

  if(theSTATrajectoryFlag) {
    event.getByLabel(theSTACollectionLabel.label(),staMuonsTraj);      
    LogTrace(metname)<<"Track Reconstruction (tracks, trajs) "<< staMuons.product()->size() << " " << staMuonsTraj.product()->size() <<endl;
  } 

  theTrackFinder->reconstruct(staMuons, staMuonsTraj, event);      

  
  LogTrace(metname)<<"Event loaded"
                   <<"================================"
                   <<endl<<endl;
    
}
