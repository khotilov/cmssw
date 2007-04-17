#include "RecoTracker/TrackProducer/interface/GsfTrackProducer.h"
// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

GsfTrackProducer::GsfTrackProducer(const edm::ParameterSet& iConfig):
  GsfTrackProducerBase(iConfig.getParameter<bool>("TrajectoryInEvent")),
  theAlgo(iConfig)
{
  setConf(iConfig);
  setSrc( iConfig.getParameter<std::string>( "src" ));
  setAlias( iConfig.getParameter<std::string>( "@module_label" ) );
//   string a = alias_;
//   a.erase(a.size()-6,a.size());
  //register your products
  produces<reco::GsfTrackCollection>().setBranchAlias( alias_ + "GsfTracks" );
  produces<reco::TrackExtraCollection>().setBranchAlias( alias_ + "TrackExtras" );
  produces<reco::GsfTrackExtraCollection>().setBranchAlias( alias_ + "GsfTrackExtras" );
  produces<TrackingRecHitCollection>().setBranchAlias( alias_ + "RecHits" );
  produces<std::vector<Trajectory> >() ;

}


void GsfTrackProducer::produce(edm::Event& theEvent, const edm::EventSetup& setup)
{
  edm::LogInfo("GsfTrackProducer") << "Analyzing event number: " << theEvent.id() << "\n";
  //
  // create empty output collections
  //
  std::auto_ptr<TrackingRecHitCollection> outputRHColl (new TrackingRecHitCollection);
  std::auto_ptr<reco::GsfTrackCollection> outputTColl(new reco::GsfTrackCollection);
  std::auto_ptr<reco::TrackExtraCollection> outputTEColl(new reco::TrackExtraCollection);
  std::auto_ptr<reco::GsfTrackExtraCollection> outputGsfTEColl(new reco::GsfTrackExtraCollection);
  std::auto_ptr<std::vector<Trajectory> >    outputTrajectoryColl(new std::vector<Trajectory>);
  //
  //declare and get stuff to be retrieved from ES
  //
  edm::ESHandle<TrackerGeometry> theG;
  edm::ESHandle<MagneticField> theMF;
  edm::ESHandle<TrajectoryFitter> theFitter;
  edm::ESHandle<Propagator> thePropagator;
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  getFromES(setup,theG,theMF,theFitter,thePropagator,theBuilder);

  //
  //declare and get TrackColection to be retrieved from the event
  //
    AlgoProductCollection algoResults;
  try{  
    edm::Handle<TrackCandidateCollection> theTCCollection;
    getFromEvt(theEvent,theTCCollection);
    
    //
    //run the algorithm  
    //
    LogDebug("GsfTrackProducer") << "run the algorithm" << "\n";
    theAlgo.runWithCandidate(theG.product(), theMF.product(), *theTCCollection, 
			     theFitter.product(), thePropagator.product(), theBuilder.product(), algoResults);
  } catch (cms::Exception &e){ edm::LogInfo("GsfTrackProducer") << "cms::Exception caught!!!" << "\n" << e << "\n";}
  //
  //put everything in the event
  putInEvt(theEvent, outputRHColl, outputTColl, outputTEColl, outputGsfTEColl,
	   outputTrajectoryColl, algoResults);
  LogDebug("GsfTrackProducer") << "end" << "\n";
}


// std::vector<reco::TransientTrack> GsfTrackProducer::getTransient(edm::Event& theEvent, const edm::EventSetup& setup)
// {
//   edm::LogInfo("GsfTrackProducer") << "Analyzing event number: " << theEvent.id() << "\n";
//   //
//   // create empty output collections
//   //
//   std::vector<reco::TransientTrack> ttks;

//   //
//   //declare and get stuff to be retrieved from ES
//   //
//   edm::ESHandle<TrackerGeometry> theG;
//   edm::ESHandle<MagneticField> theMF;
//   edm::ESHandle<TrajectoryFitter> theFitter;
//   edm::ESHandle<Propagator> thePropagator;
//   edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
//   getFromES(setup,theG,theMF,theFitter,thePropagator,theBuilder);

//   //
//   //declare and get TrackColection to be retrieved from the event
//   //
//   AlgoProductCollection algoResults;
//   try{  
//     edm::Handle<TrackCandidateCollection> theTCCollection;
//     getFromEvt(theEvent,theTCCollection);
    
//     //
//     //run the algorithm  
//     //
//     LogDebug("GsfTrackProducer") << "run the algorithm" << "\n";
//     theAlgo.runWithCandidate(theG.product(), theMF.product(), *theTCCollection, 
// 			     theFitter.product(), thePropagator.product(), theBuilder.product(), algoResults);
//   } catch (cms::Exception &e){ edm::LogInfo("GsfTrackProducer") << "cms::Exception caught!!!" << "\n" << e << "\n";}


//   for (AlgoProductCollection::iterator prod=algoResults.begin();prod!=algoResults.end(); prod++){
//     ttks.push_back( reco::TransientTrack(*((*prod).second),thePropagator.product()->magneticField() ));
//   }

//   LogDebug("GsfTrackProducer") << "end" << "\n";

//   return ttks;
// }


