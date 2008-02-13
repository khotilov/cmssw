#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
//
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/EgammaTrackReco/interface/TrackCandidateSuperClusterAssociation.h"
//
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//  Abstract classes for the conversion tracking components
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionSeedFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionTrackFinder.h"
// Class header file
#include "RecoEgamma/EgammaPhotonProducers/interface/ConversionTrackCandidateProducer.h"
//
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/NavigationSchoolRecord.h"

#include "RecoEgamma/EgammaPhotonAlgos/interface/OutInConversionSeedFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/InOutConversionSeedFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/OutInConversionTrackFinder.h"
#include "RecoEgamma/EgammaPhotonAlgos/interface/InOutConversionTrackFinder.h"

ConversionTrackCandidateProducer::ConversionTrackCandidateProducer(const edm::ParameterSet& config) : 
  conf_(config), 
  theNavigationSchool_(0), 
  theOutInSeedFinder_(0), 
  theOutInTrackFinder_(0), 
  theInOutSeedFinder_(0),
  theInOutTrackFinder_(0),
  isInitialized(0)

{


  //  LogDebug("ConversionTrackCandidateProducer") << "ConversionTrackCandidateProducer CTOR " << "\n";
 LogDebug("ConversionTrackCandidateProducer") << "ConversionTrackCandidateProducer CTOR " << "\n";
  
   
  // use onfiguration file to setup input/output collection names
 
  bcProducer_             = conf_.getParameter<std::string>("bcProducer");
  bcBarrelCollection_     = conf_.getParameter<std::string>("bcBarrelCollection");
  bcEndcapCollection_     = conf_.getParameter<std::string>("bcEndcapCollection");
  
  scHybridBarrelProducer_       = conf_.getParameter<std::string>("scHybridBarrelProducer");
  scIslandEndcapProducer_       = conf_.getParameter<std::string>("scIslandEndcapProducer");
  
  scHybridBarrelCollection_     = conf_.getParameter<std::string>("scHybridBarrelCollection");
  scIslandEndcapCollection_     = conf_.getParameter<std::string>("scIslandEndcapCollection");
  
  OutInTrackCandidateCollection_ = conf_.getParameter<std::string>("outInTrackCandidateCollection");
  InOutTrackCandidateCollection_ = conf_.getParameter<std::string>("inOutTrackCandidateCollection");


  OutInTrackSuperClusterAssociationCollection_ = conf_.getParameter<std::string>("outInTrackCandidateSCAssociationCollection");
  InOutTrackSuperClusterAssociationCollection_ = conf_.getParameter<std::string>("inOutTrackCandidateSCAssociationCollection");


  // Register the product
  produces< TrackCandidateCollection > (OutInTrackCandidateCollection_);
  produces< TrackCandidateCollection > (InOutTrackCandidateCollection_);
  //
  produces< reco::TrackCandidateSuperClusterAssociationCollection > ( OutInTrackSuperClusterAssociationCollection_);
  produces< reco::TrackCandidateSuperClusterAssociationCollection > ( InOutTrackSuperClusterAssociationCollection_);
  

}

ConversionTrackCandidateProducer::~ConversionTrackCandidateProducer() {


  delete theOutInSeedFinder_; 
  delete theOutInTrackFinder_;
  delete theInOutSeedFinder_;  
  delete theInOutTrackFinder_;


}


void  ConversionTrackCandidateProducer::beginJob (edm::EventSetup const & theEventSetup) {
  nEvt_=0;
  //get magnetic field
  edm::LogInfo("ConversionTrackCandidateProducer") << " get magnetic field" << "\n";
  theEventSetup.get<IdealMagneticFieldRecord>().get(theMF_);  


  theEventSetup .get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker_ );


  // get the measurement tracker   
  edm::ESHandle<MeasurementTracker> measurementTrackerHandle;
  theEventSetup.get<CkfComponentsRecord>().get(measurementTrackerHandle);
  theMeasurementTracker_ = measurementTrackerHandle.product();
  
  
  edm::ESHandle<NavigationSchool> nav;
  theEventSetup.get<NavigationSchoolRecord>().get("SimpleNavigationSchool", nav);
  theNavigationSchool_ = nav.product();
  


  
  // get the Out In Seed Finder  
  edm::LogInfo("ConversionTrackCandidateProducer") << " get the OutInSeedFinder" << "\n";
  theOutInSeedFinder_ = new OutInConversionSeedFinder (   &(*theMF_) ,  theMeasurementTracker_ );
  
  // get the Out In Track Finder
  edm::LogInfo("ConversionTrackCandidateProducer") << " get the OutInTrackFinder" << "\n";
  theOutInTrackFinder_ = new OutInConversionTrackFinder ( theEventSetup, conf_, &(*theMF_),  theMeasurementTracker_  );
  
  
  // get the In Out Seed Finder  
  edm::LogInfo("ConversionTrackCandidateProducer") << " get the InOutSeedFinder" << "\n";
  theInOutSeedFinder_ = new InOutConversionSeedFinder (  &(*theMF_) ,  theMeasurementTracker_ );
  
  
  
  // get the In Out Track Finder
  edm::LogInfo("ConversionTrackCandidateProducer") << " get the InOutTrackFinder" << "\n";
  theInOutTrackFinder_ = new InOutConversionTrackFinder ( theEventSetup, conf_, &(*theMF_),  theMeasurementTracker_  );
  
  
}



void ConversionTrackCandidateProducer::produce(edm::Event& theEvent, const edm::EventSetup& theEventSetup) {
  
  using namespace edm;
  nEvt_++;
  edm::LogInfo("ConversionTrackCandidateProducer") << "ConversionTrackCandidateProducer Analyzing event number: " << theEvent.id() << " Global Counter " << nEvt_ << "\n";
  LogDebug("ConversionTrackCandidateProducer") << "ConversionTrackCandidateProducer Analyzing event number " <<   theEvent.id() <<  " Global Counter " << nEvt_ << "\n";
  
  
  // Update MeasurementTracker
  theMeasurementTracker_->update(theEvent);
  // Set the navigation school  
  NavigationSetter setter(*theNavigationSchool_);  
  
  //
  // create empty output collections
  //
  //  Out In Track Candidates
  std::auto_ptr<TrackCandidateCollection> outInTrackCandidate_p(new TrackCandidateCollection); 
  //  In Out  Track Candidates
  std::auto_ptr<TrackCandidateCollection> inOutTrackCandidate_p(new TrackCandidateCollection); 
  //   Track Candidate  Super Cluster Association
  std::auto_ptr<reco::TrackCandidateSuperClusterAssociationCollection> outInAssoc_p(new reco::TrackCandidateSuperClusterAssociationCollection);
  std::auto_ptr<reco::TrackCandidateSuperClusterAssociationCollection> inOutAssoc_p(new reco::TrackCandidateSuperClusterAssociationCollection);
 
   
  // Get the basic cluster collection in the Barrel 
  edm::Handle<reco::BasicClusterCollection> bcBarrelHandle;
  theEvent.getByLabel(bcProducer_, bcBarrelCollection_, bcBarrelHandle);
  if (!bcBarrelHandle.isValid()) {
    edm::LogError("ConverionTrackCandidateProducer") << "Error! Can't get the product "<<bcBarrelCollection_.c_str();
    return;
  }
  reco::BasicClusterCollection clusterCollectionBarrel = *(bcBarrelHandle.product());
  LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer basic cluster collection size  " << clusterCollectionBarrel.size() << "\n";
  
  
  
  // Get the basic cluster collection in the Endcap 
  edm::Handle<reco::BasicClusterCollection> bcEndcapHandle;
  theEvent.getByLabel(bcProducer_, bcEndcapCollection_, bcEndcapHandle);
  if (!bcEndcapHandle.isValid()) {
    edm::LogError("CoonversionTrackCandidateProducer") << "Error! Can't get the product "<<bcEndcapCollection_.c_str();
    return;
  }
  reco::BasicClusterCollection clusterCollectionEndcap = *(bcEndcapHandle.product());
  LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer basic cluster collection size  " << clusterCollectionEndcap.size() << "\n";
  
  
  // Get the Super Cluster collection in the Barrel
  Handle<reco::SuperClusterCollection> scBarrelHandle;
  theEvent.getByLabel(scHybridBarrelProducer_,scHybridBarrelCollection_,scBarrelHandle);
  if (!scBarrelHandle.isValid()) {
    edm::LogError("CoonversionTrackCandidateProducer") << "Error! Can't get the product "<<scHybridBarrelCollection_.c_str();
    return;
  }
  reco::SuperClusterCollection scBarrelCollection = *(scBarrelHandle.product());
  LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer barrel  SC collection size  " << scBarrelCollection.size() << "\n";
  
  // Get the Super Cluster collection in the Endcap
  Handle<reco::SuperClusterCollection> scEndcapHandle;
  theEvent.getByLabel(scIslandEndcapProducer_,scIslandEndcapCollection_,scEndcapHandle);
  if (!scEndcapHandle.isValid()) {
    edm::LogError("CoonversionTrackCandidateProducer") << "Error! Can't get the product "<<scIslandEndcapCollection_.c_str();
    return;
  }
  reco::SuperClusterCollection scEndcapCollection = *(scEndcapHandle.product());
  LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Endcap SC collection size  " << scEndcapCollection.size() << "\n";
  

  /////////////
  vecOfSCRefForOutIn.clear();
  vecOfSCRefForInOut.clear();

  buildCollections(scBarrelHandle, bcBarrelHandle, *outInTrackCandidate_p,*inOutTrackCandidate_p,vecOfSCRefForOutIn,vecOfSCRefForInOut );
  buildCollections(scEndcapHandle, bcEndcapHandle, *outInTrackCandidate_p,*inOutTrackCandidate_p,vecOfSCRefForOutIn,vecOfSCRefForInOut );




  LogDebug("ConversionTrackCandidateProducer")  << "  ConversionTrackCandidateProducer vecOfSCRefForOutIn size " << vecOfSCRefForOutIn.size() << " vecOfSCRefForInOut size " << vecOfSCRefForInOut.size()  << "\n"; 
  


  // put all products in the event
 // Barrel
 LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Putting in the event " << (*outInTrackCandidate_p).size() << " Out In track Candidates " << "\n";
 edm::LogInfo("ConversionTrackCandidateProducer") << "Number of outInTrackCandidates: " <<  (*outInTrackCandidate_p).size() << "\n";
 const edm::OrphanHandle<TrackCandidateCollection> refprodOutInTrackC = theEvent.put( outInTrackCandidate_p, OutInTrackCandidateCollection_ );
 LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer  refprodOutInTrackC size  " <<  (*(refprodOutInTrackC.product())).size()  <<  "\n";
 //
 LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Putting in the event  " << (*inOutTrackCandidate_p).size() << " In Out track Candidates " <<  "\n";
 edm::LogInfo("ConversionTrackCandidateProducer") << "Number of inOutTrackCandidates: " <<  (*inOutTrackCandidate_p).size() << "\n";
 const edm::OrphanHandle<TrackCandidateCollection> refprodInOutTrackC = theEvent.put( inOutTrackCandidate_p, InOutTrackCandidateCollection_ );
 LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer  refprodInOutTrackC size  " <<  (*(refprodInOutTrackC.product())).size()  <<  "\n";
  

 LogDebug("ConversionTrackCandidateProducer") << " ConversionTrackCandidateProduce Going to fill association maps  " <<  "\n";
 for (unsigned int i=0;i< vecOfSCRefForOutIn.size(); ++i) {
   outInAssoc_p->insert(edm::Ref<TrackCandidateCollection>(refprodOutInTrackC,i), vecOfSCRefForOutIn[i]  );
 }
 for (unsigned int i=0;i< vecOfSCRefForInOut.size(); ++i) {
   inOutAssoc_p->insert(edm::Ref<TrackCandidateCollection>(refprodInOutTrackC,i), vecOfSCRefForInOut[i]  );
 }
 

  
 LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Putting in the event   OutIn track - SC association: size  " <<  (*outInAssoc_p).size() << "\n";  
 theEvent.put( outInAssoc_p, OutInTrackSuperClusterAssociationCollection_);
 
 LogDebug("ConversionTrackCandidateProducer") << "ConversionTrackCandidateProducer Putting in the event   InOut track - SC association: size  " <<  (*inOutAssoc_p).size() << "\n";  
 theEvent.put( inOutAssoc_p, InOutTrackSuperClusterAssociationCollection_);


  
}


void ConversionTrackCandidateProducer::buildCollections( const edm::Handle<reco::SuperClusterCollection> & scHandle,
                                                         const edm::Handle<reco::BasicClusterCollection> & bcHandle,
							 TrackCandidateCollection& outInTrackCandidates,
							 TrackCandidateCollection& inOutTrackCandidates,
							 std::vector<edm::Ref<reco::SuperClusterCollection> >& vecRecOI,
							 std::vector<edm::Ref<reco::SuperClusterCollection> >& vecRecIO)

{


  //  Loop over SC in the barrel and reconstruct converted photons
  
  int lSC=0; // local index for getting the right Ref to supercluster 
  reco::BasicClusterCollection clusterCollection = *(bcHandle.product());
  reco::SuperClusterCollection scCollection = *(scHandle.product());
  reco::SuperClusterCollection::iterator aClus;

  for(aClus = scCollection.begin(); aClus != scCollection.end(); ++aClus) {
  
   LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer  SC eta " <<  aClus->eta() << " phi " <<  aClus->phi() <<  " Energy " <<  aClus->energy() << "\n";

    theOutInSeedFinder_->setCandidate(*aClus);
    theOutInSeedFinder_->makeSeeds(  clusterCollection );


    
    std::vector<Trajectory> theOutInTracks= theOutInTrackFinder_->tracks(theOutInSeedFinder_->seeds(),  outInTrackCandidates);    

    theInOutSeedFinder_->setCandidate(*aClus);
    theInOutSeedFinder_->setTracks(  theOutInTracks );   
    theInOutSeedFinder_->makeSeeds(  clusterCollection);
    
    std::vector<Trajectory> theInOutTracks= theInOutTrackFinder_->tracks(theInOutSeedFinder_->seeds(),  inOutTrackCandidates); 


    // Debug
   LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer  theOutInTracks.size() " << theOutInTracks.size() << " theInOutTracks.size() " << theInOutTracks.size() <<  " Event pointer to out in track size barrel " << outInTrackCandidates.size() << " in out track size " << inOutTrackCandidates.size() <<   "\n";


    //////////// Fill vectors of Ref to SC to be used for the Track-SC association
    reco::SuperClusterRef scRefOutIn(reco::SuperClusterRef(scHandle, lSC));
    reco::SuperClusterRef scRefInOut(reco::SuperClusterRef(scHandle, lSC));
    lSC++;

    for (std::vector<Trajectory>::const_iterator it = theOutInTracks.begin(); it !=  theOutInTracks.end(); ++it) {
      vecOfSCRefForOutIn.push_back( scRefOutIn );
           
     LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Barrel OutIn Tracks Number of hits " << (*it).foundHits() << "\n"; 
    }

    for (std::vector<Trajectory>::const_iterator it = theInOutTracks.begin(); it !=  theInOutTracks.end(); ++it) {
      vecOfSCRefForInOut.push_back( scRefInOut );

     LogDebug("ConversionTrackCandidateProducer")  << "ConversionTrackCandidateProducer Barrel InOut Tracks Number of hits " << (*it).foundHits() << "\n"; 
    }




    

  }



}

