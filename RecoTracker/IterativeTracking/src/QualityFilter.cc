#include "RecoTracker/IterativeTracking/interface/QualityFilter.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
//
// class decleration
//

using namespace edm;
using namespace reco;
using namespace std;


QualityFilter::QualityFilter(const edm::ParameterSet& iConfig)
{

  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  produces<std::vector<Trajectory> >();
  produces<TrajTrackAssociationCollection>();

  tkTag  = iConfig.getParameter<edm::InputTag>("recTracks");
  string tkQuality = iConfig.getParameter<string>("TrackQuality");

  if (tkQuality=="highPurity") trackQuality_=TrackBase::highPurity;
  if (tkQuality=="tight") trackQuality_=TrackBase::tight;
  if (tkQuality=="loose") trackQuality_=TrackBase::loose;
}


QualityFilter::~QualityFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
QualityFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
  auto_ptr<TrackCollection> selTracks(new TrackCollection);
  auto_ptr<TrackingRecHitCollection> selHits(new TrackingRecHitCollection);
  auto_ptr<TrackExtraCollection> selTrackExtras(new TrackExtraCollection);
  auto_ptr<vector<Trajectory> > outputTJ(new vector<Trajectory> );
  auto_ptr<TrajTrackAssociationCollection> trajTrackMap( new TrajTrackAssociationCollection() );
  
  
  TrackExtraRefProd rTrackExtras = iEvent.getRefBeforePut<TrackExtraCollection>();
  TrackingRecHitRefProd rHits = iEvent.getRefBeforePut<TrackingRecHitCollection>();

  Handle<std::vector<Trajectory> > TrajectoryCollection;
  Handle<TrajTrackAssociationCollection> assoMap;
  
  Handle<TrackCollection> tkCollection;  
  iEvent.getByLabel( tkTag, tkCollection);
  const reco::TrackCollection*  tC = tkCollection.product();
  TrackCollection::const_iterator itxc;
  TrackCollection::const_iterator firstTrack = tC->begin();
  TrackCollection::const_iterator lastTrack = tC->end();

    
  iEvent.getByLabel(tkTag,TrajectoryCollection);
  iEvent.getByLabel(tkTag,assoMap);



  TrajTrackAssociationCollection::const_iterator it = assoMap->begin();
  TrajTrackAssociationCollection::const_iterator lastAssoc = assoMap->end();
  for( ; it != lastAssoc; ++it ) {
    const Ref<vector<Trajectory> > traj = it->key;
    const reco::TrackRef itc = it->val;
    bool goodTk = (itc->quality(trackQuality_));
 
    
    if (goodTk){
      Track track =(*itc);
      //tracks and trajectories
      selTracks->push_back( track );
      outputTJ->push_back( *traj );
      //TRACKING HITS
      trackingRecHit_iterator irhit   =(*itc).recHitsBegin();
      trackingRecHit_iterator lasthit =(*itc).recHitsEnd();
      for (; irhit!=lasthit; ++irhit)
	selHits->push_back((*irhit)->clone() );
          
    }

  }

  //PUT TRACKING HITS IN THE EVENT
  OrphanHandle<TrackingRecHitCollection> theRecoHits = iEvent.put(selHits );
  
  //PUT TRACK EXTRA IN THE EVENT
  unsigned nTracks = selTracks->size();
  selTrackExtras->reserve(nTracks);
  unsigned hits=0;

  for ( unsigned index = 0; index<nTracks; ++index ) { 

    reco::Track& aTrack = selTracks->at(index);
    TrackExtra aTrackExtra(aTrack.outerPosition(),
			   aTrack.outerMomentum(),
			   aTrack.outerOk(),
			   aTrack.innerPosition(),
			   aTrack.innerMomentum(),
			   aTrack.innerOk(),
			   aTrack.outerStateCovariance(),
			   aTrack.outerDetId(),
			   aTrack.innerStateCovariance(),
			   aTrack.innerDetId(),
			   aTrack.seedDirection(),
			   aTrack.seedRef());
	

    //unsigned nHits = aTrack.numberOfValidHits();
    unsigned nHits = aTrack.recHitsSize();
    for ( unsigned int ih=0; ih<nHits; ++ih) {
      aTrackExtra.add(TrackingRecHitRef(theRecoHits,hits++));
    }
    selTrackExtras->push_back(aTrackExtra);
  }

  //CORRECT REF TO TRACK
  OrphanHandle<TrackExtraCollection> theRecoTrackExtras = iEvent.put(selTrackExtras); 
  for ( unsigned index = 0; index<nTracks; ++index ) { 
    const reco::TrackExtraRef theTrackExtraRef(theRecoTrackExtras,index);
    (selTracks->at(index)).setExtra(theTrackExtraRef);
  }

  //TRACKS AND TRAJECTORIES
  OrphanHandle<TrackCollection> theRecoTracks = iEvent.put(selTracks);
  OrphanHandle<vector<Trajectory> > theRecoTrajectories = iEvent.put(outputTJ);

  //TRACKS<->TRAJECTORIES MAP 
  nTracks = theRecoTracks->size();
  for ( unsigned index = 0; index<nTracks; ++index ) { 
    Ref<vector<Trajectory> > trajRef( theRecoTrajectories, index );
    Ref<TrackCollection>    tkRef( theRecoTracks, index );
    trajTrackMap->insert(trajRef,tkRef);
  }
  //MAP IN THE EVENT
  iEvent.put( trajTrackMap );
}

// ------------ method called once each job just before starting event loop  ------------
void 
QualityFilter::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
QualityFilter::endJob() {
}


