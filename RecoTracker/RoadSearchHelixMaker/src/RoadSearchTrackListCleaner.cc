//
// Package:         RecoTracker/RoadSearchTrackListCleaner
// Class:           RoadSearchHelixMaker
// 
// Description:     TrackListCleaner
//
// Original Author: Steve Wagner, stevew@pizero.colorado.edu
// Created:         Sat Jan 14 22:00:00 UTC 2006
//
// $Author: gutsche $
// $Date: 2007/01/15 22:16:28 $
// $Revision: 1.3 $
//

#include <memory>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>

#include "RecoTracker/RoadSearchHelixMaker/interface/RoadSearchTrackListCleaner.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/src/classes.h"

namespace cms
{

  RoadSearchTrackListCleaner::RoadSearchTrackListCleaner(edm::ParameterSet const& conf) : 
    conf_(conf)
  {
    produces<reco::TrackCollection>();
//    produces<reco::TrackExtraCollection>();
  }


  // Virtual destructor needed.
  RoadSearchTrackListCleaner::~RoadSearchTrackListCleaner() { }  

  // Functions that gets called by framework every event
  void RoadSearchTrackListCleaner::produce(edm::Event& e, const edm::EventSetup& es)
  {
    // retrieve producer name of input SiStripRecHit2DCollection
    std::string trackProducer = conf_.getParameter<std::string>("TrackProducer");
  
    //
    // extract tracker geometry
    //
    edm::ESHandle<TrackerGeometry> theG;
    es.get<TrackerDigiGeometryRecord>().get(theG);

//    using namespace reco;

    // get Inputs 
    edm::Handle<reco::TrackCollection> trackCollection;
    e.getByLabel(trackProducer, trackCollection);

    const reco::TrackCollection tC = *(trackCollection.product());

    std::cout << "Reconstructed "<< tC.size() << " tracks" << std::endl ;

    // Step B: create empty output collection
    std::auto_ptr<reco::TrackCollection> output(new reco::TrackCollection);

  //
  //  no input tracks
  //

    if ( tC.empty() ){
//      LogDebug("RoadSearch") << "Found " << output.size() << " clouds.";
      e.put(output);
      return;  
    }

  //
  //  1 raw cloud - nothing to try merging, but one cloud to duplicate
  //

    if ( 1==tC.size() ){
      for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){
        reco::Track * theTrack = new reco::Track(track->chi2(),
						 (short unsigned)track->ndof(),
						 track->innerPosition(),
						 track->innerMomentum(),
						 track->charge(),
						 track->innerStateCovariance());    
        //fill the TrackCollection
        reco::TrackExtraRef theTrackExtraRef=track->extra();    
        theTrack->setExtra(theTrackExtraRef);    
	theTrack->setHitPattern((*theTrackExtraRef).recHits());
        output->push_back(*theTrack);
        delete theTrack;
      }//end faux loop over tracks
      e.put(output);
      return;
    }  
  //
  //  > 1 track - try merging
  //
    std::vector<int> not_dup; for (unsigned int i=0; i<tC.size(); ++i){not_dup.push_back(1);}
    int i=-1;
    for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){
      i++; std::cout << "Track number "<< i << std::endl ; if (!not_dup[i])continue;
      int j=-1;
      for (reco::TrackCollection::const_iterator track2=tC.begin(); track2!=tC.end(); track2++){
        j++;
        if ((j<=i)||(!not_dup[j])||(!not_dup[i]))continue;
        int noverlap=0;
        for (trackingRecHit_iterator it = track->recHitsBegin();  it != track->recHitsEnd(); it++){
          if ((*it)->isValid()){
            for (trackingRecHit_iterator jt = track2->recHitsBegin();  jt != track2->recHitsEnd(); jt++){
	      if ((*jt)->isValid()){
                if (((*it)->geographicalId()==(*jt)->geographicalId())&&((*it)->localPosition().x()==(*jt)->localPosition().x()))noverlap++;
              }
            }
          }
        }
        float fi=float(noverlap)/float(track->recHitsSize()); float fj=float(noverlap)/float(track2->recHitsSize());
        std::cout << " trk1 trk2 nhits1 nhits2 nover " << i << " " << j << " " << track->recHitsSize() << " " 
                  << track2->recHitsSize() << " " << noverlap << " " << fi << " " << fj  <<std::endl;
        if ((fi>0.66)||(fj>0.66)){
          if (fi<fj){
            not_dup[j]=0; std::cout << " removing 2nd trk in pair " << std::endl;
          }else{
            if (fi>fj){
              not_dup[i]=0; std::cout << " removing 1st trk in pair " << std::endl;
            }else{
              std::cout << " removing worst chisq in pair " << std::endl;
              if (track->chi2() > track2->chi2()){not_dup[i]=0;}else{not_dup[j]=0;}
            }//end fi > or = fj
          }//end fi < fj
        }//end got a duplicate
      }//end track2 loop
    }//end track loop
    i=-1;
    for (reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++){
      i++;  if (!not_dup[i])continue;
        reco::Track * theTrack = new reco::Track(track->chi2(),
						 (short unsigned)track->ndof(),
						 track->innerPosition(),
						 track->innerMomentum(),
						 track->charge(),
						 track->innerStateCovariance());    
      //fill the TrackCollection
      reco::TrackExtraRef theTrackExtraRef=track->extra();    
      theTrack->setExtra(theTrackExtraRef);    
      theTrack->setHitPattern((*theTrackExtraRef).recHits());
      output->push_back(*theTrack);
      delete theTrack;
    }//end faux loop over tracks
    e.put(output);
    return;

  }//end produce

}
