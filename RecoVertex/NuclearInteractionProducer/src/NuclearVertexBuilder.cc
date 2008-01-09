#include "RecoVertex/NuclearInteractionProducer/interface/NuclearVertexBuilder.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

void NuclearVertexBuilder::build( const reco::TrackRef& primTrack, const reco::TrackRefVector& secTracks) {

       if( secTracks.size() != 0) {
         std::vector<reco::TransientTrack> transientTracks;
         transientTracks.push_back( theTransientTrackBuilder->build(primTrack));
         // get the secondary track with the max number of hits
         for( reco::TrackRefVector::const_iterator tk = secTracks.begin(); tk != secTracks.end(); tk++) { 
                    transientTracks.push_back( theTransientTrackBuilder->build(*tk));
         }
         AdaptiveVertexFitter AVF;
         try {
            TransientVertex tv = AVF.vertex(transientTracks); 
            the_vertex = reco::Vertex(tv);
            LogDebug("NuclearInteractionMaker") << "Vertex build from AdaptiveVertexFitter with " << the_vertex.tracksSize() << " tracks";
         }
         catch(VertexException& exception){
            // AdaptivevertexFitter does not work 
            // Try here a simple intersection between the primary track
            // and the longer secondary track
            FreeTrajectoryState primTraj = getTrajectory(primTrack);

            // get the secondary track with the max number of hits
            unsigned short maxHits = 0; int indice=0;
            for( unsigned short i=0; i < secTracks.size(); i++) { 
                   // Add references to daughters
                   unsigned short nhits = secTracks[i]->numberOfValidHits();
                   if( nhits > maxHits ) { maxHits = nhits; indice=i; }
            }
            FreeTrajectoryState secTraj = getTrajectory(secTracks[indice]);

            // Closest points
            ClosestApproachInRPhi *theApproach = new ClosestApproachInRPhi();
	    bool status = theApproach->calculate(primTraj,secTraj);
            GlobalPoint crossing;
            if( status ) crossing = theApproach->crossingPoint();
            else crossing = GlobalPoint(primTrack->outerPosition().x(),
                                        primTrack->outerPosition().y(),
                                        primTrack->outerPosition().z());
            delete theApproach;

            // Create vertex (creation point)
            the_vertex = reco::Vertex(reco::Vertex::Point(crossing.x(),
                                                          crossing.y(),
                                                          crossing.z()),
                                      reco::Vertex::Error(), 0.,0.,0);
            // TODO : add a weight to the track - use Vertex::add( RefTrack, weight)
            LogDebug("NuclearInteractionMaker") << exception.what() << " - Vertex build from crossing point" << "\n"
                    << "Addition of a primary track and " << secTracks.size() << " secondary tracks to the vertex";
            the_vertex.add(reco::TrackBaseRef( primTrack ));
            for( unsigned short i=0; i < secTracks.size(); i++) { 
                 the_vertex.add(reco::TrackBaseRef( secTracks[i] ));
            }
        }
        catch(cms::Exception& exception){
           throw exception;
        }
     }
     else {
       // if no secondary tracks : vertex position = position of last rechit of the primary track 
       LogDebug("NuclearInteractionMaker") << "Vertex build from the end point of the primary track";
       the_vertex = reco::Vertex(reco::Vertex::Point(primTrack->outerPosition().x(),
                                         primTrack->outerPosition().y(),
                                         primTrack->outerPosition().z()),
                                         reco::Vertex::Error(), 0.,0.,0);
       the_vertex.add(reco::TrackBaseRef( primTrack ));
       for( unsigned short i=0; i < secTracks.size(); i++) { 
                 the_vertex.add(reco::TrackBaseRef( secTracks[i] ));
       }
     }


}
          
FreeTrajectoryState NuclearVertexBuilder::getTrajectory(const reco::TrackRef& track)
{ 
  GlobalPoint position(track->vertex().x(),
                       track->vertex().y(),
                       track->vertex().z());

  GlobalVector momentum(track->momentum().x(),
                        track->momentum().y(),
                        track->momentum().z());

  GlobalTrajectoryParameters gtp(position,momentum,
                                 track->charge(),theMagField);
  
  FreeTrajectoryState fts(gtp);

  return fts; 
} 

