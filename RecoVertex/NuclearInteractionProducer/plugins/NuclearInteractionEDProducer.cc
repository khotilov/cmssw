#include "RecoVertex/NuclearInteractionProducer/interface/NuclearInteractionEDProducer.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/NuclearInteractionFwd.h"

#include "RecoVertex/NuclearInteractionProducer/interface/NuclearVertexBuilder.h"
#include "RecoVertex/NuclearInteractionProducer/interface/NuclearLikelihood.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/Framework/interface/EventSetup.h"


NuclearInteractionEDProducer::NuclearInteractionEDProducer(const edm::ParameterSet& iConfig) : 
conf_(iConfig), 
primaryProducer_(iConfig.getParameter<std::string>("primaryProducer")),
seedsProducer_(iConfig.getParameter<std::string>("seedsProducer")),
secondaryProducer_(iConfig.getParameter<std::string>("secondaryProducer")),
minDistFromPrim_(iConfig.getParameter<double>("minDistFromPrimary")) {

  produces<reco::NuclearInteractionCollection>();
}

NuclearInteractionEDProducer::~NuclearInteractionEDProducer()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
NuclearInteractionEDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   /// Get the primary tracks
   edm::Handle<reco::TrackCollection>  primaryTrackCollection;
   iEvent.getByLabel( primaryProducer_, primaryTrackCollection );

   /// Get the primary trajectories (produced by the Refitter)
   edm::Handle< TrajectoryCollection > primaryTrajectoryCollection;
   iEvent.getByLabel( primaryProducer_, primaryTrajectoryCollection );

   /// Get the AssociationMap between primary tracks and trajectories
   edm::Handle< TrajTrackAssociationCollection > refMapH;
   iEvent.getByLabel( primaryProducer_, refMapH );
   const TrajTrackAssociationCollection& refMap = *(refMapH.product());

   /// Get the AssociationMap between seeds and primary trajectories
   edm::Handle<TrajectoryToSeedsMap>  nuclMapH;
   iEvent.getByLabel(seedsProducer_, nuclMapH);
   const TrajectoryToSeedsMap& nuclMap = *(nuclMapH.product());

   /// Get the secondary tracks
   edm::Handle<reco::TrackCollection>  secondaryTrackCollection;
   iEvent.getByLabel( secondaryProducer_, secondaryTrackCollection );

   /// Definition of the output
   std::auto_ptr<reco::NuclearInteractionCollection> theNuclearInteractions(new reco::NuclearInteractionCollection);

   typedef edm::Ref<TrajectoryCollection> TrajectoryRef;

   /// Loop on all primary trajectories
   for(unsigned int i = 0; i < primaryTrajectoryCollection->size() ; i++) {

         TrajectoryRef  trajRef( primaryTrajectoryCollection, i );

         /// 1. Get the primary track from the trajectory
         TrajTrackAssociationCollection::const_iterator itPrimTrack = refMap.find( trajRef );
         if( itPrimTrack == refMap.end() || (itPrimTrack->val).isNull() ) continue;
         const reco::TrackRef& primary_track = itPrimTrack->val;

         /// 2. Get the seeds from the trajectory
         TrajectoryToSeedsMap::const_iterator itSeeds = nuclMap.find( trajRef );
         if( itSeeds == nuclMap.end() || (itSeeds->val).isNull()) continue; 
         const TrajectorySeedRefVector& seeds = itSeeds->val;

         /// 3. Get the secondary tracks
         reco::TrackRefVector secondary_tracks;
         for( unsigned int k=0; k < secondaryTrackCollection->size(); k++) {
                     reco::TrackRef currentTrk(secondaryTrackCollection, k);
                     if( isInside( currentTrk, seeds ) ) secondary_tracks.push_back(currentTrk);
         }
              
         /// 4. Get the vertex and the likelihood
         vertexBuilder->build(primary_track, secondary_tracks);
         likelihoodCalculator->calculate( vertexBuilder->getVertex() );

         reco::NuclearInteraction nuclInter(seeds, vertexBuilder->getVertex(), likelihoodCalculator->result() );
         theNuclearInteractions->push_back( nuclInter );

         std::ostringstream  str;
         print(str, nuclInter, vertexBuilder);
         edm::LogInfo("NuclearInteractionMaker") << str.str();

   }

   LogDebug("NuclearInteractionMaker") << "End of NuclearInteractionMaker - Number of nuclear interactions found :" << theNuclearInteractions->size();
   iEvent.put(theNuclearInteractions);
}

// ------------ method called once each job just before starting event loop  ------------
void
NuclearInteractionEDProducer::beginJob(const edm::EventSetup& es)
{
   /// Get magnetic field
   edm::ESHandle<MagneticField> magField;
   es.get<IdealMagneticFieldRecord>().get(magField);

   edm::ESHandle<TransientTrackBuilder> builder;
   es.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);

   vertexBuilder = std::auto_ptr< NuclearVertexBuilder >(new NuclearVertexBuilder( magField.product(), builder.product(), minDistFromPrim_) );
   likelihoodCalculator = std::auto_ptr< NuclearLikelihood >(new NuclearLikelihood);

}

void  NuclearInteractionEDProducer::endJob() {}

// ------ method used to check whether the seed of a track belong to the vector of seeds --
bool NuclearInteractionEDProducer::isInside( const reco::TrackRef& track, const TrajectorySeedRefVector& seeds) {
    unsigned int seedKey = track->seedRef().key();
    for (unsigned int i=0; i< seeds.size(); i++) { if( seeds[i].key() == seedKey ) return true; }
    return false;
}

// -- print out
void print(std::ostringstream& out, const reco::NuclearInteraction& nucl, const std::auto_ptr< NuclearVertexBuilder >& builder) {
   out<<"Nuclear Interaction with vertex position : (";
   out<< nucl.vertex().position().x() << " , "
      << nucl.vertex().position().y() << " , "
      << nucl.vertex().position().z() << ")";
   out<<"\tLikelihood : " << nucl.likelihood() << std::endl;
   out<<"\tPrimary Track : Pt = " << nucl.primaryTrack()->pt() << "  - Nhits = "
      << nucl.primaryTrack()->numberOfValidHits() << std::endl;
   out << "\tNumber of seeds : " << nucl.seedsSize() << std::endl;
   out << "\tNumber of secondary Tracks : " << nucl.secondaryTracksSize() << std::endl;
   int it=0;
   for( reco::NuclearInteraction::trackRef_iterator itr_=nucl.secondaryTracks_begin(); itr_ != nucl.secondaryTracks_end(); itr_++, it++) {
                out << "\t\t Secondary track " << it << " : Pt = " << (*itr_)->pt() 
                    << " - Nhits = " << (*itr_)->numberOfValidHits()
                    << " - Dist = " << builder->distanceOfClosestApproach(it) 
                    << " - chi2 = " << (*itr_)->normalizedChi2() << std::endl;
      }
   out << "----------------" << std::endl;
}
