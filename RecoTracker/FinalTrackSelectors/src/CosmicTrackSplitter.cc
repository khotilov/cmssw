#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventPrincipal.h" 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"

// added by me

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h" 
#include "DataFormats/GeometrySurface/interface/Surface.h" 
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h" 
#include "DataFormats/Math/interface/Vector.h" 
#include "DataFormats/Math/interface/Error.h" 
#include "TrackingTools/TrajectoryState/interface/CopyUsingClone.h" 
#include "RecoVertex/VertexTools/interface/PerigeeLinearizedTrackState.h" 
#include "Alignment/ReferenceTrajectories/interface/TrajectoryFactoryPlugin.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <boost/regex.hpp>

/**
 * Configurables:
 *
 *   Generic: 
 *     tracks                  = InputTag of a collection of tracks to read from
 *     minimumHits             = Minimum hits that the output TrackCandidate must have to be saved
 *     replaceWithInactiveHits = instead of discarding hits, replace them with a invalid "inactive" hits, 
 *                               so multiple scattering is accounted for correctly.
 *     stripFrontInvalidHits   = strip invalid hits at the beginning of the track
 *     stripBackInvalidHits    = strip invalid hits at the end of the track
 *     stripAllInvalidHits     = remove ALL invald hits (might be a problem for multiple scattering, use with care!)
 *
 *   Per structure: 
 *      commands = list of commands, to be applied in order as they are written
 *      commands can be:
 *          "keep XYZ"  , "drop XYZ"    (XYZ = PXB, PXE, TIB, TID, TOB, TEC)
 *          "keep XYZ l", "drop XYZ n"  (XYZ as above, n = layer, wheel or disk = 1 .. 6 ; positive and negative are the same )
 *
 *   Individual modules:
 *     detsToIgnore        = individual list of detids on which hits must be discarded
 */
namespace reco { namespace modules {
class CosmicTrackSplitter : public edm::EDProducer {
    public:
       CosmicTrackSplitter(const edm::ParameterSet &iConfig) ; 
       virtual void produce(edm::Event &iEvent, const edm::EventSetup &iSetup) ;

    private:
	edm::InputTag tracks_;
    int totalTracks_;
	size_t minimumHits_;

       bool replaceWithInactiveHits_;
       bool stripFrontInvalidHits_;
       bool stripBackInvalidHits_;
       bool stripAllInvalidHits_;

       std::vector<uint32_t> detsToIgnore_;

       edm::ESHandle<TrackerGeometry> theGeometry;
       edm::ESHandle<MagneticField>   theMagField;

       TrackCandidate makeCandidate(const reco::Track &tk, std::vector<TrackingRecHit *>::iterator hitsBegin, std::vector<TrackingRecHit *>::iterator hitsEnd) ;
       
}; // class


CosmicTrackSplitter::CosmicTrackSplitter(const edm::ParameterSet &iConfig) :
    tracks_(iConfig.getParameter<edm::InputTag>("tracks")),
    minimumHits_(iConfig.getParameter<uint32_t>("minimumHits")),
    replaceWithInactiveHits_(iConfig.getParameter<bool>("replaceWithInactiveHits")),
    stripFrontInvalidHits_(iConfig.getParameter<bool>("stripFrontInvalidHits")),
    stripBackInvalidHits_( iConfig.getParameter<bool>("stripBackInvalidHits") ),
    stripAllInvalidHits_(  iConfig.getParameter<bool>("stripAllInvalidHits")  ),
    detsToIgnore_( iConfig.getParameter<std::vector<uint32_t> >("detsToIgnore") )
{
    // sanity check 
    if (stripAllInvalidHits_ && replaceWithInactiveHits_) {
        throw cms::Exception("Configuration") << "Inconsistent Configuration: you can't set both 'stripAllInvalidHits' and 'replaceWithInactiveHits' to true\n";
    }

	// sort detids to ignore
    std::sort(detsToIgnore_.begin(), detsToIgnore_.end());

	totalTracks_ = 0;
	
    // issue the produce<>
    produces<TrackCandidateCollection>();
}

void 
CosmicTrackSplitter::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) 
{
	std::cout << "Entering Producer" << std::endl;
	
    // read with View, so we can read also a TrackRefVector
    edm::Handle<std::vector<reco::Track> > tracks;
    iEvent.getByLabel(tracks_, tracks);

    // read from EventSetup
    iSetup.get<TrackerDigiGeometryRecord>().get(theGeometry);
    iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
	
    // prepare output collection
    std::auto_ptr<TrackCandidateCollection> output(new TrackCandidateCollection());
    output->reserve(tracks->size());
    
    // working area and tools
    std::vector<TrackingRecHit *> hits;

    //std::cout << "CosmicTrackSplitter: loop on tracks" << std::endl;
	totalTracks_ = totalTracks_ + tracks->size();
    // loop on tracks
    for (std::vector<reco::Track>::const_iterator itt = tracks->begin(), edt = tracks->end(); itt != edt; ++itt) {
        hits.clear(); // extra safety
		// try to find distance of closest approach
		
		reco::TransientTrack tt( *(itt), theMagField.product() );//, theGeometry);
		FreeTrajectoryState fts = tt.initialFreeState();
		TSCPBuilderNoMaterial tscpBuilder;
		TrajectoryStateClosestToPoint tsAtClosestApproach     = tscpBuilder(fts,GlobalPoint(0,0,0));//as in TrackProducerAlgorithm
		GlobalPoint v = tsAtClosestApproach.theState().position();
		GlobalVector p = tsAtClosestApproach.theState().momentum();
		//std::cout << "DCA: " << v << std::endl;
		
		
		// LOOP TWICE, ONCE FOR TOP AND ONCE FOR BOTTOM
		for (int i = 0; i < 2; ++i){
			hits.clear(); // extra safety
			int usedHitCtr = 0;
			//std::cout << "   loop on hits of track #" << (itt - tracks->begin()) << std::endl;
			for (trackingRecHit_iterator ith = itt->recHitsBegin(), edh = itt->recHitsEnd(); ith != edh; ++ith) {
				const TrackingRecHit * hit = ith->get(); // ith is an iterator on edm::Ref to rechit
				//std::cout << "         hit number " << (ith - itt->recHitsBegin()) << std::endl;
				// let's look at valid hits
				if (hit->isValid()) { 
					//std::cout << "            valid, detid = " << hit->geographicalId().rawId() << std::endl;
					DetId detid = hit->geographicalId();
					
					if ((detid.det() == DetId::Tracker)&&((detid.subdetId() == 3)||(detid.subdetId() == 5))) {  // check for tracker hits
						//std::cout << "            valid, tracker " << std::endl;
						bool  verdict = true;
						
						//trying to get the global position of the hit
						//const GeomDetUnit* geomDetUnit =  theGeometry->idToDetUnit( detid ).;
						
						const GlobalPoint pos =  theGeometry->idToDetUnit( detid )->surface().toGlobal(hit->localPosition());
						//std::cout << "hit pos: " << pos << ", dca pos: " << v << std::endl;
						
						// top half
						if ((i == 0)&&(pos.y() < v.y())){
						//if ((i == 0)&&(pos.y() < 0)){
							verdict = false;
							//std::cout << "tophalf" << std::endl;
						}
						// bottom half
						if ((i == 1)&&(pos.y() >= v.y())){
						//if ((i == 1)&&(pos.y() >= 0)){
							verdict = false;
							//std::cout << "bottomhalf" << std::endl;
						}
						
						// if the hit is good, check again at module level
						if ( verdict  && std::binary_search(detsToIgnore_.begin(), detsToIgnore_.end(), detid.rawId())) {
							verdict = false;
						}
						
						//std::cout << "                   verdict after module list: " << (verdict ? "ok" : "no") << std::endl;
						if (verdict == true) {
							// just copy the hit
							hits.push_back(hit->clone());
							usedHitCtr++;
						} 
						else {
							// still, if replaceWithInactiveHits is true we have to put a new hit
							if (replaceWithInactiveHits_) {
								hits.push_back(new InvalidTrackingRecHit(detid, TrackingRecHit::inactive));
							} 
						}
					} 
					else { // just copy non tracker hits
						hits.push_back(hit->clone());
					}
				} 
				else {
					if (!stripAllInvalidHits_) {
						hits.push_back(hit->clone());
					} 
				} // is valid hit
				//std::cout << "         end of hit " << (ith - itt->recHitsBegin()) << std::endl;
			} // loop on hits
			//std::cout << "   end of loop on hits of track #" << (itt - tracks->begin()) << std::endl;
			
			std::vector<TrackingRecHit *>::iterator begin = hits.begin(), end = hits.end();
			
			//std::cout << "   selected " << hits.size() << " hits " << std::endl;
			
			// strip invalid hits at the beginning
			if (stripFrontInvalidHits_) {
				while ( (begin != end) && ( (*begin)->isValid() == false ) ) ++begin;
			}
			
			//std::cout << "   after front stripping we have " << (end - begin) << " hits " << std::endl;
			
			// strip invalid hits at the end
			if (stripBackInvalidHits_ && (begin != end)) {
				--end;
				while ( (begin != end) && ( (*end)->isValid() == false ) ) --end;
				++end;
			}
			
			//std::cout << "   after back stripping we have " << (end - begin) << " hits " << std::endl;
			
			// if we still have some hits
			//if ((end - begin) >= int(minimumHits_)) {
			if ( usedHitCtr >= int(minimumHits_)) {
				output->push_back( makeCandidate ( *itt, begin, end ) );
				//std::cout << "we made a candidate of " << hits.size() << " hits!" << std::endl;
			} 
			// now delete the hits not used by the candidate
			for (begin = hits.begin(), end = hits.end(); begin != end; ++begin) {
				if (*begin) delete *begin;
			} 
			std::cout << "loop: " << i << " has " << usedHitCtr << " active hits and " << hits.size() << " total hits..." << std::endl;
			hits.clear();
		} // loop twice for top and bottom
    } // loop on tracks
	//std::cout << "totalTracks_ = " << totalTracks_ << std::endl;
    iEvent.put(output);
}

TrackCandidate
CosmicTrackSplitter::makeCandidate(const reco::Track &tk, std::vector<TrackingRecHit *>::iterator hitsBegin, std::vector<TrackingRecHit *>::iterator hitsEnd) {
	
	//std::cout << "Making a candidate!" << std::endl;
	
    TrajectoryStateTransform transform;
    PropagationDirection   pdir = tk.seedDirection();
    PTrajectoryStateOnDet *state;
    if ( pdir == anyDirection ) throw cms::Exception("UnimplementedFeature") << "Cannot work with tracks that have 'anyDirecton' \n";
    if ( (pdir == alongMomentum) == ( tk.p() >= tk.outerP() ) ) {
        // use inner state
        TrajectoryStateOnSurface originalTsosIn(transform.innerStateOnSurface(tk, *theGeometry, &*theMagField));
        state = transform.persistentState( originalTsosIn, DetId(tk.innerDetId()) );
    } else { 
        // use outer state
        TrajectoryStateOnSurface originalTsosOut(transform.outerStateOnSurface(tk, *theGeometry, &*theMagField));
        state = transform.persistentState( originalTsosOut, DetId(tk.outerDetId()) );
    }

    TrajectorySeed seed(*state, TrackCandidate::RecHitContainer(), pdir);
 
    TrackCandidate::RecHitContainer ownHits;
    ownHits.reserve(hitsEnd - hitsBegin);
    for ( ; hitsBegin != hitsEnd; ++hitsBegin) { ownHits.push_back( *hitsBegin ); }

    TrackCandidate cand(ownHits, seed, *state, tk.seedRef());
    delete state;

    //std::cout << "   dumping the hits now: " << std::endl;
    //for (TrackCandidate::range hitR = cand.recHits(); hitR.first != hitR.second; ++hitR.first) {
    //    std::cout << "     hit detid = " << hitR.first->geographicalId().rawId() <<
    //        ", type  = " << typeid(*hitR.first).name() << std::endl;
    //}

    return cand;
}

}} //namespaces


// ========= MODULE DEF ==============
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using reco::modules::CosmicTrackSplitter;
DEFINE_FWK_MODULE(CosmicTrackSplitter);
