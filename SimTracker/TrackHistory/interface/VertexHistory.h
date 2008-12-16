#ifndef VertexHistory_h
#define VertexHistory_h

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimTracker/TrackHistory/interface/HistoryBase.h"
#include "SimTracker/VertexAssociation/interface/VertexAssociatorBase.h"

//! This class traces the simulated and generated history of a given track.
class VertexHistory : public HistoryBase
{

public:

    //! Constructor by pset.
    /* Creates a VertexHistory with association given by pset.

       /param[in] pset with the configuration values
    */
    VertexHistory(const edm::ParameterSet &);

    //! Pre-process event information (for accessing reconstruction information)
    void newEvent(const edm::Event &, const edm::EventSetup &);

    //! Evaluate track history using a TrackingParticleRef.
    /* Return false when the history cannot be determined upto a given depth.
       If not depth is pass to the function no restriction are apply to it.

       /param[in] trackingVertexRef of a simulated track
       /param[in] depth of the vertex history
       /param[out] boolean that is true when history can be determined
    */
    bool evaluate(TrackingVertexRef tvr)
    {
        return HistoryBase::evaluate(tvr);
    }

    //! Evaluate reco::Vertex history using a given association.
    /* Return false when the track association is not possible (fake track).

       /param[in] VertexRef to a reco::track
       /param[out] boolean that is false when a fake track is detected
    */
    bool evaluate (reco::VertexRef);

private:

    bool bestMatchByMaxValue_;

    edm::InputTag trackProducer_;

    edm::InputTag vertexProducer_;

    edm::InputTag trackingTruth_;

    std::string trackAssociator_;

    std::string vertexAssociator_;

    reco::VertexRecoToSimCollection association_;

    TrackingVertexRef match ( reco::VertexRef );  
};

#endif
