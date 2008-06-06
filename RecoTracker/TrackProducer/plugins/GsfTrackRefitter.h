#ifndef GsfTrackRefitter_h
#define GsfTrackRefitter_h

/** \class GsfTrackRefitter
 *  Refit GSF Tracks. Based on the TrackRefitter.
 */

#include "RecoTracker/TrackProducer/interface/GsfTrackProducerBase.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerAlgorithm.h"

class GsfTrackRefitter : public GsfTrackProducerBase, public edm::EDProducer {
public:

  /// Constructor
  explicit GsfTrackRefitter(const edm::ParameterSet& iConfig);

  /// Implementation of produce method
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  TrackProducerAlgorithm<reco::GsfTrack> theAlgo;
//   enum Constraint { none, momentum, vertex };
//   Constraint constraint_;
};

#endif
