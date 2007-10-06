#ifndef TrackRefitter_h
#define TrackRefitter_h

/** \class TrackRefitter
 *  Refit Tracks: Produce Tracks from TrackCollection. It performs a new final fit on a TrackCollection.
 *
 *  $Date: 2007/07/16 10:07:41 $
 *  $Revision: 1.3 $
 *  \author cerati
 */

#include "RecoTracker/TrackProducer/interface/KfTrackProducerBase.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerAlgorithm.h"

class TrackRefitter : public KfTrackProducerBase, public edm::EDProducer {
public:

  /// Constructor
  explicit TrackRefitter(const edm::ParameterSet& iConfig);

  /// Implementation of produce method
  virtual void produce(edm::Event&, const edm::EventSetup&);

private:
  TrackProducerAlgorithm<reco::Track> theAlgo;
  std::string constraint_;
};

#endif
