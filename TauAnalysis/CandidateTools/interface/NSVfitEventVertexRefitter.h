#ifndef TauAnalysis_CandidateTools_NSVfitEventVertexRefitter_h
#define TauAnalysis_CandidateTools_NSVfitEventVertexRefitter_h

/** \class NSVfitEventVertexRefitter
 *
 * Class to refit position of primary event vertex,
 * excluding the tracks associated to tau decay products
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitEventVertexRefitter.h,v 1.1 2011/03/23 17:46:39 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

class NSVfitEventVertexRefitter
{
 public:
  NSVfitEventVertexRefitter(const edm::ParameterSet&);
  ~NSVfitEventVertexRefitter();

  void beginEvent(const edm::Event&, const edm::EventSetup&);

  /// Refit the primary vertex after subtracting the tau associated tracks.
  TransientVertex refit(const reco::Vertex*, const std::vector<const reco::Track*>* = 0);

  /// Fit decay vertex of three-prong tau
  TransientVertex fitSecondaryVertex(const std::vector<const reco::Track*>&) const;

 private:
  edm::InputTag srcBeamSpot_;

  const reco::BeamSpot* beamSpot_;

  const TransientTrackBuilder* trackBuilder_;

  const KalmanVertexFitter* vertexFitAlgorithm_;

  unsigned minNumTracksRefit_;
};

#endif
