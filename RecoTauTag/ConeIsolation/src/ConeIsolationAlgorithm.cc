#include "RecoTauTag/ConeIsolation/interface/ConeIsolationAlgorithm.h"
using namespace std;
using namespace reco;
using namespace edm;

ConeIsolationAlgorithm::ConeIsolationAlgorithm(void)
{ }

ConeIsolationAlgorithm::ConeIsolationAlgorithm(const ParameterSet & parameters)
{
  //FIXME: use unsigned int where needed
  m_cutPixelHits     = parameters.getParameter<int>("MinimumNumberOfPixelHits");    // not used
  m_cutTotalHits     = parameters.getParameter<int>("MinimumNumberOfHits");
  m_cutMaxTIP        = parameters.getParameter<double>("MaximumTransverseImpactParameter");
  m_cutMinPt         = parameters.getParameter<double>("MinimumTransverseMomentum");
  m_cutMaxChiSquared = parameters.getParameter<double>("MaximumChiSquared");
  dZ_vertex          = parameters.getParameter<double>("DeltaZetTrackVertex");      //  to be modified
  useVertexConstrain_ = parameters.getParameter<bool>("useVertex");

  matching_cone      = parameters.getParameter<double>("MatchingCone");
  signal_cone        = parameters.getParameter<double>("SignalCone");
  isolation_cone     = parameters.getParameter<double>("IsolationCone"); 
  pt_min_isolation   = parameters.getParameter<double>("MinimumTransverseMomentumInIsolationRing"); 
  pt_min_leadTrack   = parameters.getParameter<double>("MinimumTransverseMomentumLeadingTrack"); 
  n_tracks_isolation_ring = parameters.getParameter<int>("MaximumNumberOfTracksIsolationRing"); 
  
  useFixedSizeCone = parameters.getParameter<bool>("UseFixedSizeCone"); 
  variableConeParameter = parameters.getParameter<double>("VariableConeParameter");
  variableMaxCone = parameters.getParameter<double>("VariableMaxCone");
  variableMinCone = parameters.getParameter<double>("VariableMinCone");
}

pair<float,IsolatedTauTagInfo> ConeIsolationAlgorithm::tag(const JetTracksAssociationRef & jetTracks, const Vertex & pv) 
{
  const edm::RefVector<reco::TrackCollection> & tracks = jetTracks->second;
  edm::RefVector<reco::TrackCollection> myTracks;

  // Selection of the Tracks
  float z_pv = pv.z();
  for(edm::RefVector<reco::TrackCollection>::const_iterator it = tracks.begin(); it!= tracks.end(); ++it)
  {
    if ( (*it)->pt()                                  >  m_cutMinPt                     &&
         (*it)->normalizedChi2()                      <  m_cutMaxChiSquared             &&
         fabs((*it)->dxy(pv.position()))                            <  m_cutMaxTIP                    &&
         (*it)->recHitsSize()                         >= (unsigned int) m_cutTotalHits  &&
         (*it)->hitPattern().numberOfValidPixelHits() >= m_cutPixelHits ) 
    {
      if (useVertexConstrain_ && z_pv > -500.) {
        if (fabs((*it)->dz() - z_pv) < dZ_vertex)
          myTracks.push_back(*it);
      } else
        myTracks.push_back(*it);
    }
  }
  IsolatedTauTagInfo resultExtended(myTracks,jetTracks);

  double r_sigCone = signal_cone;
  double energyJet = jetTracks->first->energy();
  if (not useFixedSizeCone) {
    r_sigCone = std::min(variableMaxCone, variableConeParameter / energyJet);
    r_sigCone = std::max((double)r_sigCone, variableMinCone);
  }

  // now I can use it for the discriminator;
  math::XYZVector jetDir(jetTracks->first->px(), jetTracks->first->py(), jetTracks->first->pz());   
  float discriminator = 0.;
  if (useVertexConstrain_) {
    // In this case all the selected tracks comes from the same vertex, so no need to pass the dZ_vertex requirement to the discriminator 
    discriminator = resultExtended.discriminator(jetDir, matching_cone, r_sigCone, isolation_cone, pt_min_leadTrack, pt_min_isolation,  n_tracks_isolation_ring); 
  } else {
    // In this case the dZ_vertex is used to associate the tracks to the Z_imp parameter of the Leading Track
    discriminator = resultExtended.discriminator(jetDir, matching_cone, r_sigCone, isolation_cone, pt_min_leadTrack, pt_min_isolation,  n_tracks_isolation_ring, dZ_vertex); 
  }

  return std::make_pair(discriminator, resultExtended); 
}
