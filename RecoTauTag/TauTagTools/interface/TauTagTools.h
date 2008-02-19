#ifndef RecoTauTag_TauTagTools_TauTagTools_h
#define RecoTauTag_TauTagTools_TauTagTools_h

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Math/GenVector/VectorUtil.h" 

#include "RecoTauTag/TauTagTools/interface/ECALBounds.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"

using namespace std;
using namespace reco;
using namespace edm;

namespace TauTagTools{
  TrackRefVector filteredTracksByNumTrkHits(TrackRefVector theInitialTracks, int tkminTrackerHitsn);
  TrackRefVector filteredTracks(TrackRefVector theInitialTracks,double tkminPt,int tkminPixelHitsn,int tkminTrackerHitsn,double tkmaxipt,double tkmaxChi2, Vertex pV);
  TrackRefVector filteredTracks(TrackRefVector theInitialTracks,double tkminPt,int tkminPixelHitsn,int tkminTrackerHitsn,double tkmaxipt,double tkmaxChi2,double tktorefpointmaxDZ,Vertex pV,double refpoint_Z);

  PFCandidateRefVector filteredPFChargedHadrCandsByNumTrkHits(PFCandidateRefVector theInitialPFCands, int ChargedHadrCand_tkminTrackerHitsn);
  PFCandidateRefVector filteredPFChargedHadrCands(PFCandidateRefVector theInitialPFCands,double ChargedHadrCand_tkminPt,int ChargedHadrCand_tkminPixelHitsn,int ChargedHadrCand_tkminTrackerHitsn,double ChargedHadrCand_tkmaxipt,double ChargedHadrCand_tkmaxChi2, Vertex pV);
  PFCandidateRefVector filteredPFChargedHadrCands(PFCandidateRefVector theInitialPFCands,double ChargedHadrCand_tkminPt,int ChargedHadrCand_tkminPixelHitsn,int ChargedHadrCand_tkminTrackerHitsn,double ChargedHadrCand_tkmaxipt,double ChargedHadrCand_tkmaxChi2,double ChargedHadrCand_tktorefpointmaxDZ,Vertex pV, double refpoint_Z);
  PFCandidateRefVector filteredPFNeutrHadrCands(PFCandidateRefVector theInitialPFCands,double NeutrHadrCand_HcalclusminEt);
  PFCandidateRefVector filteredPFGammaCands(PFCandidateRefVector theInitialPFCands,double GammaCand_EcalclusminEt);
  math::XYZPoint propagTrackECALSurfContactPoint(const MagneticField*,TrackRef); 
}

#endif

