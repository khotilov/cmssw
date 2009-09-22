#ifndef RecoEgamma_EgammaTools_ConversionFinder_h
#define RecoEgamma_EgammaTools_ConversionFinder_h

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/VectorUtil.h"

class ConversionFinder {
 public:
  ConversionFinder();
  ~ConversionFinder();
  //bField has to be supplied in Tesla
  static bool isElFromConversion(const reco::GsfElectron& gsfElectron, 
				 const edm::View<reco::Track>& v_tracks, 
				 const float bFieldAtOrigin,
				 const float maxAbsDist = 0.02,
				 const float maxAbsDCot = 0.02,
				 const float minFracSharedHits = 0.45);
  
  
  static std::pair<double, double> getConversionInfo(math::XYZTLorentzVector trk1_p4, 
						     int trk1_q, float trk1_d0, 
						     math::XYZTLorentzVector trk2_p4,
						     int trk2_q, float trk2_d0,
						     float bFieldAtOrigin);
  
};
#endif
