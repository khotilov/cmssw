#ifndef DataFormats_ParticleFlowReco_PFRecTrack_h
#define DataFormats_ParticleFlowReco_PFRecTrack_h

#include "DataFormats/ParticleFlowReco/interface/PFTrack.h"
#include <iostream>

namespace reco {

  /**\class PFRecTrack
     \brief reconstructed track used as an input to particle flow    

     Additional information w/r to PFTrack: 
     - algorithm used to reconstruct the track
     - track ID, soon to be replaced by a RefToBase to the corresponding Track

     \author Renaud Bruneliere, Michele Pioppi
     \date   July 2006
  */
  class PFRecTrack : public PFTrack {

  public:
    
    /// different types of fitting algorithms
    enum AlgoType_t {
      Unknown = 0,
      KF = 1, // Kalman filter 
      GSF = 2,
      KF_ELCAND=3// Gaussian sum filter
    };

    PFRecTrack();
  
    PFRecTrack(double charge, AlgoType_t algoType, int trackId);

    PFRecTrack(double charge, AlgoType_t algoType);

    PFRecTrack(const PFRecTrack& other);

    /// get type of algorithm
    unsigned int algoType() const { return algoType_; }

    friend  std::ostream& operator<<(std::ostream& out, 
				     const PFRecTrack& track);
    int recTrackId() const {return trackId_;}

  private:

    /// type of fitting algorithm used to reconstruct the track
    AlgoType_t algoType_;
    
    int trackId_;
  };

}

#endif
