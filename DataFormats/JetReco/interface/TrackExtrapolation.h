#ifndef DataFormats_TrackReco_interface_TrackExtrapolation_h
#define DataFormats_TrackReco_interface_TrackExtrapolation_h
/* \class reco::TrackExtrapolation TrackExtrapolation.h DataFormats/TrackReco/interface/TrackExtrapolation.h
*
* This class represents the track state at several radii (specified by user in producer).
* It stores a TrackRef to the original track, as well as vectors of the positions and  momenta
* of the track at the various radii. 
*
* \author Salvatore Rappoccio, JHU
*
*
*/

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include <vector>

namespace reco {
  class TrackExtrapolation {
    typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<Double32_t> > Point;
    typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double32_t> > Vector;
  public:


    TrackExtrapolation() {}
    TrackExtrapolation( reco::TrackRef const & track,
			std::vector<Point> const & pos,
			std::vector<Vector> const & mom) :
    track_(track)
    {
      pos_.resize( pos.size() );
      copy( pos.begin(), pos.end(), pos_.begin() );
      mom_.resize( mom.size() );
      copy( mom.begin(), mom.end(), mom_.begin() );
    }

    ~TrackExtrapolation() {}

    reco::TrackRef const &                 track()      const { return track_;}
    std::vector<Point> const &  positions()  const { return pos_;}
    std::vector<Vector> const & momenta()    const { return mom_;}

  protected:
    reco::TrackRef    track_;
    std::vector<Point>      pos_;
    std::vector<Vector>     mom_;
  };
}


#endif
