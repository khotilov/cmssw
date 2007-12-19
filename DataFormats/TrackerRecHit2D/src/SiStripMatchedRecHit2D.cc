#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"


SiStripMatchedRecHit2D::SiStripMatchedRecHit2D( const LocalPoint& pos, const LocalError& err,
								const DetId& id , const SiStripRecHit2D* rMono,const SiStripRecHit2D* rStereo): BaseSiTrackerRecHit2DLocalPos(pos, err, id), componentMono_(*rMono),componentStereo_(*rStereo){}

bool 
SiStripMatchedRecHit2D::sharesInput( const TrackingRecHit* other, 
				     SharedInputType what) const
{
  if (geographicalId() != other->geographicalId()) return false;
  
  const SiStripMatchedRecHit2D* otherMatched = 
    dynamic_cast<const SiStripMatchedRecHit2D*>(other);
  if ( otherMatched == 0 )  return false;

  if ( what == all) {
    return (monoHit()->sharesInput( otherMatched->monoHit(),what) && 
	    stereoHit()->sharesInput( otherMatched->stereoHit(),what));
  }
  else {
    return (monoHit()->sharesInput( otherMatched->monoHit(),what) || 
	    stereoHit()->sharesInput( otherMatched->stereoHit(),what));
  }
}

std::vector<const TrackingRecHit*>
SiStripMatchedRecHit2D::recHits()const {
  std::vector<const TrackingRecHit*> rechits(2);
  rechits[0]=&componentMono_;
  rechits[1]=&componentStereo_;
  return rechits;
}

std::vector<TrackingRecHit*>
SiStripMatchedRecHit2D::recHits() {
  std::vector<TrackingRecHit*> rechits(2);
  rechits[0]=&componentMono_;
  rechits[1]=&componentStereo_;
  return rechits;
}


