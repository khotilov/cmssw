#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"



namespace {
  void
  throwExceptionUninitialized(const char *where)
  {
    throw cms::Exception("SiStripRecHit1D") << 
      "Trying to access " << where << " for a RecHit that was read from disk, but since CMSSW_2_1_X local positions are transient.\n" <<
      "If you want to get coarse position/error estimation from disk, please set: ComputeCoarseLocalPositionFromDisk = True \n " <<
      " to the TransientTrackingRecHitBuilder you are using from RecoTracker/TransientTrackingRecHit/python/TTRHBuilders_cff.py";
  }
  
  void obsolete() {
    throw cms::Exception("SiStripRecHit1D") << "CLHEP is obsolete for Tracker Hits";
  }
}


SiStripRecHit1D::SiStripRecHit1D( const LocalPoint& pos, const LocalError& err,
				  const DetId& id,
				  ClusterRef const & cluster):
  Base(id),
  sigmaPitch_(-1.),pos_(pos),err_(err),
  cluster_(cluster){}


SiStripRecHit1D::SiStripRecHit1D( const LocalPoint& pos, const LocalError& err,
				  const DetId& id,
				  ClusterRegionalRef const& cluster): 
  Base(id),
  sigmaPitch_(-1.),pos_(pos),err_(err),
  cluster_(cluster){}


SiStripRecHit1D::SiStripRecHit1D(const SiStripRecHit2D* hit2D):
  Base(hit2D->geographicalId()),
  sigmaPitch_(-1),pos_(hit2D->localPosition()),
  err_(hit2D->localPositionError().xx(),0.,DBL_MAX),
  cluster_(hit2D->omniCluster())
{}


bool SiStripRecHit1D::hasPositionAndError() const {
    return (err_.xx() != 0) || (err_.yy() != 0) || (err_.xy() != 0) ||
           (pos_.x()  != 0) || (pos_.y()  != 0) || (pos_.z()  != 0);
}

LocalPoint SiStripRecHit1D::localPosition() const {
    if (!hasPositionAndError()) throwExceptionUninitialized("localPosition");
    return pos_;
}

LocalError SiStripRecHit1D::localPositionError() const{ 
    if (!hasPositionAndError()) throwExceptionUninitialized("localPositionError");
    return err_;
}

void 
SiStripRecHit1D::getKfComponents( KfComponentsHolder & holder ) const 
{
   if (!hasPositionAndError()) throwExceptionUninitialized("getKfComponents");
   AlgebraicVector1 & pars = holder.params<1>();
   pars[0] = pos_.x(); 

   AlgebraicSymMatrix11 & errs = holder.errors<1>();
   errs(0,0) = err_.xx();

   AlgebraicMatrix15 & proj = holder.projection<1>();
   proj(0,3) = 1;

   holder.measuredParams<1>() = AlgebraicVector1( holder.tsosLocalParameters().At(3) );
   holder.measuredErrors<1>() = holder.tsosLocalErrors().Sub<AlgebraicSymMatrix11>( 3, 3 );
}



bool 
SiStripRecHit1D::sharesInput( const TrackingRecHit* other, 
			      SharedInputType what) const
{
  //here we exclude non si-strip subdetectors
  if( ((geographicalId().rawId()) >> (DetId::kSubdetOffset) ) != ( (other->geographicalId().rawId())>> (DetId::kSubdetOffset)) ) return false;

  //Protection against invalid hits
  if(! other->isValid()) return false;

  const std::type_info & otherType = typeid(*other);
  if (otherType == typeid(SiStripRecHit2D)) {
    const SiStripRecHit2D* otherCast = static_cast<const SiStripRecHit2D*>(other);
    // as 'null == null' is true, we can't just "or" the two equality tests: one of the two refs is always null! (gpetrucc)
    if (cluster().isNonnull()) {
      return (cluster() == otherCast->cluster());
    } else {
      return (cluster_regional() == otherCast->cluster_regional());
    }
  } else if (otherType == typeid(SiStripRecHit1D)) {
    const SiStripRecHit1D* otherCast = static_cast<const SiStripRecHit1D*>(other);
    // as 'null == null' is true, we can't just "or" the two equality tests: one of the two refs is always null! (gpetrucc)
    if (cluster().isNonnull()) {
      return (cluster() == otherCast->cluster());
    } else {
      return (cluster_regional() == otherCast->cluster_regional());
    }
  } else if (otherType == typeid(ProjectedSiStripRecHit2D)) {
    const SiStripRecHit2D* otherCast = & (static_cast<const ProjectedSiStripRecHit2D*>(other)->originalHit());
    // as 'null == null' is true, we can't just "or" the two equality tests: one of the two refs is always null! (gpetrucc)
    if (cluster().isNonnull()) {
      return (cluster() == otherCast->cluster());
    } else {
      return (cluster_regional() == otherCast->cluster_regional());
    }
  } else if ((otherType == typeid(SiStripMatchedRecHit2D)) && (what == all)) {
    return false; 
  } else {
    // last resort, recur to 'recHits()', even if it returns a vector by value
    std::vector<const TrackingRecHit*> otherHits = other->recHits();
    int ncomponents=otherHits.size();
    if(ncomponents==0)return false;
    else if(ncomponents==1)return sharesInput(otherHits.front(),what);
    else if (ncomponents>1){
      if(what == all )return false;
      else{
	for(int i=0;i<ncomponents;i++){
	  if(sharesInput(otherHits[i],what))return true;
	}
	return false;
      }
    }
    return false;
  }
}


std::vector<const TrackingRecHit*> SiStripRecHit1D::recHits() const {
  std::vector<const TrackingRecHit*> nullvector;
  return nullvector; 
}
std::vector<TrackingRecHit*> SiStripRecHit1D::recHits() {
  std::vector<TrackingRecHit*> nullvector;
  return nullvector; 
}

 // obsolete (for what tracker is concerned...) interface
AlgebraicVector SiStripRecHit1D::parameters() const {
  obsolete();
  return AlgebraicVector();
}

AlgebraicSymMatrix SiStripRecHit1D::parametersError() const {
  obsolete();
  return AlgebraicSymMatrix();
}


AlgebraicMatrix SiStripRecHit1D::projectionMatrix() const {
  obsolete();
  return AlgebraicMatrix();
}
