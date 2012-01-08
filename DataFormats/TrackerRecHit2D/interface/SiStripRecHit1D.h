#ifndef SiStripRecHit1D_H
#define SiStripRecHit1D_H



#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

class SiStripRecHit1D : public TrackerSingleRecHit { 
public:

 
  SiStripRecHit1D(): sigmaPitch_(-1.){}
  
  
  typedef OmniClusterRef::ClusterRef         ClusterRef;
  typedef OmniClusterRef::ClusterRegionalRef ClusterRegionalRef;


  SiStripRecHit1D( const LocalPoint& p, const LocalError& e,
		   const DetId& id, 
		   ClusterRef const&  clus) : TrackerSingleRecHit(p,e,id,clus), sigmaPitch_(-1.){}

  SiStripRecHit1D( const LocalPoint& p, const LocalError& e,
		   const DetId& id, 
		   ClusterRegionalRef const& clus) : TrackerSingleRecHit(p,e,id,clus), sigmaPitch_(-1.){}
  
  /// method to facilitate the convesion from 2D to 1D hits
  SiStripRecHit1D(const SiStripRecHit2D*);

  virtual SiStripRecHit1D * clone() const {return new SiStripRecHit1D( * this); }
  

  virtual int dimension() const {return 1;}
  virtual void getKfComponents( KfComponentsHolder & holder ) const {getKfComponents1D(holder);}

 
  double sigmaPitch() const { return sigmaPitch_;}
  void setSigmaPitch(double sigmap) const { sigmaPitch_=sigmap;}

private:
 
 /// cache for the matcher....
  mutable double sigmaPitch_;  // transient.... 
};

#endif
