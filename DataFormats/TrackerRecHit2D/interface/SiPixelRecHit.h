#ifndef DataFormats_SiPixelRecHit_h
#define DataFormats_SiPixelRecHit_h 1

//---------------------------------------------------------------------------
//!  \class SiPixelRecHit
//!  \brief Pixel Reconstructed Hit
//!
//!  A pixel hit is a 2D position and error in a given
//!  pixel sensor. It contains a persistent reference edm::Ref to the
//!  pixel cluster. 
//!
//!  \author porting from ORCA: Petar Maksimovic (JHU), 
//!          DetSetVector and persistent references: V.Chiochia (Uni Zurich)
//---------------------------------------------------------------------------

//! Our base class
#include "DataFormats/TrackerRecHit2D/interface/BaseSiTrackerRecHit2DLocalPos.h"
//! Quality word packing
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitQuality.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Ref.h"



class SiPixelRecHit : public  BaseSiTrackerRecHit2DLocalPos {
public:

  typedef edm::Ref<edmNew::DetSetVector<SiPixelCluster>, SiPixelCluster > ClusterRef;

  SiPixelRecHit(): BaseSiTrackerRecHit2DLocalPos(), qualWord_(0), cluster_()  {}

  ~SiPixelRecHit() {}
  
  SiPixelRecHit( const LocalPoint&, const LocalError&,
		 const DetId&, 
		 ClusterRef const&  cluster);  

  virtual SiPixelRecHit * clone() const {return new SiPixelRecHit( * this); }
  
  ClusterRef const& cluster() const { return cluster_;}
  void setClusterRef(const ClusterRef &ref) { cluster_  = ref; }

  virtual bool sharesInput( const TrackingRecHit* other, SharedInputType what) const;


  //--------------------------------------------------------------------------
  //--- Accessors of other auxiliary quantities
  //--- Added Oct 07 by Petar for 18x.
private:
  // *************************************************************************
  //
  SiPixelRecHitQuality::QualWordType  qualWord_ ;   // unsigned int 32-bit wide
  //
  // *************************************************************************

public:
  //--- The overall probability.  flags is the 32-bit-packed set of flags that
  //--- our own concrete implementation of clusterProbability() uses to direct
  //--- the computation based on the information stored in the quality word
  //--- (and which was computed by the CPE).  The default of flags==0 returns
  //--- probabilityY() only (as that's the safest thing to do).
  //--- Flags are static and kept in the transient rec hit.
  float clusterProbability(unsigned int flags = 0) const;


  //--- Allow direct access to the packed quality information.
  inline SiPixelRecHitQuality::QualWordType rawQualityWord() const { 
    return qualWord_ ; 
  }
  inline void setRawQualityWord( SiPixelRecHitQuality::QualWordType w ) { 
    qualWord_ = w; 
  }


  //--- Template fit probability, in X and Y directions
  inline float probabilityX() const     {
    return SiPixelRecHitQuality::thePacking.probabilityX( qualWord_ );
  }
  inline float probabilityY() const     {
    return SiPixelRecHitQuality::thePacking.probabilityY( qualWord_ );
  }

	//--- Charge `bin' (values 0, 1, 2, 3) according to Morris's template
  //--- code. qBin==4 is unphysical, qBin=5,6,7 are yet unused)
  //
  inline int qBin() const     {
    return SiPixelRecHitQuality::thePacking.qBin( qualWord_ );
  }

  //--- Quality flags (true or false):

  //--- The cluster is on the edge of the module, or straddles a dead ROC
  inline bool isOnEdge() const     {
    return SiPixelRecHitQuality::thePacking.isOnEdge( qualWord_ );
  }
  //--- The cluster contains bad pixels, or straddles bad column or two-column
  inline bool hasBadPixels() const     {
    return SiPixelRecHitQuality::thePacking.hasBadPixels( qualWord_ );
  }
  //--- The cluster spans two ROCS (and thus contains big pixels)
  inline bool spansTwoROCs() const     {
    return SiPixelRecHitQuality::thePacking.spansTwoROCs( qualWord_ );
  }

	//--- Quality flag for whether the probability is filled
	inline bool hasFilledProb() const {
		return SiPixelRecHitQuality::thePacking.hasFilledProb( qualWord_ );
	}
  
  //--- Setters for the above
	inline void setProbabilityX( float prob ) {
    SiPixelRecHitQuality::thePacking.setProbabilityX( prob, qualWord_ );
  }
  inline void setProbabilityY( float prob ) {
    SiPixelRecHitQuality::thePacking.setProbabilityY( prob, qualWord_ );
  }  
  inline void setQBin( int qbin ) {
    SiPixelRecHitQuality::thePacking.setQBin( qbin, qualWord_ );
  }
  inline void setIsOnEdge( bool flag ) {
    SiPixelRecHitQuality::thePacking.setIsOnEdge( flag, qualWord_ );
  }
  inline void setHasBadPixels( bool flag ) {
    SiPixelRecHitQuality::thePacking.setHasBadPixels( flag, qualWord_ );
  }
  inline void setSpansTwoROCs( bool flag ) {
    SiPixelRecHitQuality::thePacking.setSpansTwoROCs( flag, qualWord_ );
  }
	inline void setHasFilledProb( bool flag ) {
		SiPixelRecHitQuality::thePacking.setHasFilledProb( flag, qualWord_ );
	}

private:

  SiPixelClusterRefNew cluster_;

};

// Comparison operators
inline bool operator<( const SiPixelRecHit& one, const SiPixelRecHit& other) {
  if ( one.geographicalId() < other.geographicalId() ) {
    return true;
  } else {
    return false;
  }
}

#endif
