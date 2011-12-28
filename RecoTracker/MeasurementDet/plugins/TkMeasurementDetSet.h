#ifndef TMeasurementDetSet_H
#define TkMeasurementDetSet_H

#include<vector>
class TkGluedMeasurementDet;
class SiStripRecHitMatcher;
class StripClusterParameterEstimator;
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/* Struct of arrays supporting "members of Tk...MeasurementDet
 * implemented with vectors, to be optimized...
 */


class TkMeasurementDetSet {
public:
  struct BadStripCuts {
     BadStripCuts() : maxBad(9999), maxConsecutiveBad(9999) {}
     BadStripCuts(const edm::ParameterSet &pset) :
        maxBad(pset.getParameter<uint32_t>("maxBad")),
        maxConsecutiveBad(pset.getParameter<uint32_t>("maxConsecutiveBad")) {}
     uint16_t maxBad, maxConsecutiveBad;
  };

  struct BadStripBlock {
      short first;
      short last;
      BadStripBlock(const SiStripBadStrip::data &data) : first(data.firstStrip), last(data.firstStrip+data.range-1) { }
  };


  TkMeasurementDetSet(const SiStripRecHitMatcher* matcher,
		      const StripClusterParameterEstimator* cpe,
		      bool regional):
    theMatcher(matcher), theCPE(cpe), skipClusters_(0), regional_(regional){}
  

  void init() {
    // assume vector is full and ordered!
    int size = theStripDets.size();

    empty_.resize(size,true);
    activeThisEvent_.resize(size,true);
    activeThisPeriod_.resize(size,true);
    id_.resize(size);
    subId_.resize(size);
    totalStrips_.resize[size];

    bad128Strip_.resize(size*6);
    hasAny128StripBad_.resize(size);
    maskBad128StripBlocks_.resize(size);

    if (isRegional()) {
      detset_.resize(size);
    }  else {
    beginClusterI_.resize(size);
    endClusterI_.resize(size);
    }

    for (int i=0; i!=size; ++i) {
      auto & mdet =  *theStripDets[i]; 
      mdet.setIndex(i);
      //intialize the detId !
      id_[i] = mdet->gdet->geographicalId().rawId();
      subId_[i]=SiStripDetId(id_[i]).subdetId()-3;
      //initalize the total number of strips
      totalStrips_[i] =  mdet->specificGeomDet().specificTopology().nstrips();
    }
  }

  const std::vector<TkStripMeasurementDet*>& stripDets() const {return  theStripDets;}


  std::vector<bool> const & clusterToSkip() const { return theStripsToSkip; }
 

  void update(int i,
	      const detset &detSet ) { 
    detSet_[i] = detSet; 

    empty_[i] = false;
  }

  void update(int i,
	      std::vector<SiStripCluster>::const_iterator begin ,std::vector<SiStripCluster>::const_iterator end) { 
    beginClusterI_[i] = begin - regionalHandle_->begin_record();
    endClusterI_[i] = end - regionalHandle_->begin_record();

    empty_[i] = false;
    activeThisEvent_[i] = true;
  }

  edm::Handle<edmNew::DetSetVector<SiStripCluster> > & handle() {  return handle_;}
  edm::Handle<edm::LazyGetter<SiStripCluster> > & regionalHandle() { return regionalHandle_;}



  bool isRegional() const { return regional_;}
  bool empty(int i) const { return empty_[i];}  
  bool isActive(int i) const { return activeThisEvent_[i] && activeThisPeriod_[i]; }
  void setEmpty(int i) {empty_[i] = true; activeThisEvent_[i] = true; }

  void setEmpty() {
    std::fill(empty_.begin(),empty_.end(),true;);
    std::fill(activeThisEvent_.begin(),activeThisEvent_.end(),true;);
  }

  /** \brief Turn on/off the module for reconstruction, for the full run or lumi (using info from DB, usually).
             This also resets the 'setActiveThisEvent' to true */
  void setActive(int it, bool active) { activeThisPeriod_[i] = active; activeThisEvent_[i] = true; if (!active) empty_[i] = true; }
  /** \brief Turn on/off the module for reconstruction for one events.
             This per-event flag is cleared by any call to 'update' or 'setEmpty'  */
  void setActiveThisEvent(int i, bool active) { activeThisEvent_[i] = active;  if (!active) empty[i] = true; }



  unsigned int beginClusterI(int i) const {return beginClusterI[i];}
  unsigned int endClusterI(int i) const {return endClusterI[i];}

  int totalStrips(int) const { return totalStrips_[i];}

  BadStripCuts & badStripCuts(int i) { return  badStripCuts(subId_[i]);}

  bool hasAny128StripBad(int i) const { return  hasAny128StripBad_[i];}

  bool isMasked(int i, const SiStripCluster &cluster) const {
    int offset =  nbad128*i;
    if ( bad128Strip_[offset+cluster.firstStrip() >> 7] ) {
      if ( bad128Strip_[offset+(cluster.firstStrip()+cluster.amplitudes().size())  >> 7] ||
	   bad128Strip_[offset+static_cast<int32_t>(cluster.barycenter()-0.499999) >> 7] ) {
	return true;
      }
    } else {
      if ( bad128Strip_[offset+(cluster.firstStrip()+cluster.amplitudes().size())  >> 7] &&
	   bad128Strip_[offset+static_cast<int32_t>(cluster.barycenter()-0.499999) >> 7] ) {
	return true;
      }
    }
    return false;
  }
 

  void set128StripStatus(int i, bool good, int idx) { 
    int offset =  nbad128*i;
    if (idx == -1) {
      std::fill(bad128Strip_[offset], bad128Strip_[offset+6], !good);
      hasAny128StripBad_[i] = !good;
    } else {
      bad128Strip_[offset+idx] = !good;
      if (good == false) {
	hasAny128StripBad_[i] = false;
      } else { // this should not happen, as usually you turn on all fibers
	// and then turn off the bad ones, and not vice-versa,
	// so I don't care if it's not optimized
	hasAny128StripBad_[i] = true;
	for (int j = 0; i < (totalStrips_[j] >> 7); j++) {
	  if (bad128Strip_[j+offset] == false) hasAny128StripBad_[i] = false;
	}
      }    
    } 
  }

private:

  friend class  MeasurementTrackerImpl;
  mutable vector<TkStripMeasurementDet*> theStripDets;
  
  // globals
  const SiStripRecHitMatcher*       theMatcher;
  const StripClusterParameterEstimator* theCPE;

  edm::Handle<edmNew::DetSetVector<SiStripCluster> > handle_;
  edm::Handle<edm::LazyGetter<SiStripCluster> > regionalHandle_;

  mutable std::vector<bool> theStripsToSkip;

  bool regional_;

  BadStripCuts badStripCuts[4];

  // members of TkStripMeasurementDet
  std::vector<unsigned int> id_;
  std::vector<unsigned char> subId_;

  std::vector<int> totalStrips_;

  const int nbad128 = 6;
  std::vector<bool> bad128Strip_;
  std::vector<bool> hasAny128StripBad_, maskBad128StripBlocks_;

  std::vector<bool> empty_;

  std::vector<bool> activeThisEvent_,activeThisPeriod_;

  // full reco
  std::vector<detset> detSet_;

  // --- regional unpacking

  std::vector<unsigned int> beginClusterI_;
  std::vector<unsigned int> endClusterI_;


};


#endif
