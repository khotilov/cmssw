#ifndef RecoTracker_MeasurementDet_OnDemandMeasurementTracker_H
#define RecoTracker_MeasurementDet_OnDemandMeasurementTracker_H

#include "TrackingTools/MeasurementDet/interface/MeasurementDetSystem.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "DataFormats/SiStripCommon/interface/SiStripRefGetter.h" 
#include "CalibFormats/SiStripObjects/interface/SiStripRegionCabling.h"

 
class OnDemandMeasurementTracker : public MeasurementTracker {
public:
  /// constructor
  OnDemandMeasurementTracker(const edm::ParameterSet&              conf,
			     const PixelClusterParameterEstimator* pixelCPE,
			     const StripClusterParameterEstimator* stripCPE,
			     const SiStripRecHitMatcher*  hitMatcher,
			     const TrackerGeometry*  trackerGeom,
			     const GeometricSearchTracker* geometricSearchTracker,
			     const SiStripDetCabling *stripCabling,
			     const SiStripNoises *stripNoises,
			     const SiStripRegionCabling * stripRegionCabling,
			     bool  isRegional=false);
  /// destructor
  virtual ~OnDemandMeasurementTracker() {}
 
  /// MeasurementTracker overloaded function
  void update( const edm::Event&) const;

  typedef edm::SiStripLazyGetter<SiStripCluster> LazyGetter;
  typedef edm::SiStripRefGetter<SiStripCluster> RefGetter;

  /// OnDemandMeasurementTracker specific function to be called to define the region in the RefGetter according to MeasurementDet content
  void define(const edm::Handle< edm::SiStripLazyGetter<SiStripCluster> > & ,
	      std::auto_ptr< RefGetter > &  ) const;

  /// MeasurementDetSystem interface
  virtual const MeasurementDet*       idToDet(const DetId& id) const;
    
 private:
  /// log category
  std::string category_;
  /// internal flag to avoid unpacking things with LogDebug on
  bool StayPacked_;

  /// internal flag to do strip on demand (not configurable) true by default
  bool StripOnDemand_;
  /// internal flag to do pixel on demand (not configurable) false by default
  bool PixelOnDemand_;
  
  /// the cabling region tool to update a RefGetter
  const  SiStripRegionCabling * theStripRegionCabling;
  
  /// the handle is retrieved from the event to make reference to cluster in it
  mutable edm::Handle< edm::SiStripRefGetter<SiStripCluster> > theGetterH;

  /// a class that holds flags, region_range (in RefGetter) for a given MeasurementDet
  class DetODStatus {
  public:
    DetODStatus(MeasurementDet * m):defined(false),updated(false),mdet(m){ region_range = std::pair<uint,uint>(0,0);}
      bool defined;
      bool updated;
      std::pair<uint, uint> region_range;
      MeasurementDet * mdet;
  };

  typedef std::map<DetId, DetODStatus> DetODContainer;
  /// mapping of detid -> MeasurementDet+flags+region_range
  mutable DetODContainer theDetODMap;

  /// mapping of elementIndex -> iterator to the DetODMap: to know what are the regions that needs to be defined in the ref getter
  mutable std::map<SiStripRegionCabling::ElementIndex, std::vector< DetODContainer::iterator> > region_mapping;

  /// assigne the cluster iterator to the TkStipMeasurementDet (const_cast in the way)
    void assign(const  TkStripMeasurementDet * csmdet,
	      DetODContainer::iterator * alreadyFound=0) const;

  /// some printouts, exclusively under LogDebug
  std::string dumpCluster(const std::vector<SiStripCluster> ::const_iterator & begin, const  std::vector<SiStripCluster> ::const_iterator& end)const;
  std::string dumpRegion(std::pair<uint,uint> indexes,
			 const RefGetter & theGetter,
			 bool stayUnpacked = false)const;

};

#endif
