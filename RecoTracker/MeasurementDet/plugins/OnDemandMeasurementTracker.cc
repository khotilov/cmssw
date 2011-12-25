#include "OnDemandMeasurementTracker.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/Common/interface/ContainerMask.h"

#include "TrackingTools/MeasurementDet/interface/MeasurementDetException.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/Records/interface/TrackerCPERecord.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitMatcher.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPE.h"  

#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "TkStripMeasurementDet.h"
#include "TkPixelMeasurementDet.h"
#include "TkGluedMeasurementDet.h"

#include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
#include "CondFormats/DataRecord/interface/SiStripNoisesRcd.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibTracker/Records/interface/SiStripRegionCablingRcd.h"

#include <iostream>
#include <typeinfo>
#include <map>

#include <DataFormats/GeometrySurface/interface/BoundPlane.h>
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/UpdaterService.h"

using namespace std;

OnDemandMeasurementTracker::OnDemandMeasurementTracker(const edm::ParameterSet&              conf,
						       const PixelClusterParameterEstimator* pixelCPE,
						       const StripClusterParameterEstimator* stripCPE,
						       const SiStripRecHitMatcher*  hitMatcher,
						       const TrackerGeometry*  trackerGeom,
						       const GeometricSearchTracker* geometricSearchTracker,
						       const SiStripQuality *stripQuality,
                                                       int   stripQualityFlags,
                                                       int   stripQualityDebugFlags,
                                                       const SiPixelQuality *pixelQuality,
                                                       const SiPixelFedCabling *pixelCabling,
                                                       int   pixelQualityFlags,
                                                       int   pixelQualityDebugFlags,
						       const SiStripRegionCabling * stripRegionCabling,
						       bool isRegional):
  MeasurementTrackerImpl(conf,pixelCPE,stripCPE,hitMatcher,trackerGeom,geometricSearchTracker,
        stripQuality,stripQualityFlags,stripQualityDebugFlags,
        pixelQuality,pixelCabling,pixelQualityFlags,pixelQualityDebugFlags,
        isRegional)
  , category_("OnDemandMeasurementTracker")
  , StayPacked_(true)
  , StripOnDemand_(true)
  , PixelOnDemand_(false)
  , theStripRegionCabling(stripRegionCabling)
{
  //  the constructor does construct the regular MeasurementTracker
  //  then a smart copy of the DetMap is made into DetODMap: this could be avoided with modification to MeasurementDet interface
  //  the elementIndex to be defined in the refgetter is mapped to the detId
  //  flags are set to initialize the DetODMap
  
  std::map<SiStripRegionCabling::ElementIndex, std::vector< DetODContainer::const_iterator> > local_mapping;

  for (DetContainer::iterator it=theDetMap.begin(); it!= theDetMap.end();++it)
    {
      DetODContainer::iterator inserted = theDetODMap.insert(make_pair(it->first,DetODStatus(const_cast<MeasurementDet*>(it->second)))).first;


      GeomDet::SubDetector subdet = it->second->geomDet().subDetector();
      if (subdet == GeomDetEnumerators::PixelBarrel  || subdet == GeomDetEnumerators::PixelEndcap ){
	//special flag treatement for pixels
	//one can never be updated and not defined. except if we don't need to care about them: pixels
	inserted->second.defined=false;
	inserted->second.updated=true;
      }//pixel module
      else if (subdet == GeomDetEnumerators::TIB || subdet == GeomDetEnumerators::TOB ||
	  subdet == GeomDetEnumerators::TID || subdet == GeomDetEnumerators::TEC )
	{
	  //set flag to false
	  inserted->second.defined=false;
	  inserted->second.updated=false;

	  //what will be the element index in the refgetter
	  GlobalPoint center = it->second->geomDet().position();
	  double eta = center.eta();
	  double phi = center.phi();
	  uint32_t id = it->first;
	  SiStripRegionCabling::ElementIndex eIndex = theStripRegionCabling->elementIndex(SiStripRegionCabling::Position(eta,phi),
											  SiStripRegionCabling::subdetFromDetId(id),
											  SiStripRegionCabling::layerFromDetId(id));
	  LogDebug(category_)<<"region selected (from "<<id<<" center) is:\n"
			     <<"position: "<<center
			     <<"\n center absolute index: "<<theStripRegionCabling->region(theStripRegionCabling->positionIndex(SiStripRegionCabling::Position(eta,phi)))
			     <<"\n center position index: "<<theStripRegionCabling->positionIndex(SiStripRegionCabling::Position(eta,phi)).first<<
	    " "<<theStripRegionCabling->positionIndex(SiStripRegionCabling::Position(eta,phi)).second
			     <<"\n center postion: "<<theStripRegionCabling->position(theStripRegionCabling->positionIndex(SiStripRegionCabling::Position(eta,phi))).first<<
	    " "<<theStripRegionCabling->position(theStripRegionCabling->positionIndex(SiStripRegionCabling::Position(eta,phi))).second
			     <<"\n eta: "<<eta
			     <<"\n phi: "<<phi
			     <<"\n subedet: "<<SiStripRegionCabling::subdetFromDetId(id)
			     <<" layer: "<<SiStripRegionCabling::layerFromDetId(id);

	  //	  register those in a map
	  //to be able to know what are the detid in a given elementIndex
	  local_mapping[eIndex].push_back(inserted);
	}//strip module
      else{
	//abort
	edm::LogError(category_)<<"not a tracker geomdet in constructor: "<<it->first;
	throw MeasurementDetException("OnDemandMeasurementTracker dealing with a non tracker GeomDet.");
      }//abort
    }//loop over DetMap
  if (theInactiveStripDetectorLabels.size()!=0)
    theRawInactiveStripDetIds.reserve(200);

  //move into a vector
  region_mapping.reserve(local_mapping.size());
  for( auto eIt= local_mapping.begin();
       eIt!=local_mapping.end();++eIt)
    region_mapping.push_back(std::make_pair((*eIt).first,(*eIt).second));
}


void OnDemandMeasurementTracker::define( const edm::Handle< LazyGetter> & theLazyGetterH,
					 std::auto_ptr< RefGetter > &  theGetter ) const
{
  //  define is supposed to be call by an EDProducer module, which wil put the RefGetter in the event
  //  so that reference can be made to it.
  //  the lazy getter is retrieved by the calling module and passed along with the event
  //  the map is cleared, except for pixel
  //  then the known elementIndex are defined to the RefGetter. no unpacking is done at this time
  //  the defined region range is registered in the DetODMap for further use.

  //clear all defined tags to start from scratch (except for pixel)
  for (DetODContainer::iterator it=theDetODMap.begin(); it!= theDetODMap.end();++it)
    {
      if (it->second.updated && !it->second.defined) continue; //special treatement for pixels
      it->second.defined= false; 
      it->second.updated = false;
    }

  //define all the elementindex in the refgetter
  for(auto eIt= region_mapping.begin();
       eIt!=region_mapping.end();++eIt){
    std::pair<unsigned int, unsigned int> region_range; 
    
    //before update of the refgetter
    region_range.first = theGetter->size();
    //update the refegetter with the elementindex
    theStripRegionCabling->updateSiStripRefGetter<SiStripCluster> (*theGetter, theLazyGetterH, eIt->first);
    //after update of the refgetter
    region_range.second = theGetter->size();

    LogDebug(category_)<<"between index: "<<region_range.first<<" "<<region_range.second
		       <<"\n"<<dumpRegion(region_range,*theGetter,StayPacked_);
    
    //now assign to each measurement det for that element index
    for (auto dIt=eIt->second.begin();
	 dIt!=eIt->second.end();++dIt){
      DetODStatus & elem = const_cast<DetODStatus &>((*dIt)->second);
      elem.region_range = region_range;
      elem.defined=true;
      LogDebug(category_)<<"detId: "<<(*dIt)->first<<" in region range: "<<region_range.first<<" "<<region_range.second;
    }//loop over MeasurementDet attached to that elementIndex
  }//loop over know elementindex
}

void OnDemandMeasurementTracker::updateStrips( const edm::Event& event) const 
{
  bool oncePerEvent= edm::Service<UpdaterService>()->checkOnce("OnDemandMeasurementTracker::updateStrips::"+name_);
  bool failedToGet = false;
  if (!oncePerEvent)
    failedToGet = theRefGetterH.failedToGet() || theLazyGetterH.failedToGet();

  if (oncePerEvent || failedToGet)
    {
      LogDebug(category_)<<"Updating siStrip on event: "<< (unsigned int) event.id().run() <<" : "<<(unsigned int) event.id().event();
      
      //get the ref getter back from the event
      std::string stripClusterProducer = pset_.getParameter<std::string>("stripClusterProducer");
      event.getByLabel(stripClusterProducer,theRefGetterH);
      
      std::string stripLazyGetter = pset_.getParameter<std::string>("stripLazyGetterProducer");
      event.getByLabel(stripLazyGetter,theLazyGetterH);

      //get the skip clusters
      if (selfUpdateSkipClusters_){
        theSkipClusterRefs=true;
        event.getByLabel(pset_.getParameter<edm::InputTag>("skipClusters"),theStripClusterMask);
        theStripClusterMask->copyMaskTo(theStripsToSkip);
      } else {
        theStripsToSkip.clear();
      }

      //get the detid that are inactive
      theRawInactiveStripDetIds.clear();
      getInactiveStrips(event,theRawInactiveStripDetIds);
    }
}

void OnDemandMeasurementTracker::update( const edm::Event& event) const
{
  //  update is supposed to be called by any module that is useing the MeasurementTracker
  //  after checking the the event has not yet been seen
  //  update the pixel using MeasurementTracekr specific function
  //  retreive the RefGetter from the event: the very one that has been pass to define(...) and put into the event

  if (!PixelOnDemand_) {
    LogDebug(category_)<<"pixel are not OnDemand. updating them a la MeasurmentTracker.";
    MeasurementTrackerImpl::updatePixels(event);}
  else{
    edm::LogError(category_)<<"trying to update siPixel as on-demand. Not Implemented yet.";
  }

  if (!StripOnDemand_) {
    LogDebug(category_)<<"strip are not OnDemand. updating them a la MeasurmentTracker.";
    MeasurementTrackerImpl::updateStrips(event);}
  else{
    LogDebug(category_)<<"strip are OnDemand. updating them a la OnDemandMeasurmentTracker."; 
    updateStrips(event);
  }
}

#include <sstream>

std::string OnDemandMeasurementTracker::dumpCluster(const std::vector<SiStripCluster> ::const_iterator & begin,const  std::vector<SiStripCluster> ::const_iterator & end)const
{
  //  dumpCluster is a printout of all the clusters between the iterator. returns a string
  std::string tab="      ";
  std::stringstream ss;
  std::vector<SiStripCluster> ::const_iterator it = begin;
  unsigned int i=0;
  for (;it!=end;++it){
    ss<<tab<<i++<<") center: "<<it->barycenter()<<",id: "<<it->geographicalId()<<" with: "<<it->amplitudes().size()<<" strips\n"<<tab<<tab<<"{";
    for (unsigned int is=0;is!=it->amplitudes().size();++is){
      ss<<it->amplitudes()[is]<<" ";
    }ss<<"}\n";
  }
  return ss.str();
}

std::string OnDemandMeasurementTracker::dumpRegion(std::pair<unsigned int,unsigned int> indexes,
					     const RefGetter & theGetter,
					     bool stayPacked)const
{
  //  dumpRegion is a printout of all the clusters in a region defined on the RefGetter. returns a string
  std::stringstream ss;
  ss<<"cluster between: "<<indexes.first<<" and: "<<indexes.second<<"\n";
  for (unsigned int iRegion = indexes.first; iRegion != indexes.second; ++iRegion){    
    uint32_t reg = SiStripRegionCabling::region((theGetter)[iRegion].region());
    SiStripRegionCabling::Position pos = theStripRegionCabling->position(reg);
    SiStripRegionCabling::PositionIndex posI = theStripRegionCabling->positionIndex(reg);
    
    ss<<"Clusters for region:["<<iRegion<<"]"
      <<"\n element index: "<<(theGetter)[iRegion].region()
      <<"\n region absolute index: "<<reg
      <<"\n region position index: "<<posI.first<<" "<<posI.second
      <<"\n region position: "<<pos.first<<" "<<pos.second
      <<"\n"<< (stayPacked? " hidden to avoid unpacking." : dumpCluster((theGetter)[iRegion].begin(),(theGetter)[iRegion].end()));
  }
  return ss.str();
}

void OnDemandMeasurementTracker::assign(const TkStripMeasurementDet * csmdet,
				  DetODContainer::iterator * alreadyFound)const {
  //  assign is using the handle to the refgetter and the region index range to update the MeasurementDet with their clusters
  
  TkStripMeasurementDet * smdet = const_cast<TkStripMeasurementDet *>(csmdet);
  DetId id = smdet->geomDet().geographicalId();
  
  LogDebug(category_)<<"assigning: "<<id.rawId();

  // what is the iterator. do not look again if already found and provided with
  DetODContainer::iterator elementInMap;
  if (alreadyFound){ elementInMap=*alreadyFound;}
  else{ elementInMap = theDetODMap.find(id);}
  
  if  ( elementInMap != theDetODMap.end()){
    //flag it as updated
    elementInMap->second.updated = true;

    if (!theRawInactiveStripDetIds.empty() && std::binary_search(theRawInactiveStripDetIds.begin(), theRawInactiveStripDetIds.end(), id)) {
      smdet->setActiveThisEvent(false); 
      return;
    }

    //retrieve the region range index for this module
    std::pair<unsigned int,unsigned int> & indexes =elementInMap->second.region_range;

    //this printout will trigger unpacking. no problem. it is done on the next regular line (find(id.rawId())
    LogDebug(category_)<<"between index: "<<indexes.first<<" and: "<<indexes.second
		       <<"\nretrieved for module: "<<id.rawId()
		       <<"\n"<<dumpRegion(indexes,*theRefGetterH);
    
    //look for iterator range in the regions defined for that module
    for (unsigned int iRegion = indexes.first; iRegion != indexes.second; ++iRegion){
      RefGetter::record_pair range = (*theRefGetterH)[iRegion].find(id.rawId());
      if (range.first!=range.second){
	//	found something not empty
	//update the measurementDet
	smdet->update(range.first, range.second, theLazyGetterH, id);
	LogDebug(category_)<<"Valid clusters for: "<<id.rawId()
			   <<"\nnumber of regions defined here: "<< indexes.second-indexes.first
			   <<"\n"<<dumpCluster(range.first,range.second);
	/* since theStripsToSkip is a "static" pointer of the MT, no need to set it at all time.
	  if (selfUpdateSkipClusters_){
	  //assign skip clusters
	  smdet->setClusterToSkip(&theStripsToSkip);
	}
	*/
	//and you are done
	return;}
    }//loop over regions, between indexes

    //if reached. no cluster are found. set the TkStripMeasurementDet to be empty
    smdet->setEmpty();

  }//found in the map
  else{
    //throw excpetion
    edm::LogError(category_)<<"failed to find the MeasurementDet for: "<<id.rawId();
    throw MeasurementDetException("failed to find the MeasurementDet for: <see message logger>");
  }
}



const MeasurementDet* 
OnDemandMeasurementTracker::idToDet(const DetId& id) const
{
  //  overloaded from MeasurementTracker
  //  find the detid. if not found throw exception
  //  if already updated: always for pixel or strip already queried. return it
  DetODContainer::iterator it = theDetODMap.find(id);
  if ( it != theDetODMap.end()) {
    
    //it has already been queried. and so already updated: nothing to be done here (valid for pixels too)
    if (it->second.updated){LogDebug(category_)<<"found id: "<<id.rawId()<<" as aleardy updated."; return it->second.mdet;}
    
    if (StripOnDemand_){
      //check glued or single
      std::vector< const GeomDet*> comp = it->second.mdet->geomDet().components();
      if (!comp.empty()){
	//glued det
	LogDebug(category_)<<"updating glued id: "<<id.rawId()<<" ("<<comp.size()<<").";
	//cast to specific type
	TkGluedMeasurementDet*  theConcreteDet = static_cast<TkGluedMeasurementDet*>(it->second.mdet);
		
	//	update the two components
	//update the mono
	assign(theConcreteDet->monoDet());
	//update the stereo
	assign(theConcreteDet->stereoDet());
	      
	//flag the glued det as updated (components are flagged in assign)
	it->second.updated=true;
      }//glued det
      else{
	//single layer 
	LogDebug(category_)<<"updating singel id: "<<id.rawId();
	//cast to specific type
	TkStripMeasurementDet*  theConcreteDet = static_cast<TkStripMeasurementDet*>(it->second.mdet);

	//update the single layer
	assign(theConcreteDet,&it);
      }//single det
    }
    //eventually return it
    return it->second.mdet;
  }//found in the DetODMap
  else{
    //throw excpetion
    edm::LogError(category_)<<"failed to find the MeasurementDet for: "<<id.rawId();
    throw MeasurementDetException("failed to find the MeasurementDet for: <see message logger>");
  }
  return 0;
}


