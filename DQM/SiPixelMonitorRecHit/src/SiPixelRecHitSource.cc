// -*- C++ -*-
//
// Package:    SiPixelMonitorRecHits
// Class:      SiPixelRecHitSource
// 
/**\class 

 Description: Pixel DQM source for RecHits

 Implementation:
     Originally based on the code for Digis, adapted
	to read RecHits and create relevant histograms
*/
//
// Original Author:  Vincenzo Chiochia
//         Created:  
// $Id: SiPixelRecHitSource.cc,v 1.2 2007/10/19 12:02:32 krose Exp $
//
//
// Adapted by:  Keith Rose
//  	For use in SiPixelMonitorClient for RecHits

#include "DQM/SiPixelMonitorRecHit/interface/SiPixelRecHitSource.h"
// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
// DQM Framework
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQM/SiPixelCommon/interface/SiPixelFolderOrganizer.h"
// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
// DataFormats
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"


#include "DataFormats/Common/interface/DetSetVector.h"
//
#include <boost/cstdint.hpp>
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;
using namespace edm;

SiPixelRecHitSource::SiPixelRecHitSource(const edm::ParameterSet& iConfig) :
  conf_(iConfig),
  src_( conf_.getParameter<edm::InputTag>( "src" ) )
{
   theDMBE = edm::Service<DaqMonitorBEInterface>().operator->();
   LogInfo ("PixelDQM") << "SiPixelRecHitSource::SiPixelRecHitSource: Got DQM BackEnd interface"<<endl;
}


SiPixelRecHitSource::~SiPixelRecHitSource()
{
   // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  LogInfo ("PixelDQM") << "SiPixelRecHitSource::~SiPixelRecHitSource: Destructor"<<endl;
}


void SiPixelRecHitSource::beginJob(const edm::EventSetup& iSetup){

  LogInfo ("PixelDQM") << " SiPixelRecHitSource::beginJob - Initialisation ... " << std::endl;
  eventNo = 0;
	
  // Build map
  buildStructure(iSetup);
	
  // Book Monitoring Elements
  bookMEs();

}


void SiPixelRecHitSource::endJob(void){
  cout << "here" << endl;
  cout << " SiPixelDigiSource::endJob - Saving Root File " << std::endl;
  std::string outputFile = conf_.getParameter<std::string>("outputFile");
  cout << "ending" << endl;
  theDMBE->save( outputFile );
  cout << "last" << endl;
}

//------------------------------------------------------------------
// Method called for every event
//------------------------------------------------------------------
void SiPixelRecHitSource::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventNo++;
  //cout << eventNo << endl;
  // get input data
  edm::Handle<SiPixelRecHitCollection>  recHitColl;
  iEvent.getByLabel( src_, recHitColl );

  std::map<uint32_t,SiPixelRecHitModule*>::iterator struct_iter;
  for (struct_iter = thePixelStructure.begin() ; struct_iter != thePixelStructure.end() ; struct_iter++) {
    uint32_t TheID = (*struct_iter).first;
	
    SiPixelRecHitCollection::range pixelrechitRange = (recHitColl.product())->get(TheID);
    SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorBegin = pixelrechitRange.first;
    
	SiPixelRecHitCollection::const_iterator pixelrechitRangeIteratorEnd = pixelrechitRange.second;
      SiPixelRecHitCollection::const_iterator pixeliter = pixelrechitRangeIteratorBegin;

      // if( pixelrechitRangeIteratorBegin == pixelrechitRangeIteratorEnd) {cout << "oops" << endl;}
      float rechit_x = 0;
      float rechit_y = 0;

  for ( ; pixeliter != pixelrechitRangeIteratorEnd; pixeliter++) 
	{
	  


	  //cout << TheID << endl;
	  edm::Ref<edm::DetSetVector<SiPixelCluster>, SiPixelCluster> const& clust = pixeliter->cluster();
	  int sizeX = (*clust).sizeX();
	  //cout << sizeX << endl;
	  int sizeY = (*clust).sizeY();
	  //cout << sizeY << endl;
	  LocalPoint lp = pixeliter->localPosition();
	  rechit_x = lp.x();
	  rechit_y = lp.y();
	  
	  LocalError lerr = pixeliter->localPositionError();
	  //  float lerr_x = sqrt(lerr.xx());
	  //  float lerr_y = sqrt(lerr.yy());
	  //cout << "hh" << endl;
	  (*struct_iter).second->fill(rechit_x, rechit_y, sizeX, sizeY);
	  //cout << "ii" << endl;
	
	}
    
  }

  // slow down...
  //usleep(100000);
  
}

//------------------------------------------------------------------
// Build data structure
//------------------------------------------------------------------
void SiPixelRecHitSource::buildStructure(const edm::EventSetup& iSetup){

  LogInfo ("PixelDQM") <<" SiPixelRecHitSource::buildStructure" ;
  edm::ESHandle<TrackerGeometry> pDD;
  iSetup.get<TrackerDigiGeometryRecord>().get( pDD );

  LogVerbatim ("PixelDQM") << " *** Geometry node for TrackerGeom is  "<<&(*pDD)<<std::endl;
  LogVerbatim ("PixelDQM") << " *** I have " << pDD->dets().size() <<" detectors"<<std::endl;
  LogVerbatim ("PixelDQM") << " *** I have " << pDD->detTypes().size() <<" types"<<std::endl;
  
  for(TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++){
    
    if(dynamic_cast<PixelGeomDetUnit*>((*it))!=0){

      DetId detId = (*it)->geographicalId();
      // const GeomDetUnit      * geoUnit = pDD->idToDetUnit( detId );
      //const PixelGeomDetUnit * pixDet  = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);

     	  
	  
	      // SiPixelRecHitModule *theModule = new SiPixelRecHitModule(id, rechit_x, rechit_y, x_res, y_res, x_pull, y_pull);
	
	
	      if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {
		LogDebug ("PixelDQM") << " ---> Adding Barrel Module " <<  detId.rawId() << endl;
		uint32_t id = detId();
	
	SiPixelRecHitModule* theModule = new SiPixelRecHitModule(id);
		thePixelStructure.insert(pair<uint32_t,SiPixelRecHitModule*> (id,theModule));
		
	      }	else if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
		LogDebug ("PixelDQM") << " ---> Adding Endcap Module " <<  detId.rawId() << endl;
		uint32_t id = detId();
		SiPixelRecHitModule* theModule = new SiPixelRecHitModule(id);
		thePixelStructure.insert(pair<uint32_t,SiPixelRecHitModule*> (id,theModule));
	      }
      
	}	    
  }

  LogInfo ("PixelDQM") << " *** Pixel Structure Size " << thePixelStructure.size() << endl;
}
//------------------------------------------------------------------
// Book MEs
//------------------------------------------------------------------
void SiPixelRecHitSource::bookMEs(){
  
  std::map<uint32_t,SiPixelRecHitModule*>::iterator struct_iter;
  theDMBE->setVerbose(0);
    
  SiPixelFolderOrganizer theSiPixelFolder;
  
  for(struct_iter = thePixelStructure.begin(); struct_iter != thePixelStructure.end(); struct_iter++){
    
    /// Create folder tree and book histograms 
    if(theSiPixelFolder.setModuleFolder((*struct_iter).first)){
      (*struct_iter).second->book( conf_ );
    } else {
      throw cms::Exception("LogicError")
	<< "[SiPixelDigiSource::bookMEs] Creation of DQM folder failed";
    }
    
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelRecHitSource);
