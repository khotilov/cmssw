// -*- C++ -*-
//
// Package:    SiPixelMonitorCluster
// Class:      SiPixelClusterSource
// 
/**\class 

 Description: Pixel DQM source for Clusters

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Vincenzo Chiochia
//                   
//         Created:  
// $Id: SiPixelClusterSource.cc,v 1.7 2007/02/02 08:53:26 chiochia Exp $
//
//
#include "DQM/SiPixelMonitorCluster/interface/SiPixelClusterSource.h"
// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
// DQM Framework
#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
// DataFormats
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
//#include "DataFormats/SiPixelCluster/interface/PixelClusterCollection.h"
//#include "DataFormats/SiPixelCluster/interface/PixelCluster.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVector.h"
//
#include <boost/cstdint.hpp>
#include <string>
#include <stdlib.h>
using namespace std;

SiPixelClusterSource::SiPixelClusterSource(const edm::ParameterSet& iConfig)
{

   //now do what ever initialization is needed
   theDMBE = edm::Service<DaqMonitorBEInterface>().operator->();
   cout<<endl<<"SiPixelClusterSource::SiPixelClusterSource: BackEnd interface"<<endl;
}


SiPixelClusterSource::~SiPixelClusterSource()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


void SiPixelClusterSource::beginJob(const edm::EventSetup& iSetup){

  std::cout << " SiPixelClusterSource::beginJob - Initialisation ..... " << std::endl;
  eventNo = 0;
  // Build map
  buildStructure(iSetup);
  // Book Monitoring Elements
  bookMEs();

}


void SiPixelClusterSource::endJob(void){
  std::cout << " SiPixelClusterSource::endJob - Saving Root File " << std::endl;
  theDMBE->save("sourceoutputfile.root");
}

//------------------------------------------------------------------
// Method called for every event
//------------------------------------------------------------------
void
SiPixelClusterSource::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  eventNo++;
  
  std::cout << " Processing event: " << eventNo << std::endl;

  // retrieve producer name of input SiPixelClusterCollection
  std::string clusterProducer = conf_.getUntrackedParameter<std::string>("ClusterProducer","siPixelClusters");

  // get input data
  edm::Handle< edm::DetSetVector<SiPixelCluster> >  input;
  iEvent.getByLabel(clusterProducer, input);

  std::map<uint32_t,SiPixelClusterModule*>::iterator struct_iter;
  for (struct_iter = thePixelStructure.begin() ; struct_iter != thePixelStructure.end() ; struct_iter++) {
    
    (*struct_iter).second->fill(*input);
    
  }


  // slow down...
  usleep(100000);
 
  
}

//------------------------------------------------------------------
// Build data structure
//------------------------------------------------------------------
void SiPixelClusterSource::buildStructure(const edm::EventSetup& iSetup){

  //
  // Later on the structure will be built from cabling and
  // *not* from geometry
  //
  std::cout <<" *** SiPixelClusterSource::buildStructure" << std::endl;
  edm::ESHandle<TrackerGeometry> pDD;
  iSetup.get<TrackerDigiGeometryRecord>().get( pDD );
                                                                                                                                                       
  std::cout <<" *** Geometry node for TrackerGeom is  "<<&(*pDD)<<std::endl;
  std::cout <<" *** I have " << pDD->dets().size() <<" detectors"<<std::endl;
  std::cout <<" *** I have " << pDD->detTypes().size() <<" types"<<std::endl;
  
  for(TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++){
    if(dynamic_cast<PixelGeomDetUnit*>((*it))!=0){
      DetId detId = (*it)->geographicalId();
      
      if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {
        //cout << " ---> Adding Barrel Module " <<  detId.rawId() << endl;
	uint32_t id = detId();
	SiPixelClusterModule* pippo = new SiPixelClusterModule(id);
	thePixelStructure.insert(pair<uint32_t,SiPixelClusterModule*> (id,pippo));

      }	else if(detId.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
	//cout << " ---> Adding Endcap Module " <<  detId.rawId() << endl;
	uint32_t id = detId();
	SiPixelClusterModule* pippo = new SiPixelClusterModule(id);
	thePixelStructure.insert(pair<uint32_t,SiPixelClusterModule*> (id,pippo));
      }

    }
  }
  cout << " *** Pixel Structure Size " << thePixelStructure.size() << endl;
}
//------------------------------------------------------------------
// Book MEs
//------------------------------------------------------------------
void SiPixelClusterSource::bookMEs(){

  std::map<uint32_t,SiPixelClusterModule*>::iterator struct_iter;
  string rootDir = "Tracker";
  theDMBE->setVerbose(0);

  for(struct_iter = thePixelStructure.begin(); struct_iter != thePixelStructure.end(); struct_iter++){
    if(DetId::DetId((*struct_iter).first).subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {
      int layer  = PXBDetId::PXBDetId((*struct_iter).first).layer();
      int ladder = PXBDetId::PXBDetId((*struct_iter).first).ladder();
      int module = PXBDetId::PXBDetId((*struct_iter).first).module();
      
      string ssubdet = "PixelBarrel"; 
      char slayer[80];  sprintf(slayer, "Layer_%i",layer);
      char sladder[80]; sprintf(sladder,"Ladder_%02i",ladder);
      char smodule[80]; sprintf(smodule,"Module_%i",module);
      string sfolder = rootDir + "/" + ssubdet + "/" + slayer + "/" + sladder + "/" + smodule;
      theDMBE->setCurrentFolder(sfolder.c_str());
      (*struct_iter).second->book();

    } else 
      if(DetId::DetId((*struct_iter).first).subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
      string ssubdet = "PixelEndcap";
      int side   =  PXFDetId::PXFDetId((*struct_iter).first).side();
      int disk   =  PXFDetId::PXFDetId((*struct_iter).first).disk();
      int blade  =  PXFDetId::PXFDetId((*struct_iter).first).blade();
      int panel  =  PXFDetId::PXFDetId((*struct_iter).first).panel();
      int module =  PXFDetId::PXFDetId((*struct_iter).first).module();
      char sside[80];  sprintf(sside, "Side_%i",side);
      char sdisk[80];  sprintf(sdisk, "Disk_%i",disk);
      char sblade[80]; sprintf(sblade, "Blade_%02i",blade);
      char spanel[80]; sprintf(spanel, "Panel_%i",panel);
      char smodule[80];sprintf(smodule,"Module_%i",module);
      string sfolder = rootDir + "/" + ssubdet + "/" + sside + "/" + sdisk + "/" + sblade + "/" + spanel + "/" + smodule;
      theDMBE->setCurrentFolder(sfolder.c_str());
      (*struct_iter).second->book();
    }
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelClusterSource);
