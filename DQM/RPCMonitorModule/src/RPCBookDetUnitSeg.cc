#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <map>


#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <DQM/RPCMonitorModule/interface/MuonSegmentEff.h>
#include <DQMOffline/Muon/interface/RPCBookFolderStructure.h>
#include "DQMServices/Core/interface/MonitorElement.h"

#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
std::map<std::string, MonitorElement*> MuonSegmentEff::bookDetUnitSeg(RPCDetId & detId,int nstrips,float stripw,float stripl) {
  
  std::map<std::string, MonitorElement*> meMap;
   
  RPCBookFolderStructure *  folderStr = new RPCBookFolderStructure(); //Anna
  std::string folder = "Muons/MuonSegEff/" +  folderStr->folderStructure(detId);

  dbe->setCurrentFolder(folder);

  RPCGeomServ RPCname(detId);
  std::string nameRoll = RPCname.name();

  char detUnitLabel[128];
  char layerLabel[128];

  sprintf(detUnitLabel ,"%s",nameRoll.c_str());
  sprintf(layerLabel ,"%s",nameRoll.c_str());

  char meId [128];
  char meTitle [128];

  float scale2D = 0.6;

  //Begin booking DT
  if(detId.region()==0) {
    
    sprintf(meId,"ExpectedOccupancyFromDT_%s",detUnitLabel);
    sprintf(meTitle,"ExpectedOccupancyFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"RealDetectedOccupancyFromDT_%s",detUnitLabel);
    sprintf(meTitle,"RealDetectedOccupancyFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"RPCDataOccupancyFromDT_%s",detUnitLabel);
    sprintf(meTitle,"RPCDataOccupancyFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"BXDistribution_%s",detUnitLabel);
    sprintf(meTitle,"BXDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, 11,-5.5, 5.5);

    sprintf(meId,"CLSDistribution_%s",detUnitLabel);
    sprintf(meTitle,"CLSDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, 11,0.5,10.5);

    sprintf(meId,"BXYDistribution_%s",detUnitLabel);
    sprintf(meTitle,"BXYDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle,7,-3.5,3.5,20,0,stripl);   

    //New 2D and more

    sprintf(meId,"ExpectedOccupancy2DFromDT_%s",detUnitLabel);
    sprintf(meTitle,"ExpectedOccupancy2DFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"RPCDataOccupancy2DFromDT_%s",detUnitLabel);
    sprintf(meTitle,"RPCDataOccupancy2DFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"Inefficiency2DFromDT_%s",detUnitLabel);
    sprintf(meTitle,"Inefficiency2DFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"RPCResidualsFromDT_%s",detUnitLabel);
    sprintf(meTitle,"RPCResidualsFromDT_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle,101,-7.*stripw,7*stripw);
  
  }else{
    //std::cout<<"Booking for the EndCap"<<detUnitLabel<<std::endl;

    sprintf(meId,"ExpectedOccupancyFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"ExpectedOccupancyFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"RealDetectedOccupancyFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"RealDetectedOccupancyFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"RPCDataOccupancyFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"RPCDataOccupancyFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, nstrips, 0.5, nstrips+0.5);
    
    sprintf(meId,"BXDistribution_%s",detUnitLabel);
    sprintf(meTitle,"BXDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, 11,-5.5, 5.5);

    sprintf(meId,"CLSDistribution_%s",detUnitLabel);
    sprintf(meTitle,"CLSDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle, 11,0.5,10.5);

    sprintf(meId,"BXYDistribution_%s",detUnitLabel);
    sprintf(meTitle,"BXYDistribution_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle,7,-3.5,3.5,20,0,stripl);   
    
    //New 2D and more

    sprintf(meId,"ExpectedOccupancy2DFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"ExpectedOccupancy2DFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"RPCDataOccupancy2DFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"RPCDataOccupancy2DFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"Inefficiency2DFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"Inefficiency2DFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book2D(meId, meTitle, 2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

    sprintf(meId,"RPCResidualsFromCSC_%s",detUnitLabel);
    sprintf(meTitle,"RPCResidualsFromCSC_for_%s",layerLabel);
    meMap[meId] = dbe->book1D(meId, meTitle,101,-7.*stripw,7*stripw);
   
  }
  return meMap;
}



