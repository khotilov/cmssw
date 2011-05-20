#include <DQM/RPCMonitorDigi/interface/RPCMonitorDigi.h>
#include <DQM/RPCMonitorDigi/interface/RPCBookFolderStructure.h>
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <DQM/RPCMonitorDigi/interface/utils.h>
#include <iomanip>

void RPCMonitorDigi::bookRollME(RPCDetId & detId, const edm::EventSetup & iSetup, const std::string & recHitType, std::map<std::string, MonitorElement*>  & meMap) {
  //std::map<std::string, MonitorElement*> RPCMonitorDigi::bookRollME(RPCDetId & detId, const edm::EventSetup & iSetup, std::string recHitType) {
  //std::map<std::string, MonitorElement*> meMap;  

  RPCBookFolderStructure *  folderStr = new RPCBookFolderStructure();
  std::string folder = subsystemFolder_+ "/"+ recHitType +"/"+folderStr->folderStructure(detId);

  dbe->setCurrentFolder(folder);
  
  //get number of strips in current roll
  int nstrips = this->stripsInRoll(detId, iSetup);
  if (nstrips == 0 ) nstrips = 1;

  /// Name components common to current RPCDetId  
  RPCGeomServ RPCname(detId);
  std::string nameRoll = RPCname.name();
  
  std::stringstream os;
  os.str("");
  os<<"Occupancy_"<<nameRoll;
  meMap[os.str()] = dbe->book1D(os.str(), os.str(), nstrips, 0.5, nstrips+0.5);
  dbe->tag( meMap[os.str()],  rpcdqm::OCCUPANCY);

  os.str("");
  os<<"BXDistribution_"<<nameRoll;
  meMap[os.str()] = dbe->book1D(os.str(), os.str(), 9, -4.5, 4.5);
  
  os.str("");
  os<<"ClusterSize_"<<nameRoll;
  meMap[os.str()] = dbe->book1D(os.str(), os.str(), 20, 0.5, 20.5);
  dbe->tag( meMap[os.str()],  rpcdqm::CLUSTERSIZE);
  
  os.str("");
  os<<"Multiplicity_"<<nameRoll;
  meMap[os.str()] = dbe->book1D(os.str(), os.str(), 50, 0.5, 50.5);
  dbe->tag( meMap[os.str()],  rpcdqm::MULTIPLICITY);
  
//   os.str("");
//   os<<"BXWithData_"<<nameRoll;
//   meMap[os.str()] = dbe->book1D(os.str(), os.str(), 10, 0.5, 10.5);
  
  os.str("");
  os<<"NumberOfClusters_"<<nameRoll;
  meMap[os.str()] = dbe->book1D(os.str(), os.str(),10,0.5,10.5);
  

  delete folderStr;
  //  return meMap;
}


void RPCMonitorDigi::bookSectorRingME(const std::string &recHitType, std::map<std::string, MonitorElement*> & meMap) {  
  //std::map<std::string, MonitorElement*> RPCMonitorDigi::bookSectorRingME(std::string recHitType) {  

  //  std::map<std::string, MonitorElement*> meMap;  
  std::stringstream os;
 
  for(int wheel = -2 ; wheel <= 2; wheel++){
      os.str("");     
      os<< subsystemFolder_<< "/"<<recHitType<<"/Barrel/Wheel_"<<wheel<<"/SummaryBySectors";
      dbe->setCurrentFolder(os.str());
      
      for (int sector = 1 ; sector <= 12 ; sector++){
	
	os.str("");
	os<<"Occupancy_Wheel_"<<wheel<<"_Sector_"<<sector;
    
	if (sector==9 || sector==11)
	  meMap[os.str()] = dbe->book2D(os.str(), os.str(), 96, 0.5,96.5, 15, 0.5, 15.5);
	else  if (sector==4) 
	  meMap[os.str()] = dbe->book2D(os.str(), os.str(),  96, 0.5, 96.5, 21, 0.5, 21.5);
	else
	  meMap[os.str()] = dbe->book2D(os.str(), os.str(), 96, 0.5,  96.5, 17, 0.5, 17.5);
	
	meMap[os.str()]->setAxisTitle("strip", 1);
	rpcdqm::utils rpcUtils;
	rpcUtils.labelYAxisRoll( meMap[os.str()], 0, wheel);
	
// 	os.str("");
// 	os<<"BxDistribution_Wheel_"<<wheel<<"_Sector_"<<sector;
// 	meMap[os.str()] = dbe->book1D(os.str(), os.str(), 11, -5.5, 5.5);

      }
  }


  for (int region = -1 ; region <= 1; region++){
    if( region == 0 ) continue;

    std::string regionName = "Endcap-";
    if(region == 1) regionName = "Endcap+";

    for (int disk = 1; disk <=  RPCMonitorDigi::numberOfDisks_; disk++) {
      os.str("");
      os<< subsystemFolder_<< "/"<<recHitType<<"/"<<regionName<<"/Disk_"<<(region * disk)<<"/SummaryByRings/";
     
      dbe->setCurrentFolder(os.str());

      for (int ring = RPCMonitorDigi::numberOfInnerRings_  ; ring <= 3; ring ++) {

	os.str("");
	os<<"Occupancy_Disk_"<<(region * disk)<<"_Ring_"<<ring<<"_CH01-CH12";

	meMap[os.str()] = dbe->book2D(os.str(), os.str(), 96, 0.5, 96.5, 12 , 0.5,  12.5);
	meMap[os.str()]->setAxisTitle("strip", 1);

	std::stringstream yLabel;
	for (int i = 1 ; i<=12; i++) {
	  yLabel.str("");
	  yLabel<<"R"<<ring<<"_CH"<<std::setw(2)<<std::setfill('0')<<i;
	  meMap[os.str()]->setBinLabel(i, yLabel.str(), 2);
	}
      
      
	for(int i = 1; i <= 96 ; i++) {
	  if (i ==1) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==16) meMap[os.str()]->setBinLabel(i, "RollA", 1);
	  else if (i==32) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else if (i==33) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==48) meMap[os.str()]->setBinLabel(i, "RollB", 1);
	  else if (i==64) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else if (i==65) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==80) meMap[os.str()]->setBinLabel(i, "RollC", 1);
	  else if (i==96) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else  meMap[os.str()]->setBinLabel(i, "", 1);
	}
  

	os.str("");
	os<<"Occupancy_Disk_"<<(region * disk)<<"_Ring_"<<ring<<"_CH13-CH36";

	meMap[os.str()] = dbe->book2D(os.str(), os.str(), 96, 0.5, 96.5, 12 , 12.5,  36.5);
	meMap[os.str()]->setAxisTitle("strip", 1);
	
	for (int i = 1 ; i<= 12; i++) {
	  yLabel.str("");
	  yLabel<<"R"<<ring<<"_CH"<<i+12;
	  meMap[os.str()]->setBinLabel(i, yLabel.str(), 2);
	}
	
	
	for(int i = 1; i <= 96 ; i++) {
	  if (i ==1) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==16) meMap[os.str()]->setBinLabel(i, "RollA", 1);
	  else if (i==32) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else if (i==33) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==48) meMap[os.str()]->setBinLabel(i, "RollB", 1);
	  else if (i==64) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else if (i==65) meMap[os.str()]->setBinLabel(i, "1", 1);
	  else if (i==80) meMap[os.str()]->setBinLabel(i, "RollC", 1);
	  else if (i==96) meMap[os.str()]->setBinLabel(i, "32", 1);
	  else  meMap[os.str()]->setBinLabel(i, "", 1);
	}
   
        
// 	os.str("");
// 	os<<"BxDistribution_Disk_"<<(region * disk)<<"_Ring_"<<ring;
// 	meMap[os.str()] = dbe->book1D(os.str(), os.str(), 11, -5.5, 5.5);
	
      }  //loop ring
    } //loop disk
  } //loop region

  // return meMap;
} 


void RPCMonitorDigi::bookWheelDiskME(const std::string &recHitType, std::map<std::string, MonitorElement*> &meMap) {  
  //std::map<std::string, MonitorElement*> RPCMonitorDigi::bookWheelDiskME(std::string recHitType) {  

  //  std::map<std::string, MonitorElement*> meMap;  
  dbe->setCurrentFolder(subsystemFolder_ +"/"+recHitType+"/"+ globalFolder_);

  std::stringstream os, label, name, title ;
  rpcdqm::utils rpcUtils;

  for (int wheel = -2 ; wheel<= 2; wheel++ ) {//Loop on wheel

    //    os<<"OccupancyXY_"<<ringType<<"_"<<ring;
    //    meMap[os.str()] = dbe->book2D(os.str(), os.str(),63, -800, 800, 63, -800, 800);
    //    meMap[os.str()] = dbe->book2D(os.str(), os.str(),1000, -800, 800, 1000, -800, 800);
    
    
    os.str("");
    os<<"1DOccupancy_Wheel_"<<wheel;
    meMap[os.str()] = dbe->book1D(os.str(), os.str(), 12, 0.5, 12.5);
    for(int i=1; i<12; i++) {
      label.str("");
      label<<"Sec"<<i;
      meMap[os.str()] ->setBinLabel(i, label.str(), 1); 
    }
    
    os.str("");
    os<<"Occupancy_Roll_vs_Sector_Wheel_"<<wheel;                                   
    meMap[os.str()] = dbe->book2D(os.str(), os.str(), 12, 0.5,12.5, 21, 0.5, 21.5);
    rpcUtils.labelXAxisSector(meMap[os.str()]);
    rpcUtils.labelYAxisRoll( meMap[os.str()], 0, wheel);

    os.str("");
    os<<"BxDistribution_Wheel_"<<wheel;
    meMap[os.str()] = dbe->book1D(os.str(), os.str(), 9, -4.5, 4.5);
    

    for(int layer = 1 ; layer <= 6 ; layer ++){
      name.str("");
      title.str("");
      name<<"ClusterSize_Wheel_"<<wheel<<"_Layer"<< layer;
      title<< "ClusterSize - Wheel "<<wheel<<" Layer"<<layer;
      meMap[name.str()] = dbe->book1D(name.str(), title.str(),  20, 0.5, 20.5);
    }


    
  }//end loop on wheel 


  for (int disk = - RPCMonitorDigi::numberOfDisks_; disk <=  RPCMonitorDigi::numberOfDisks_; disk++){
    
    if(disk == 0) continue;
  

    os.str("");
    os<<"Occupancy_Ring_vs_Segment_Disk_"<<disk;                                  
    meMap[os.str()] = dbe->book2D(os.str(), os.str(), 36, 0.5,36.5, 6, 0.5, 6.5);
    
    rpcUtils.labelXAxisSegment(meMap[os.str()]);
    rpcUtils.labelYAxisRing(meMap[os.str()], 2);

    os.str("");
    os<<"BxDistribution_Disk_"<<disk;
    meMap[os.str()] = dbe->book1D(os.str(), os.str(), 9, -4.5, 4.5);


    for(int ring = RPCMonitorDigi::numberOfInnerRings_  ; ring <= 3 ; ring ++){
    
      name.str("");
      title.str("");
      name<<"ClusterSize_Disk_"<<disk<<"_Ring"<< ring;
      title<< "ClusterSize - Disk"<<disk<<" Ring"<<ring;
      meMap[name.str()] = dbe->book1D(name.str(), title.str(),  20, 0.5, 20.5);
      
    }
    
  }

   for(int ring = RPCMonitorDigi::numberOfInnerRings_  ; ring <= 3 ; ring ++){
     os.str("");
     os<<"1DOccupancy_Ring_"<<ring;
     meMap[os.str()] = dbe->book1D(os.str(), os.str(), RPCMonitorDigi::numberOfDisks_ + 1, -(RPCMonitorDigi::numberOfDisks_ + 0.5), (RPCMonitorDigi::numberOfDisks_ + 0.5));
     for(int xbin= 1 ; xbin<=(RPCMonitorDigi::numberOfDisks_ *2) ; xbin++) {
       label.str("");
       if ((xbin - 4)!=0) label<<"Disk "<< (xbin - 4);
       else label<<"-";
       meMap[os.str()] ->setBinLabel(xbin, label.str(), 1); 
     }
   }



      
  //  return meMap; 
}



//returns the number of strips in each roll
int  RPCMonitorDigi::stripsInRoll(RPCDetId & id, const edm::EventSetup & iSetup) {
  edm::ESHandle<RPCGeometry> rpcgeo;
  iSetup.get<MuonGeometryRecord>().get(rpcgeo);

  const RPCRoll * rpcRoll = rpcgeo->roll(id);

  if (rpcRoll)
    return  rpcRoll->nstrips();
  else 
    return 1;
}


void RPCMonitorDigi::bookRegionME(const std::string & recHitType, std::map<std::string, MonitorElement*>  & meMap) {
  //std::map<std::string, MonitorElement*>   RPCMonitorDigi::bookRegionME(std::string recHitType) {

  //  std::map<std::string, MonitorElement*> meMap;  

  std::string currentFolder = subsystemFolder_ +"/"+recHitType+"/"+ globalFolder_;
  dbe->setCurrentFolder(currentFolder);  
  
  MonitorElement * me = NULL;
  std::stringstream name;
  std::stringstream title;
  for(int r = 0; r < 3; r++){ //RPC regions are E-, B, and E+
    
    std::string regionName = RPCMonitorDigi::regionNames_[r];
    //Cluster size
    name.str("");
    title.str("");
    name<<"ClusterSize_"<< regionName;
    title<< "ClusterSize - "<<regionName;
    me = dbe->get(currentFolder+ "/" + name.str());
    if (me) dbe->removeElement(me->getName());
     meMap[name.str()] = dbe->book1D(name.str(), title.str(),  20, 0.5, 20.5);
    
    //Number of Cluster
    name.str("");
    title.str("");
    name<<"NumberOfClusters_"<< regionName;
    title<< "Number of Clusters per Event - "<< regionName;
    me = dbe->get(currentFolder+ "/" + name.str());
    if (me) dbe->removeElement(me->getName());
    meMap[name.str()]  = dbe->book1D(name.str(), title.str(),  30, 0.5, 30.5);
    
    //Number of Digis
    name.str("");
    title.str("");
    name<<"Multiplicity_"<< regionName;
    title<< "Multiplicity per Event per Roll - "<< regionName;
    me = dbe->get(currentFolder+ "/" + name.str());
    if (me) dbe->removeElement(me->getName());
    meMap[name.str()] = dbe->book1D(name.str(), title.str(), 50, 0.5, 50.5);   
           
  }//end loop on regions



  for(int layer = 1 ; layer <= 6 ; layer ++){
    
    name.str("");
    title.str("");
    name<<"ClusterSize_Layer"<< layer;
    title<< "ClusterSize - Layer"<<layer;
    me = dbe->get(currentFolder+ "/" + name.str());
    if (me) dbe->removeElement(me->getName());
    meMap[name.str()] = dbe->book1D(name.str(), title.str(),  20, 0.5, 20.5);
  }

  for(int ring = RPCMonitorDigi::numberOfInnerRings_  ; ring <= 3 ; ring ++){
    
    name.str("");
    title.str("");
    name<<"ClusterSize_Ring"<< ring;
    title<< "ClusterSize - Ring"<<ring;
    me = dbe->get(currentFolder+ "/" + name.str());
    if (me) dbe->removeElement(me->getName());
    meMap[name.str()] = dbe->book1D(name.str(), title.str(),  20, 0.5, 20.5);
  
  }


  me = dbe->get(currentFolder+ "/Occupancy_for_Endcap");
  if (me) dbe->removeElement(me->getName());
  meMap["Occupancy_for_Endcap"] = dbe -> book2D("Occupancy_for_Endcap", "Occupancy Endcap", 6, 0.5, 6.5, 2, 1.5, 3.5);
  meMap["Occupancy_for_Endcap"] ->setAxisTitle("Disk", 1);
  meMap["Occupancy_for_Endcap"] ->setAxisTitle("Ring", 2);



  me = dbe->get(currentFolder+ "/Occupancy_for_Barrel");
  if (me) dbe->removeElement(me->getName());
  meMap["Occupancy_for_Barrel"]  = dbe -> book2D("Occupancy_for_Barrel", "Occupancy Barrel", 12, 0.5 , 12.5, 5, -2.5, 2.5 );
  meMap["Occupancy_for_Barrel"] ->setAxisTitle("Sec", 1);
  meMap["Occupancy_for_Barrel"] ->setAxisTitle("Wheel", 2);


  //  return meMap; 

}
