/*
 *  \author Anna Cimmino
 */
#include "DQM/RPCMonitorDigi/interface/utils.h"
#include <DQM/RPCMonitorClient/interface/RPCDeadChannelTest.h>

// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <FWCore/Framework/interface/ESHandle.h>

// Geometry
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include <sstream>

using namespace edm;
using namespace std;

RPCDeadChannelTest::RPCDeadChannelTest(const ParameterSet& ps ){
 
  LogVerbatim ("deadChannel") << "[RPCDeadChannelTest]: Constructor";

 std::string prefixDir = ps.getUntrackedParameter<std::string>("RPCFolder", "RPC");
  std::string recHitType =  ps.getUntrackedParameter<std::string>("NoiseOrMuons", "Noise");
  std::string gFolder = ps.getUntrackedParameter<std::string>("GlobalFolder", "SummaryHistograms");

  globalFolder_ =  prefixDir + "/" +  recHitType +"/" + gFolder;

  prescaleFactor_ = ps.getUntrackedParameter<int>("DiagnosticPrescale", 1);
  numberOfDisks_ = ps.getUntrackedParameter<int>("NumberOfEndcapDisks", 3);
  numberOfRings_ = ps.getUntrackedParameter<int>("NumberOfEndcapRings", 2);
}

RPCDeadChannelTest::~RPCDeadChannelTest(){dbe_ = 0;}

void RPCDeadChannelTest::beginJob(DQMStore *  dbe ){
  dbe_=dbe;
}

void RPCDeadChannelTest::endRun(const Run& r, const EventSetup& iSetup,vector<MonitorElement *> meVector, vector<RPCDetId> detIdVector){

 edm::LogVerbatim ("deadChannel") << "[RPCDeadChannelTest]: End run";

 MonitorElement* me;
 dbe_->setCurrentFolder( globalFolder_);

 stringstream histoName;

 rpcdqm::utils rpcUtils;
  
 int limit = numberOfDisks_;
 if(numberOfDisks_ < 2) limit = 2;
  
 for (int i = -1 * limit; i<= limit;i++ ){//loop on wheels and disks
   if (i>-3 && i<3){//wheels
     histoName.str("");
     histoName<<"DeadChannelFraction_Roll_vs_Sector_Wheel"<<i;
     me = 0;
     me = dbe_->get(globalFolder_ +"/"+ histoName.str());
     if (0!=me ) {
       dbe_->removeElement(me->getName());
     }
     DEADWheel[i+2] = dbe_->book2D(histoName.str().c_str(), histoName.str().c_str(), 12, 0.5, 12.5, 21, 0.5, 21.5);

     for (int x = 1; x<=12; x++)
       for(int y=1; y<=21; y++)
	 DEADWheel[i+2]->setBinContent(x,y,-1);

     rpcUtils.labelXAxisSector( DEADWheel[i+2]);
     rpcUtils.labelYAxisRoll( DEADWheel[i+2], 0, i);
   }//end wheels
     
   if (i == 0  || i > numberOfDisks_ || i< (-1 * numberOfDisks_))continue;
  
   int offset = numberOfDisks_;
   if (i>0) offset --; //used to skip case equale to zero
  
   histoName.str("");
   histoName<<"DeadChannelFraction_Ring_vs_Segment_Disk"<<i;
   me = 0;
   me = dbe_->get(globalFolder_ +"/"+ histoName.str());
   if (0!=me ) {
     dbe_->removeElement(me->getName());
   }
  
   DEADDisk[i+offset] = dbe_->book2D(histoName.str().c_str(), histoName.str().c_str(),36, 0.5, 36.5, 3*numberOfRings_, 0.5,3*numberOfRings_+ 0.5);
   
   rpcUtils.labelXAxisSegment(DEADDisk[i+offset]);
   rpcUtils.labelYAxisRing(DEADDisk[i+offset], numberOfRings_);

  
 }//end loop on wheels and disks

 //Get Occuoancy ME for each roll
  
 for (unsigned int i = 0 ; i<meVector.size(); i++){

   DQMNet::TagList tagList;
   tagList = meVector[i]->getTags();
   DQMNet::TagList::iterator tagItr = tagList.begin();

   while (tagItr != tagList.end() ) {
     if((*tagItr) ==  rpcdqm::OCCUPANCY){
      myOccupancyMe_.push_back(meVector[i]);
      myDetIds_.push_back(detIdVector[i]);
      break;
     }
       
     tagItr++;
   }
   
 }  
}

void RPCDeadChannelTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {}

void RPCDeadChannelTest::analyze(const edm::Event& iEvent, const edm::EventSetup& c){}

void RPCDeadChannelTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& iSetup){}

void RPCDeadChannelTest::clientOperation( EventSetup const& iSetup){
 
  edm::LogVerbatim ("deadChannel") <<"[RPCDeadChannelTest]:Client Operation";
  
  //Loop on chambers
    for (unsigned int  i = 0 ; i<myOccupancyMe_.size();i++){
      this->CalculateDeadChannelPercentage(myDetIds_[i],myOccupancyMe_[i],iSetup);
    }//End loop on rolls in given chambers

}
 
void RPCDeadChannelTest::beginRun(const Run& r, const EventSetup& c){}

void RPCDeadChannelTest::endJob(){}

//
//User Defined methods
//
void  RPCDeadChannelTest::CalculateDeadChannelPercentage(RPCDetId & detId, MonitorElement * myMe, EventSetup const& iSetup){

  ESHandle<RPCGeometry> rpcgeo;
  iSetup.get<MuonGeometryRecord>().get(rpcgeo); 

  const RPCRoll * rpcRoll = rpcgeo->roll(detId);      

  unsigned int nstrips =rpcRoll->nstrips();

  MonitorElement * DEAD = NULL;

   const QReport * theOccupancyQReport = myMe->getQReport("DeadChannel_0");  
   if(theOccupancyQReport) {
     
     vector<dqm::me_util::Channel> badChannels = theOccupancyQReport->getBadChannels();
     
     if (detId.region()==0)   DEAD = DEADWheel[detId.ring() + 2] ;
     else{
       if(-detId.station()+ numberOfDisks_ >= 0 ){
	 
	 if(detId.region()<0){
	   DEAD  = DEADDisk[-detId.station() + numberOfDisks_];
	 }else{
	   DEAD = DEADDisk[detId.station() + numberOfDisks_-1];
	 }
       }
     }

     if (DEAD){
       int xBin,yBin;
       if(detId.region()==0){//Barrel
	 xBin= detId.sector();
	 rpcdqm::utils rollNumber;
	 yBin = rollNumber.detId2RollNr(detId);
       }else{//Endcap
	 //get segment number
	 RPCGeomServ RPCServ(detId);
	 xBin = RPCServ.segment();
	 (numberOfRings_ == 3 ? yBin= detId.ring()*3-detId.roll()+1 : yBin= (detId.ring()-1)*3-detId.roll()+1);
     }
       DEAD->setBinContent(xBin,yBin, (float)badChannels.size()/nstrips );

     }


     
 
   }
}


