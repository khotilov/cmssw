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

  globalFolder_ = ps.getUntrackedParameter<string>("RPCGlobalFolder", "RPC/RecHits/SummaryHistograms");
  prescaleFactor_ = ps.getUntrackedParameter<int>("DiagnosticPrescale", 1);
  numberOfDisks_ = ps.getUntrackedParameter<int>("NumberOfEndcapDisks", 3);
}

RPCDeadChannelTest::~RPCDeadChannelTest(){dbe_ = 0;}

void RPCDeadChannelTest::beginJob(DQMStore *  dbe ){
  dbe_=dbe;
}

void RPCDeadChannelTest::beginRun(const Run& r, const EventSetup& iSetup,vector<MonitorElement *> meVector, vector<RPCDetId> detIdVector){

 edm::LogVerbatim ("deadChannel") << "[RPCDeadChannelTest]: Begin run";

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
     if ( me = dbe_->get(globalFolder_ +"/"+ histoName.str()) ) {
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
   histoName<<"DeadChannelFraction_Roll_vs_Sector_Disk"<<i;
   if ( me = dbe_->get(globalFolder_ +"/"+ histoName.str()) ) {
     dbe_->removeElement(me->getName());
   }
  
   DEADDisk[i+offset] = dbe_->book2D(histoName.str().c_str(), histoName.str().c_str(), 6, 0.5, 6.5, 54, 0.5, 54.5);
   
   for (int x = 1; x<=6; x++)
     for(int y=1; y<=54; y++)
       DEADDisk[i+offset]->setBinContent(x,y,-1);
   rpcUtils.labelXAxisSector( DEADDisk[i+offset]);
   rpcUtils.labelYAxisRoll( DEADDisk[i+offset], 1, i);
  
 }//end loop on wheels and disks

 //Get Occuoancy ME for each roll
  
 for (unsigned int i = 0 ; i<meVector.size(); i++){

   bool flag= false;
   
   DQMNet::TagList tagList;
   tagList = meVector[i]->getTags();
   DQMNet::TagList::iterator tagItr = tagList.begin();

   while (tagItr != tagList.end() && !flag ) {
     if((*tagItr) ==  rpcdqm::OCCUPANCY)
       flag= true;
   
     tagItr++;
   }
   
   if(flag){
      myOccupancyMe_.push_back(meVector[i]);
      myDetIds_.push_back(detIdVector[i]);
   }
 }  
}

void RPCDeadChannelTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {}

void RPCDeadChannelTest::analyze(const edm::Event& iEvent, const edm::EventSetup& c){}

void RPCDeadChannelTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& iSetup) {
 
  edm::LogVerbatim ("deadChannel") <<"[RPCDeadChannelTest]: End of LS transition, performing the DQM client operation";
  
  // counts number of lumiSegs 
  int nLumiSegs = lumiSeg.id().luminosityBlock();

  //check some statements and prescale Factor
  if(nLumiSegs%prescaleFactor_ != 0) return;

    //Loop on chambers
    for (unsigned int  i = 0 ; i<myOccupancyMe_.size();i++){
      this->CalculateDeadChannelPercentage(myDetIds_[i],myOccupancyMe_[i],iSetup);
    }//End loop on rolls in given chambers

}
 
void RPCDeadChannelTest::endRun(const Run& r, const EventSetup& c){}

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
       if(((detId.station() * detId.region() ) + numberOfDisks_) >= 0 ){
	 
	 if(detId.region()<0){
	   DEAD  = DEADDisk[(detId.station() * detId.region() ) + numberOfDisks_];
	 }else{
	   DEAD = DEADDisk[(detId.station() * detId.region() ) + numberOfDisks_-1];
	 }
       }
     }
     
     if (DEAD){
       rpcdqm::utils rollNumber;
       int nr = rollNumber.detId2RollNr(detId);
       DEAD->setBinContent(detId.sector(),nr, (float)badChannels.size()/nstrips );
     }
   }
}


