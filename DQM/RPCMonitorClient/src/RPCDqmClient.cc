// Package:    RPCDqmClient
// Original Author:  Anna Cimmino

#include "DQM/RPCMonitorClient/interface/RPCDqmClient.h"
#include "DQM/RPCMonitorDigi/interface/RPCBookFolderStructure.h"
#include "DQM/RPCMonitorDigi/interface/utils.h"

//include client headers
#include  "DQM/RPCMonitorClient/interface/RPCDeadChannelTest.h"
#include "DQM/RPCMonitorClient/interface/RPCMultiplicityTest.h"
#include "DQM/RPCMonitorClient/interface/RPCClusterSizeTest.h"
#include "DQM/RPCMonitorClient/interface/RPCOccupancyTest.h"
#include "DQM/RPCMonitorClient/interface/RPCNoisyStripTest.h"

//Geometry
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

//Framework
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//DQMServices
#include "DQMServices/Core/interface/MonitorElement.h"

RPCDqmClient::RPCDqmClient(const edm::ParameterSet& iConfig)

{
 edm::LogVerbatim ("rpcdqmclient") << "[RPCDqmClient]: Constructor";

  parameters_ = iConfig;
  
  //check enabling
  enableDQMClients_ =parameters_.getUntrackedParameter<bool> ("EnableRPCDqmClients",true); 
  minimumEvents_= parameters_.getUntrackedParameter<int>("MinimumRPCEvents", 10000);

  //get prescale factor
  prescaleGlobalFactor_ = parameters_.getUntrackedParameter<int>("DiagnosticPrescale", 5);

  prefixDir_ = parameters_.getUntrackedParameter<std::string>("RPCFolder", "RPC");
  recHitType_ =  parameters_.getUntrackedParameter<std::string>("NoiseOrMuons", "Noise");
  std::string gFolder = parameters_.getUntrackedParameter<std::string>("GlobalFolder", "SummaryHistograms");

  globalFolder_ =  prefixDir_+ "/" +  recHitType_ +"/" + gFolder;


  //make default client list  
  clientList_.push_back("RPCMultiplicityTest");
  clientList_.push_back("RPCDeadChannelTest");
  clientList_.push_back("RPCClusterSizeTest");
  clientList_= parameters_.getUntrackedParameter<std::vector<std::string> >("RPCDqmClientList",clientList_);

  //get all the possible RPC DQM clients 
  this->makeClientMap();
}

RPCDqmClient::~RPCDqmClient(){dbe_ = 0;}

void RPCDqmClient::beginJob(){
 edm::LogVerbatim ("rpcdqmclient") << "[RPCDqmClient]: Begin Job";
  if (!enableDQMClients_) return;                 ;

  dbe_ = edm::Service<DQMStore>().operator->();
  dbe_->setVerbose(0);
  
  //Do whatever the begin jobs of all client modules do
  for(std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
   (*it)->beginJob(dbe_);
  
}

void  RPCDqmClient::endRun(const edm::Run& r, const edm::EventSetup& c){
   edm::LogVerbatim ("rpcdqmclient") << "[RPCDqmClient]: End Run";
  if (!enableDQMClients_) return;

  init_ = false;

  std::vector<MonitorElement *>  myMeVect;
  std::vector<RPCDetId>   myDetIds;

  edm::ESHandle<RPCGeometry> rpcGeo;
  c.get<MuonGeometryRecord>().get(rpcGeo);
 
  dbe_->setCurrentFolder(prefixDir_);

 
  //loop on all geometry and get all histos
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
    if( dynamic_cast< RPCChamber* >( *it ) != 0 ){
      
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
      std::vector< const RPCRoll*> roles = (ch->rolls());
      //Loop on rolls in given chamber
      for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
	RPCDetId detId = (*r)->id();
	
	//Get Occupancy ME for roll
	RPCGeomServ RPCname(detId);	   
	RPCBookFolderStructure *  folderStr = new RPCBookFolderStructure();

	//loop on clients
	for( unsigned int cl = 0; cl<clientModules_.size(); cl++ ){

 	  MonitorElement * myMe = dbe_->get(prefixDir_+"/"+ folderStr->folderStructure(detId, recHitType_)+"/"+clientHisto_[cl]+ "_"+RPCname.name()); 

	  if (!myMe || find(myMeVect.begin(), myMeVect.end(), myMe)!=myMeVect.end())continue;

	  dbe_->tag(myMe, clientTag_[cl]);

	  myMeVect.push_back(myMe);
	  myDetIds.push_back(detId);
	}//end loop on clients
      }//end loop on roll in given chamber
    }
  }//end loop on all geometry and get all histos  
  
  for (std::vector<RPCClient*>::iterator  it= clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->endRun(r,c,myMeVect, myDetIds);

  MonitorElement * RPCEvents = dbe_->get(prefixDir_+"/RPCEvents");  

  float   rpcevents = 0;
  if(RPCEvents) rpcevents = RPCEvents -> getIntValue();

  if(rpcevents < minimumEvents_) return;

  for (std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->clientOperation(c);

}

void RPCDqmClient::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& context) {
  if (!enableDQMClients_) return;

  for ( std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->beginLuminosityBlock(lumiSeg,context);
}

void RPCDqmClient::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

 if (!enableDQMClients_) return;

 for ( std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->analyze( iEvent,iSetup);
}


void RPCDqmClient::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& c){

  edm::LogVerbatim ("rpcdqmclient") <<"[RPCDqmClient]: End of LS ";
 
  if (!enableDQMClients_ ) return;
    
  for (std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->endLuminosityBlock( lumiSeg, c);
}


void  RPCDqmClient::beginRun(const edm::Run& r, const edm::EventSetup& c){

 if (!enableDQMClients_) return;
   for ( std::vector<RPCClient*>::iterator it = clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->beginRun(r,c);
}


void RPCDqmClient::endJob() {
 if (!enableDQMClients_) return;
 
 for ( std::vector<RPCClient*>::iterator it= clientModules_.begin(); it!=clientModules_.end(); it++ )
    (*it)->endJob();
}


void RPCDqmClient::makeClientMap() {
  
  std::vector< std::string>  clientList,clientNames,clientHisto; 
  std::vector<RPCClient*> clientModules;
  std::vector<int> clientTag;
  
  //clear global vectors;
  clientNames_.clear();
  clientHisto_.clear();
  clientTag_.clear();
  clientModules_.clear();

  if (clientList_.size()==0) return; //if no client is selected by user, return
 
  //Fill vectors with all possible RPC DQM clients , source histos names, and tag values
  //RPCMultiplicityTest
  clientNames.push_back("RPCMultiplicityTest");
  clientHisto.push_back("Multiplicity");
  clientTag.push_back(rpcdqm::MULTIPLICITY);
  clientModules.push_back( new RPCMultiplicityTest(parameters_));
  //RPCDeadChannelTest
  clientNames.push_back("RPCDeadChannelTest");
  clientHisto.push_back("Occupancy");
  clientModules.push_back( new RPCDeadChannelTest(parameters_));
  clientTag.push_back(rpcdqm::OCCUPANCY);
  //RPCClusterSizeTest
  clientNames.push_back("RPCClusterSizeTest");
  clientHisto.push_back("ClusterSize");
  clientTag.push_back(rpcdqm::CLUSTERSIZE);
  clientModules.push_back( new RPCClusterSizeTest(parameters_));
  //RPCOccupancyTest
  clientNames.push_back("RPCOccupancyTest");
  clientHisto.push_back("Occupancy");
  clientTag.push_back(rpcdqm::OCCUPANCY);
  clientModules.push_back( new RPCOccupancyTest(parameters_));
 //RPCNoisyStripTest
  clientNames.push_back("RPCNoisyStripTest");
  clientHisto.push_back("Occupancy");
  clientTag.push_back(rpcdqm::OCCUPANCY);
  clientModules.push_back( new RPCNoisyStripTest(parameters_));


  //take only user specified clients and associate its source histograms to it
  for(unsigned int i = 0; i<clientNames.size(); i++){

    if(find(clientList_.begin(),clientList_.end(),clientNames[i])!=clientList_.end()) {

      clientHisto_.push_back(clientHisto[i]);
      clientTag_.push_back(clientTag[i]);
      clientNames_.push_back(clientNames[i]);
      clientModules_.push_back(clientModules[i]);
    }
  }
  
  return;
}
