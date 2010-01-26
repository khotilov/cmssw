#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DQM/SiStripMonitorClient/interface/SiStripDaqInfo.h"
#include "DQM/SiStripMonitorClient/interface/SiStripUtility.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"

//Run Info
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
#include "CondFormats/RunInfo/interface/RunSummary.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"
// FED cabling and numbering
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/SiStripObjects/interface/SiStripFedCabling.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;
//
// -- Contructor
//
SiStripDaqInfo::SiStripDaqInfo(edm::ParameterSet const& pSet) {
  // Create MessageSender
  edm::LogInfo( "SiStripDaqInfo") << "SiStripDaqInfo::Deleting SiStripDaqInfo ";

  // get back-end interface
  dqmStore_ = edm::Service<DQMStore>().operator->();
  nFedTotal = 0;
  bookedStatus_ = false;
}
SiStripDaqInfo::~SiStripDaqInfo() {
  edm::LogInfo("SiStripDaqInfo") << "SiStripDaqInfo::Deleting SiStripDaqInfo ";

}
//
// -- Begin Job
//
void SiStripDaqInfo::beginJob() {
 
}
//
// -- Book MEs for SiStrip Daq Fraction
//
void SiStripDaqInfo::bookStatus() {
   edm::LogInfo( "SiStripDcsInfo") << " SiStripDaqInfo::bookStatus " << bookedStatus_;
  if (!bookedStatus_) {
    string strip_dir = "";
    SiStripUtility::getTopFolderPath(dqmStore_, "SiStrip", strip_dir);
    if (strip_dir.size() > 0) {

      dqmStore_->setCurrentFolder(strip_dir+"/EventInfo");
    
      DaqFraction_= dqmStore_->bookFloat("DAQSummary");  
      
      dqmStore_->setCurrentFolder(strip_dir+"/EventInfo/DAQContents");
      
      vector<string> det_type;
      det_type.push_back("TIB");
      det_type.push_back("TOB");
      det_type.push_back("TIDF");
      det_type.push_back("TIDB");
      det_type.push_back("TECF");
      det_type.push_back("TECB");
      
      for ( vector<string>::iterator it = det_type.begin(); it != det_type.end(); it++) {
	string det = (*it);
	
	SubDetMEs local_mes;	
	string me_name;
	me_name = "SiStrip_" + det;    
	local_mes.DaqFractionME = dqmStore_->bookFloat(me_name);  	
	local_mes.ConnectedFeds = 0;
	SubDetMEsMap.insert(make_pair(det, local_mes));
      } 
      bookedStatus_ = true;
      dqmStore_->cd();
    }
  }
}
//
// -- Fill with Dummy values
//
void SiStripDaqInfo::fillDummyStatus() {
  if (!bookedStatus_) bookStatus();
  if (bookedStatus_) {
    for (map<string, SubDetMEs>::iterator it = SubDetMEsMap.begin(); it != SubDetMEsMap.end(); it++) {
      it->second.DaqFractionME->Reset();
      it->second.DaqFractionME->Fill(-1.0);
    }
    DaqFraction_->Reset();
    DaqFraction_->Fill(-1.0);
  }
}
//
// -- Begin Run
//
void SiStripDaqInfo::beginRun(edm::Run const& run, edm::EventSetup const& eSetup) {
  edm::LogInfo ("SiStripDaqInfo") <<"SiStripDaqInfo:: Begining of Run";

  // Check latest Fed cabling and create TrackerMapCreator
  unsigned long long cacheID = eSetup.get<SiStripFedCablingRcd>().cacheIdentifier();
  if (m_cacheID_ != cacheID) {
    m_cacheID_ = cacheID;       

    eSetup.get<SiStripFedCablingRcd>().get(fedCabling_); 

    readFedIds(fedCabling_);
  }
  if (!bookedStatus_) bookStatus();  
}
//
// -- Analyze
//
void SiStripDaqInfo::analyze(edm::Event const& event, edm::EventSetup const& eSetup) {
}

//
// -- End Luminosity Block
//
void SiStripDaqInfo::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& iSetup) {
  edm::LogInfo( "SiStripDaqInfo") << "SiStripDaqInfo::endLuminosityBlock";
}
//
// -- End Run
//
void SiStripDaqInfo::endRun(edm::Run const& run, edm::EventSetup const& eSetup){
  edm::LogInfo ("SiStripDaqInfo") <<"SiStripDaqInfo::EndRun";
  fillDummyStatus();

  if (nFedTotal == 0) {
   edm::LogInfo ("SiStripDaqInfo") <<" SiStripDaqInfo::No FEDs Connected!!!";
   return;
  }

  if (!bookedStatus_) {
    edm::LogError("SiStripDaqInfo") << " SiStripDaqInfo::endRun : MEs missing ";
    return;
  }

  float nFEDConnected = 0.0;

  const FEDNumbering numbering;
  const int siStripFedIdMin = numbering.getSiStripFEDIds().first;
  const int siStripFedIdMax = numbering.getSiStripFEDIds().second; 

  edm::eventsetup::EventSetupRecordKey recordKey(edm::eventsetup::EventSetupRecordKey::TypeTag::findType("RunInfoRcd"));
  if( eSetup.find( recordKey ) != 0) {

    edm::ESHandle<RunInfo> sumFED;
    eSetup.get<RunInfoRcd>().get(sumFED);    
    
    if ( sumFED.isValid() ) {
      vector<int> FedsInIds= sumFED->m_fed_in;   
      for(unsigned int it = 0; it < FedsInIds.size(); ++it) {
	int fedID = FedsInIds[it];     
	
	if(fedID>=siStripFedIdMin &&  fedID<=siStripFedIdMax)  ++nFEDConnected;
      }
      edm::LogInfo ("SiStripDaqInfo") <<" SiStripDaqInfo::Total # of FEDs " << nFedTotal 
                                      << " Connected FEDs " << nFEDConnected;
      if (nFEDConnected > 0) {
	DaqFraction_->Reset();
	DaqFraction_->Fill(nFEDConnected/nFedTotal);
	readSubdetFedFractions(FedsInIds);
      }
    }
  } 
}
//
// -- Read Sub Detector FEDs
//
void SiStripDaqInfo::readFedIds(const edm::ESHandle<SiStripFedCabling>& fedcabling) {
  const vector<uint16_t>& feds = fedCabling_->feds(); 

  nFedTotal = feds.size();
  for(std::vector<unsigned short>::const_iterator ifed = feds.begin(); ifed != feds.end(); ifed++){
    const std::vector<FedChannelConnection> fedChannels = fedCabling_->connections( *ifed );
    for (std::vector<FedChannelConnection>::const_iterator iconn = fedChannels.begin(); iconn < fedChannels.end(); iconn++){
      if (!iconn->isConnected()) continue;
      uint32_t detId = iconn->detId();
      if (detId == 0 || detId == 0xFFFFFFFF)  continue;
      std::string subdet_tag;
      SiStripUtility::getSubDetectorTag(detId,subdet_tag);         
      subDetFedMap[subdet_tag].push_back(*ifed); 
      break;
    }
  }  
}
//
// -- Fill Subdet FEDIds 
//
void SiStripDaqInfo::readSubdetFedFractions(std::vector<int>& fed_ids) {

  const FEDNumbering numbering;
  const int siStripFedIdMin = numbering.getSiStripFEDIds().first;
  const int siStripFedIdMax = numbering.getSiStripFEDIds().second; 

  // initialiase 
  for (std::map<std::string, std::vector<unsigned short> >::const_iterator it = subDetFedMap.begin();
       it != subDetFedMap.end(); it++) {
    std::string name = it->first;
    map<string, SubDetMEs>::iterator iPos = SubDetMEsMap.find(name);
    if (iPos == SubDetMEsMap.end()) continue;
    iPos->second.ConnectedFeds = 0;
  }
  // count sub detector feds

  
  for (std::map<std::string, std::vector<unsigned short> >::const_iterator it = subDetFedMap.begin();
	   it != subDetFedMap.end(); it++) {
    std::string name = it->first;
    std::vector<unsigned short> subdetIds = it->second; 
    map<string, SubDetMEs>::iterator iPos = SubDetMEsMap.find(name);
    if (iPos == SubDetMEsMap.end()) continue;
    iPos->second.ConnectedFeds = 0;
    for (std::vector<unsigned short>::iterator iv = subdetIds.begin();
	 iv != subdetIds.end(); iv++) {
      bool fedid_found = false;
      for(unsigned int it = 0; it < fed_ids.size(); ++it) {
	unsigned short fedID = fed_ids[it];     
        if(fedID < siStripFedIdMin ||  fedID > siStripFedIdMax) continue;
	if ((*iv) == fedID) {
	  fedid_found = true;
          iPos->second.ConnectedFeds++;
	  break;
	}
      }
      if (!fedid_found) findExcludedModule((*iv));   
    }
    int nFedsConnected   = iPos->second.ConnectedFeds;
    int nFedSubDet       = subdetIds.size();
    if (nFedSubDet > 0) {
      iPos->second.DaqFractionME->Reset();
      iPos->second.DaqFractionME->Fill(nFedsConnected*1.0/nFedSubDet);
    }
  }
}
//
// -- find Excluded Modules
//
void SiStripDaqInfo::findExcludedModule(unsigned short fed_id) {
  dqmStore_->cd();
  string mdir = "MechanicalView";
  if (!SiStripUtility::goToDir(dqmStore_, mdir)) return;
  string mechanical_dir = dqmStore_->pwd();
  const std::vector<FedChannelConnection> fedChannels = fedCabling_->connections(fed_id);
  int ichannel = 0;
  string tag = "ExcludedFedChannel";
  string bad_module_folder;
  for (std::vector<FedChannelConnection>::const_iterator iconn = fedChannels.begin(); 
                                                         iconn < fedChannels.end(); iconn++){
    if (!iconn->isConnected()) continue;
    uint32_t detId = iconn->detId();
    if (detId == 0 || detId == 0xFFFFFFFF)  continue;
    
    ichannel++;
    if (ichannel == 1) {
      string subdet_folder ;
      SiStripFolderOrganizer folder_organizer;
      folder_organizer.getSubDetFolder(detId,subdet_folder);
      if (!dqmStore_->dirExists(subdet_folder)) {
	subdet_folder = mechanical_dir + subdet_folder.substr(subdet_folder.find("MechanicalView")+14);
	if (!dqmStore_->dirExists(subdet_folder)) continue;
      }
      bad_module_folder = subdet_folder + "/" + "BadModuleList";
      dqmStore_->setCurrentFolder(bad_module_folder);    
    }
    ostringstream detid_str;
    detid_str << detId;
    string full_path = bad_module_folder + "/" + detid_str.str();
    MonitorElement* me = dqmStore_->get(full_path);
    uint16_t flag = 0; 
    if (me) {
      flag = me->getIntValue();
      me->Reset();
    } else me = dqmStore_->bookInt(detid_str.str());
    SiStripUtility::setBadModuleFlag(tag, flag);
    me->Fill(flag);
  }
  dqmStore_->cd();
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SiStripDaqInfo);
