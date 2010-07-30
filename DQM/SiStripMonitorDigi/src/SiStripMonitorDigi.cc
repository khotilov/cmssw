// -*- c++ -*-
/**\class SiStripMonitorDigi SiStripMonitorDigi.cc DQM/SiStripMonitorDigi/src/SiStripMonitorDigi.cc
 */
// Original Author:  Dorian Kcira
//         Created:  Sat Feb  4 20:49:10 CET 2006
// $Id: SiStripMonitorDigi.cc,v 1.61 2010/05/06 08:23:15 dutta Exp $
#include<fstream>
#include "TNamed.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
#include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
#include "DataFormats/Common/interface/DetSetNew.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/SiStripDigi/interface/SiStripDigi.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "DQM/SiStripCommon/interface/SiStripHistoId.h"
#include "DQM/SiStripMonitorDigi/interface/SiStripMonitorDigi.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DPGAnalysis/SiStripTools/interface/APVCyclePhaseCollection.h"
#include "DPGAnalysis/SiStripTools/interface/EventWithHistory.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDCSStatus.h"

#include "TMath.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "iostream"


//--------------------------------------------------------------------------------------------
SiStripMonitorDigi::SiStripMonitorDigi(const edm::ParameterSet& iConfig) : dqmStore_(edm::Service<DQMStore>().operator->()), conf_(iConfig), show_mechanical_structure_view(true), show_readout_view(false), show_control_view(false), select_all_detectors(true), reset_each_run(false), folder_organizer(), m_cacheID_(0) 
{
  firstEvent = -1;
  eventNb = 0;

  // Detector Partitions
  SubDetPhasePartMap["TIB"]        = "TI";
  SubDetPhasePartMap["TID_side_1"] = "TI";
  SubDetPhasePartMap["TID_side_2"] = "TI";
  SubDetPhasePartMap["TOB"]        = "TO";
  SubDetPhasePartMap["TEC_side_1"] = "TM";
  SubDetPhasePartMap["TEC_side_2"] = "TP";

  // get Digi Producer List   
  digiProducerList = conf_.getParameter<std::vector<edm::InputTag> >("DigiProducersList");

  //get on/off option for every cluster from cfi
  edm::ParameterSet ParametersNumberOfDigis =  conf_.getParameter<edm::ParameterSet>("TH1NumberOfDigis");
  layerswitchnumdigison = ParametersNumberOfDigis.getParameter<bool>("layerswitchon");
  moduleswitchnumdigison = ParametersNumberOfDigis.getParameter<bool>("moduleswitchon");
  
  edm::ParameterSet ParametersNumberOfDigisPerStrip =  conf_.getParameter<edm::ParameterSet>("TH1NumberOfDigisPerStrip");
  moduleswitchnumdigispstripon = ParametersNumberOfDigisPerStrip.getParameter<bool>("moduleswitchon");

  edm::ParameterSet ParametersADCsHottestStrip =  conf_.getParameter<edm::ParameterSet>("TH1ADCsHottestStrip");
  layerswitchadchotteston = ParametersADCsHottestStrip.getParameter<bool>("layerswitchon");
  moduleswitchadchotteston = ParametersADCsHottestStrip.getParameter<bool>("moduleswitchon");
  
  edm::ParameterSet ParametersADCsCoolestStrip =  conf_.getParameter<edm::ParameterSet>("TH1ADCsCoolestStrip");
  layerswitchadccooleston = ParametersADCsCoolestStrip.getParameter<bool>("layerswitchon");
  moduleswitchadccooleston = ParametersADCsCoolestStrip.getParameter<bool>("moduleswitchon");
  
  edm::ParameterSet ParametersDigiADCs =  conf_.getParameter<edm::ParameterSet>("TH1DigiADCs");
  layerswitchdigiadcson = ParametersDigiADCs.getParameter<bool>("layerswitchon");
  moduleswitchdigiadcson = ParametersDigiADCs.getParameter<bool>("moduleswitchon");
   
  edm::ParameterSet ParametersStripOccupancy =  conf_.getParameter<edm::ParameterSet>("TH1StripOccupancy");
  layerswitchstripoccupancyon = ParametersStripOccupancy.getParameter<bool>("layerswitchon");
  moduleswitchstripoccupancyon = ParametersStripOccupancy.getParameter<bool>("moduleswitchon");

  edm::ParameterSet ParametersDigiProf = conf_.getParameter<edm::ParameterSet>("TProfNumberOfDigi");
  layerswitchnumdigisprofon = ParametersDigiProf.getParameter<bool>("layerswitchon");
  edm::ParameterSet ParametersDigiADC = conf_.getParameter<edm::ParameterSet>("TProfDigiADC");
  layerswitchdigiadcprofon = ParametersDigiProf.getParameter<bool>("layerswitchon");

  edm::ParameterSet ParametersTotDigiProf = conf_.getParameter<edm::ParameterSet>("TProfTotalNumberOfDigis");
  subdetswitchtotdigiprofon = ParametersTotDigiProf.getParameter<bool>("subdetswitchon");

  edm::ParameterSet ParametersDigiApvProf = conf_.getParameter<edm::ParameterSet>("TProfDigiApvCycle");
  subdetswitchapvcycleprofon = ParametersDigiApvProf.getParameter<bool>("subdetswitchon");

  edm::ParameterSet ParametersDigiApvTH2 = conf_.getParameter<edm::ParameterSet>("TH2DigiApvCycle");
  subdetswitchapvcycleth2on = ParametersDigiApvTH2.getParameter<bool>("subdetswitchon");

  digitkhistomapon = conf_.getParameter<bool>("TkHistoMap_On"); 
  
  createTrendMEs = conf_.getParameter<bool>("CreateTrendMEs");
  Mod_On_ = conf_.getParameter<bool>("Mod_On");

  // Event History Producer
  historyProducer_ = conf_.getParameter<edm::InputTag>("HistoryProducer");
  // Apv Phase Producer
  apvPhaseProducer_ = conf_.getParameter<edm::InputTag>("ApvPhaseProducer");

  // Create DCS Status
  bool checkDCS    = conf_.getParameter<bool>("UseDCSFiltering");
  if (checkDCS) dcsStatus_ = new SiStripDCSStatus();
  else dcsStatus_ = 0; 
}
//------------------------------------------------------------------------------------------
SiStripMonitorDigi::~SiStripMonitorDigi() { 
  if (dcsStatus_) delete dcsStatus_;
}

//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::beginRun(const edm::Run& run, const edm::EventSetup& es){

  if (show_mechanical_structure_view) {
    unsigned long long cacheID = es.get<SiStripDetCablingRcd>().cacheIdentifier();
    if (m_cacheID_ != cacheID) {
      m_cacheID_ = cacheID;       
      edm::LogInfo("SiStripMonitorDigi") <<"SiStripMonitorDigi::beginRun: " 
					 << " Creating MEs for new Cabling ";     
      createMEs(es);
    } 
  } else if (reset_each_run) {
    edm::LogInfo("SiStripMonitorDigi") <<"SiStripMonitorDigi::beginRun: " 
				       << " Resetting MEs ";        
    for (std::map<uint32_t, ModMEs >::const_iterator idet = DigiMEs.begin() ; idet!=DigiMEs.end() ; idet++) {
      ResetModuleMEs(idet->first);
    }
  }
}

//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::endRun(const edm::Run&, const edm::EventSetup&){
}



//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::beginJob(){
}


//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::createMEs(const edm::EventSetup& es){

  if ( show_mechanical_structure_view ){
    // take from eventSetup the SiStripDetCabling object - here will use SiStripDetControl later on
    es.get<SiStripDetCablingRcd>().get(SiStripDetCabling_);
    
    // get list of active detectors from SiStripDetCabling
    std::vector<uint32_t> activeDets; 
    activeDets.clear(); // just in case
    SiStripDetCabling_->addActiveDetectorsRawIds(activeDets);

    SiStripSubStructure substructure;

    // remove any eventual zero elements - there should be none, but just in case
    for(std::vector<uint32_t>::iterator idets = activeDets.begin(); idets != activeDets.end(); idets++){
      if(*idets == 0) activeDets.erase(idets);
    }
    
    // create SiStripFolderOrganizer
    SiStripFolderOrganizer folder_organizer;

    // Create TkHistoMap for Digi
    if (digitkhistomapon) tkmapdigi = new TkHistoMap("SiStrip/TkHisto","TkHMap_NumberOfDigi",0.0,1);
    
    std::vector<uint32_t> tibDetIds;
    // loop over detectors and book MEs
    edm::LogInfo("SiStripTkDQM|SiStripMonitorDigi")<<"nr. of activeDets:  "<<activeDets.size();
    for(std::vector<uint32_t>::const_iterator detid_iterator = activeDets.begin(); detid_iterator!=activeDets.end(); detid_iterator++){

      uint32_t detid = (*detid_iterator);

      ModMEs local_modmes;
      local_modmes.NumberOfDigis = 0;
      local_modmes.NumberOfDigisPerStrip = 0;
      local_modmes.ADCsHottestStrip = 0;
      local_modmes.ADCsCoolestStrip = 0;
      local_modmes.DigiADCs = 0;
      local_modmes.StripOccupancy = 0;

      if (Mod_On_) {
	// set appropriate folder using SiStripFolderOrganizer
	folder_organizer.setDetectorFolder(detid); // pass the detid to this method
	if (reset_each_run) ResetModuleMEs(detid);
	createModuleMEs(local_modmes, detid);
	// append to DigiMEs
	DigiMEs.insert( std::make_pair(detid, local_modmes));
      }

      // Create Layer Level MEs if they are not created already
      std::pair<std::string,int32_t> det_layer_pair = folder_organizer.GetSubDetAndLayer(detid);
      SiStripHistoId hidmanager;
      std::string label = hidmanager.getSubdetid(detid,false);
      
      // get detids for the layer
      std::map<std::string, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(label);
      if(iLayerME==LayerMEsMap.end()) {
        int32_t lnumber = det_layer_pair.second;
        std::vector<uint32_t> layerDetIds;
        if (det_layer_pair.first == "TIB") {
          substructure.getTIBDetectors(activeDets,layerDetIds,lnumber,0,0,0);
        } else if (det_layer_pair.first == "TOB") {
          substructure.getTOBDetectors(activeDets,layerDetIds,lnumber,0,0);
        } else if (det_layer_pair.first == "TID" && lnumber > 0) {
          substructure.getTIDDetectors(activeDets,layerDetIds,2,abs(lnumber),0,0);
        } else if (det_layer_pair.first == "TID" && lnumber < 0) {
          substructure.getTIDDetectors(activeDets,layerDetIds,1,abs(lnumber),0,0);
        } else if (det_layer_pair.first == "TEC" && lnumber > 0) {
          substructure.getTECDetectors(activeDets,layerDetIds,2,abs(lnumber),0,0,0,0);
        } else if (det_layer_pair.first == "TEC" && lnumber < 0) {
          substructure.getTECDetectors(activeDets,layerDetIds,1,abs(lnumber),0,0,0,0);
        }

        LayerDetMap[label] = layerDetIds;

        // book Layer plots      
	folder_organizer.setLayerFolder(detid,det_layer_pair.second); 
	createLayerMEs(label, layerDetIds.size());
      }
      
      // book sub-detector plots
      std::pair<std::string,std::string> sdet_pair = folder_organizer.getSubDetFolderAndTag(detid);
      if (SubDetMEsMap.find(sdet_pair.second) == SubDetMEsMap.end()){
	dqmStore_->setCurrentFolder(sdet_pair.first);
	createSubDetMEs(sdet_pair.second);        
      }

    }//end of loop over detectors

  }//end of if

}//end of method

//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  // Filter out events if DCS Event if requested
  if (dcsStatus_ && !dcsStatus_->getStatus(iEvent, iSetup)) return;

  runNb   = iEvent.id().run();
  eventNb++;
  float iOrbitSec      = iEvent.orbitNumber()/11223.0;

  digi_detset_handles.clear();
  for(std::vector<edm::InputTag>::iterator itDigiProducerList = digiProducerList.begin();
     itDigiProducerList != digiProducerList.end(); ++itDigiProducerList ) {

    edm::Handle< edm::DetSetVector<SiStripDigi> > digi_handle;

    iEvent.getByLabel((*itDigiProducerList),digi_handle);
    if (digi_handle.isValid()) digi_detset_handles.push_back(digi_handle.product());
  }    

  // initialise # of clusters to zero
  for (std::map<std::string, SubDetMEs>::iterator iSubdet  = SubDetMEsMap.begin();
       iSubdet != SubDetMEsMap.end(); iSubdet++) {
    iSubdet->second.totNDigis = 0;
  }

  for (std::map<std::string, std::vector< uint32_t > >::const_iterator iterLayer = LayerDetMap.begin();
       iterLayer != LayerDetMap.end(); iterLayer++) {
    
    std::string layer_label = iterLayer->first;    
    std::vector< uint32_t > layer_dets = iterLayer->second;
    std::map<std::string, LayerMEs>::iterator iLayerME = LayerMEsMap.find(layer_label);
      
    //get Layer MEs 
    LayerMEs local_layermes;
    if(iLayerME == LayerMEsMap.end()) continue;
    else local_layermes = iLayerME->second; 
    int largest_adc_layer= 0;
    int smallest_adc_layer= 99999;
    int ndigi_layer = 0;
    
    uint16_t iDet = 0;
    std::string subdet_label = ""; 
    // loop over all modules in the layer
    for (std::vector< uint32_t >::const_iterator iterDets = layer_dets.begin() ; 
	 iterDets != layer_dets.end() ; iterDets++) {
      iDet++;
      // detid and type of ME
      uint32_t detid = (*iterDets);
	
      // Get SubDet label once
      if (subdet_label.size() == 0) subdet_label = folder_organizer.getSubDetFolderAndTag(detid).second;

      // DetId and corresponding set of MEs
      std::map<uint32_t, ModMEs >::iterator pos = DigiMEs.find(detid);
      ModMEs local_modmes = pos->second;
	
      // search  digis of detid
      int loc = getDigiSourceIndex(detid); 
      
      int ndigi_det = 0;
      if (loc > -1) ndigi_det = (*(digi_detset_handles[loc]))[detid].size();


      // no digis for this detector module, so fill histogram with 0
      if(Mod_On_ && moduleswitchnumdigison && (local_modmes.NumberOfDigis != NULL))
	(local_modmes.NumberOfDigis)->Fill(ndigi_det); 
      
      if (layerswitchnumdigisprofon) 
	local_layermes.LayerNumberOfDigisProfile->Fill(iDet*1.0,ndigi_det);

      if (digitkhistomapon) tkmapdigi->fill(detid,ndigi_det);

      if (ndigi_det == 0) continue; // no digis for this detid => jump to next step of loop
      
      const edm::DetSet<SiStripDigi> & digi_detset = (*(digi_detset_handles[loc]))[detid]; 
     
      ndigi_layer += ndigi_det;	
      // ADCs
      int largest_adc=(digi_detset.data.begin())->adc();
      int smallest_adc=(digi_detset.data.begin())->adc();
      

      // Check if these parameters are really needed
      float det_occupancy = 0.0;
      
      for(edm::DetSet<SiStripDigi>::const_iterator digiIter = digi_detset.data.begin(); 
	  digiIter!= digi_detset.data.end(); digiIter++ ){
	
	int this_adc = digiIter->adc();
	if (this_adc > 0.0) det_occupancy++;
	
	if(this_adc>largest_adc) largest_adc  = this_adc; 
	if(this_adc<smallest_adc) smallest_adc  = this_adc; 

	if(Mod_On_ && moduleswitchnumdigispstripon && (local_modmes.NumberOfDigisPerStrip != NULL) && (this_adc > 0.0) )
          (local_modmes.NumberOfDigisPerStrip)->Fill(digiIter->strip());

	if(Mod_On_ && moduleswitchdigiadcson && (local_modmes.DigiADCs != NULL) )
	  (local_modmes.DigiADCs)->Fill(static_cast<float>(this_adc));
	
	//Fill #ADCs for this digi at layer level
	if(layerswitchdigiadcson) {
	  fillME(local_layermes.LayerDigiADCs , this_adc);
	  if (createTrendMEs) fillTrend(local_layermes.LayerDigiADCsTrend, this_adc, iOrbitSec);
	}
	
	if (layerswitchdigiadcprofon) 
	  local_layermes.LayerDigiADCProfile->Fill(iDet*1.0,this_adc);
	
      }//end of loop over digis in this det
      
      // Occupancy
      short nstrips = SiStripDetCabling_->nApvPairs(detid) * 2 * 128;
      if (nstrips > 0 && det_occupancy > 0 ) {
	det_occupancy = det_occupancy/nstrips;
	if (Mod_On_ && moduleswitchstripoccupancyon && (local_modmes.StripOccupancy != NULL))
	  (local_modmes.StripOccupancy)->Fill(det_occupancy);
	if (layerswitchstripoccupancyon) {
	  fillME(local_layermes.LayerStripOccupancy, det_occupancy);
	  if (createTrendMEs) fillTrend(local_layermes.LayerStripOccupancyTrend, det_occupancy, iOrbitSec);
	}
      }
      
      if  (largest_adc > largest_adc_layer) largest_adc_layer = largest_adc;
      if  (smallest_adc < smallest_adc_layer) smallest_adc_layer = smallest_adc;
      
      // nr. of adcs for hottest strip
      if( Mod_On_ && moduleswitchadchotteston && (local_modmes.ADCsHottestStrip != NULL)) 
	(local_modmes.ADCsHottestStrip)->Fill(static_cast<float>(largest_adc));
      
      // nr. of adcs for coolest strip	
      if(Mod_On_ && moduleswitchadccooleston && (local_modmes.ADCsCoolestStrip != NULL)) 
	(local_modmes.ADCsCoolestStrip)->Fill(static_cast<float>(smallest_adc));
      
    }//end of loop over DetIds
    
    if(layerswitchnumdigison) {
      fillME(local_layermes.LayerNumberOfDigis,ndigi_layer);
      if (createTrendMEs) fillTrend(local_layermes.LayerNumberOfDigisTrend, ndigi_layer, iOrbitSec);
    }
    if(layerswitchadchotteston) {
      fillME(local_layermes.LayerADCsHottestStrip,largest_adc_layer);
      if (createTrendMEs) fillTrend(local_layermes.LayerADCsHottestStripTrend, largest_adc_layer, iOrbitSec);
    }
    if(layerswitchadccooleston) {
      fillME(local_layermes.LayerADCsCoolestStrip ,smallest_adc_layer);
      if (createTrendMEs) fillTrend(local_layermes.LayerADCsCoolestStripTrend, smallest_adc_layer, iOrbitSec);
    }
    std::map<std::string, SubDetMEs>::iterator iSubdet  = SubDetMEsMap.find(subdet_label);
    if(iSubdet != SubDetMEsMap.end()) iSubdet->second.totNDigis += ndigi_layer;
  }

  // get EventHistory 
  edm::Handle<EventWithHistory> event_history;
  iEvent.getByLabel(historyProducer_,event_history);  

  // get Phase of APV
  edm::Handle<APVCyclePhaseCollection> apv_phase_collection;
  iEvent.getByLabel(apvPhaseProducer_,apv_phase_collection);

  if (event_history.isValid() 
      && !event_history.failedToGet()
      && apv_phase_collection.isValid() 
      && !apv_phase_collection.failedToGet()) {

    
    long long tbx = event_history->absoluteBX();
    float iOrbitSec = iEvent.orbitNumber()/11223.0;

    for (std::map<std::string, SubDetMEs>::iterator it = SubDetMEsMap.begin();
	 it != SubDetMEsMap.end(); it++) {

      SubDetMEs subdetmes;
      std::string subdet = it->first;
      subdetmes = it->second;
 
      int the_phase = APVCyclePhaseCollection::invalid;
      long long tbx_corr = tbx;

      if (SubDetPhasePartMap.find(subdet) != SubDetPhasePartMap.end()) the_phase = apv_phase_collection->getPhase(SubDetPhasePartMap[subdet]);
      if(the_phase==APVCyclePhaseCollection::nopartition ||
         the_phase==APVCyclePhaseCollection::multiphase ||
         the_phase==APVCyclePhaseCollection::invalid) the_phase=30;
      tbx_corr  -= the_phase;
      
      if (subdetswitchtotdigiprofon)subdetmes.SubDetTotDigiProf->Fill(iOrbitSec,subdetmes.totNDigis);  
      if (subdetswitchapvcycleprofon)subdetmes.SubDetDigiApvProf->Fill(tbx_corr%70,subdetmes.totNDigis);
      if (subdetswitchapvcycleth2on) subdetmes.SubDetDigiApvTH2->Fill(tbx_corr%70,subdetmes.totNDigis);
    }
  }
}//end of method analyze
//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::endJob(void){
  bool outputMEsInRootFile = conf_.getParameter<bool>("OutputMEsInRootFile");
  std::string outputFileName = conf_.getParameter<std::string>("OutputFileName");

  // save histograms in a file
  if(outputMEsInRootFile)     dqmStore_->save(outputFileName);
  
}//end of method
//--------------------------------------------------------------------------------------------
void SiStripMonitorDigi::ResetModuleMEs(uint32_t idet){
  std::map<uint32_t, ModMEs >::iterator pos = DigiMEs.find(idet);
  ModMEs mod_me = pos->second;

  if(Mod_On_ && moduleswitchnumdigison) mod_me.NumberOfDigis->Reset();
  if(Mod_On_ && moduleswitchnumdigispstripon) mod_me.NumberOfDigisPerStrip ->Reset();
  if(Mod_On_ && moduleswitchadchotteston) mod_me.ADCsHottestStrip->Reset();
  if(Mod_On_ && moduleswitchadccooleston) mod_me.ADCsCoolestStrip->Reset();
  if(Mod_On_ && moduleswitchdigiadcson) mod_me.DigiADCs->Reset();
  if(Mod_On_ && moduleswitchstripoccupancyon) mod_me.StripOccupancy->Reset();

}
//------------------------------------------------------------------------------------------
MonitorElement* SiStripMonitorDigi::bookMETrend(const char* ParameterSetLabel, const char* HistoName)
{
  edm::ParameterSet ParametersTrend =  conf_.getParameter<edm::ParameterSet>("Trending");
  MonitorElement* me = dqmStore_->bookProfile(HistoName,HistoName,
					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      // 					      0,
					      ParametersTrend.getParameter<double>("xmin"),
					      ParametersTrend.getParameter<double>("xmax"),
					      // 					      ParametersTrend.getParameter<int32_t>("Nbins"),
					      100, //that parameter should not be there !?
					      ParametersTrend.getParameter<double>("ymin"),
					      ParametersTrend.getParameter<double>("ymax"),
					      "" );
  if(!me) return me;

  me->setAxisTitle("Event Time in Seconds",1);
  if (me->kind() == MonitorElement::DQM_KIND_TPROFILE) me->getTH1()->SetBit(TH1::kCanRebin);
  return me;
}

//------------------------------------------------------------------------------------------
MonitorElement* SiStripMonitorDigi::bookME1D(const char* ParameterSetLabel, const char* HistoName)
{
  edm::ParameterSet Parameters =  conf_.getParameter<edm::ParameterSet>(ParameterSetLabel);
  return dqmStore_->book1D(HistoName,HistoName,
			   Parameters.getParameter<int32_t>("Nbinx"),
			   Parameters.getParameter<double>("xmin"),
			   Parameters.getParameter<double>("xmax")
			   );
}

//--------------------------------------------------------------------------------
void SiStripMonitorDigi::fillTrend(MonitorElement* me ,float value, float timeinorbit)
{
  if(!me) return;
  me->Fill(timeinorbit,value);
}

//
// -- Create Module Level MEs
//
void SiStripMonitorDigi::createModuleMEs(ModMEs& mod_single, uint32_t detid) {

  // use SistripHistoId for producing histogram id (and title)
  SiStripHistoId hidmanager;
  std::string hid;
  
  //nr. of digis per module
  if(moduleswitchnumdigison) {
    hid = hidmanager.createHistoId("NumberOfDigis","det",detid);
    mod_single.NumberOfDigis = dqmStore_->book1D(hid, hid, 21, -0.5, 20.5);
    dqmStore_->tag(mod_single.NumberOfDigis, detid);
    mod_single.NumberOfDigis->setAxisTitle("number of digis in one detector module");
    mod_single.NumberOfDigis->getTH1()->StatOverflows(kTRUE);  // over/underflows in Mean calculation
  }
  
  //nr. of digis per strip in module
  if(moduleswitchnumdigispstripon){
    hid = hidmanager.createHistoId("NumberOfDigisPerStrip","det",detid);
    short nstrips = SiStripDetCabling_->nApvPairs(detid) * 2 * 128; 
    mod_single.NumberOfDigisPerStrip = dqmStore_->book1D(hid, hid, nstrips, -0.5, nstrips+0.5);
    dqmStore_->tag(mod_single.NumberOfDigisPerStrip, detid);
    mod_single.NumberOfDigisPerStrip->setAxisTitle("number of (digis > 0) per strip");
    mod_single.NumberOfDigisPerStrip->getTH1()->StatOverflows(kTRUE);  // over/underflows in Mean calculation
  }
  //#ADCs for hottest strip
  if(moduleswitchadchotteston) {
    hid = hidmanager.createHistoId("ADCsHottestStrip","det",detid);
    mod_single.ADCsHottestStrip = bookME1D("TH1ADCsHottestStrip", hid.c_str());
    dqmStore_->tag(mod_single.ADCsHottestStrip, detid); // 6 APVs -> 768 strips
    mod_single.ADCsHottestStrip->setAxisTitle("number of ADCs for hottest strip");
  }
  
  //#ADCs for coolest strip
  if(moduleswitchadccooleston) {
    hid = hidmanager.createHistoId("ADCsCoolestStrip","det",detid);
    mod_single.ADCsCoolestStrip = bookME1D("TH1ADCsCoolestStrip", hid.c_str());
    dqmStore_->tag(mod_single.ADCsCoolestStrip, detid);
    mod_single.ADCsCoolestStrip->setAxisTitle("number of ADCs for coolest strip");
  }
  
  //#ADCs for each digi
  if(moduleswitchdigiadcson) {
    hid = hidmanager.createHistoId("DigiADCs","det",detid);
    mod_single.DigiADCs = bookME1D("TH1DigiADCs", hid.c_str());
    dqmStore_->tag(mod_single.DigiADCs, detid);
    mod_single.DigiADCs->setAxisTitle("number of ADCs for each digi");
  }
  
  //Strip occupancy
  if(moduleswitchstripoccupancyon) {
    hid = hidmanager.createHistoId("StripOccupancy","det",detid);
    mod_single.StripOccupancy = bookME1D("TH1StripOccupancy", hid.c_str());
    dqmStore_->tag(mod_single.StripOccupancy, detid);
    mod_single.StripOccupancy->setAxisTitle("strip occupancy");
  }
  
}
  
//
// -- Create Module Level MEs
//  

void SiStripMonitorDigi::createLayerMEs(std::string label, int ndets) {

  std::map<std::string, LayerMEs>::iterator iLayerME  = LayerMEsMap.find(label);
  if(iLayerME==LayerMEsMap.end()){
    SiStripHistoId hidmanager;
    LayerMEs layerMEs; 
    layerMEs.LayerNumberOfDigis = 0;
    layerMEs.LayerNumberOfDigisTrend = 0;
    layerMEs.LayerADCsHottestStrip = 0;
    layerMEs.LayerADCsHottestStripTrend = 0;
    layerMEs.LayerADCsCoolestStrip = 0;
    layerMEs.LayerADCsCoolestStripTrend = 0;
    layerMEs.LayerDigiADCs = 0;
    layerMEs.LayerDigiADCsTrend = 0;
    layerMEs.LayerStripOccupancy = 0;
    layerMEs.LayerStripOccupancyTrend = 0;
    layerMEs.LayerNumberOfDigisProfile = 0;
    layerMEs.LayerDigiADCProfile = 0;
    

    //#Digis
    if(layerswitchnumdigison) {
      layerMEs.LayerNumberOfDigis=bookME1D("TH1NumberOfDigis", hidmanager.createHistoLayer("Summary_TotalNumberOfDigis","layer",label,"").c_str()); 
      if (createTrendMEs) layerMEs.LayerNumberOfDigisTrend=bookMETrend("TH1NumberOfDigis", hidmanager.createHistoLayer("Trend_NumberOfDigis","layer",label,"").c_str()); 
    }

    //#ADCs for hottest strip
    if(layerswitchadchotteston) {
      layerMEs.LayerADCsHottestStrip=bookME1D("TH1ADCsHottestStrip", hidmanager.createHistoLayer("Summary_ADCsHottestStrip","layer",label,"").c_str()); 
      if (createTrendMEs) layerMEs.LayerADCsHottestStripTrend=bookMETrend("TH1ADCsHottestStrip", hidmanager.createHistoLayer("Trend_ADCsHottestStrip","layer",label,"").c_str()); 
    }

    //#ADCs for coolest strip
    if(layerswitchadccooleston) {
      layerMEs.LayerADCsCoolestStrip=bookME1D("TH1ADCsCoolestStrip", hidmanager.createHistoLayer("Summary_ADCsCoolestStrip","layer",label,"").c_str());
      if (createTrendMEs) layerMEs.LayerADCsCoolestStripTrend=bookMETrend("TH1ADCsCoolestStrip", hidmanager.createHistoLayer("Trend_ADCsCoolestStrip","layer",label,"").c_str());
    }

    //#ADCs for each digi
    if(layerswitchdigiadcson) {
      layerMEs.LayerDigiADCs=bookME1D("TH1DigiADCs", hidmanager.createHistoLayer("Summary_DigiADCs","layer",label,"").c_str());
      if (createTrendMEs) layerMEs.LayerDigiADCsTrend=bookMETrend("TH1DigiADCs", hidmanager.createHistoLayer("Trend_DigiADCs","layer",label,"").c_str());
    }

    //Strip Occupancy
    if(layerswitchstripoccupancyon) {
      layerMEs.LayerStripOccupancy=bookME1D("TH1StripOccupancy", hidmanager.createHistoLayer("Summary_StripOccupancy","layer",label,"").c_str());  
      if (createTrendMEs) layerMEs.LayerStripOccupancyTrend=bookMETrend("TH1StripOccupancy", hidmanager.createHistoLayer("Trend_StripOccupancy","layer",label,"").c_str());  
      
    }
    // # of Digis 
    if(layerswitchnumdigisprofon) {
      std::string hid = hidmanager.createHistoLayer("NumberOfDigiProfile","layer",label,"");
      layerMEs.LayerNumberOfDigisProfile = dqmStore_->bookProfile(hid, hid, ndets, 0.5, ndets+0.5,21, -0.5, 200.5);
    }

    // # of Digis 
    if(layerswitchdigiadcprofon) {
      std::string hid = hidmanager.createHistoLayer("DigiADCProfile","layer",label,"");      
      layerMEs.LayerDigiADCProfile = dqmStore_->bookProfile(hid, hid, ndets, 0.5, ndets+0.5, 64, -0.5, 255.5);      
    }

    LayerMEsMap[label]=layerMEs;
  }
}
//
// -- Create SubDetector MEs
//
void SiStripMonitorDigi::createSubDetMEs(std::string label) {

  SubDetMEs subdetMEs; 
  subdetMEs.totNDigis         = 0;
  subdetMEs.SubDetTotDigiProf = 0;
  subdetMEs.SubDetDigiApvProf = 0;
  subdetMEs.SubDetDigiApvTH2  = 0;
  
  std::string HistoName;
  
  // Total Number of Digi - Profile
  if(subdetswitchtotdigiprofon){
    edm::ParameterSet Parameters =  conf_.getParameter<edm::ParameterSet>("TProfTotalNumberOfDigis");
    HistoName = "TotalNumberOfDigiProfile__" + label;
    subdetMEs.SubDetTotDigiProf=dqmStore_->bookProfile(HistoName,HistoName,
						       Parameters.getParameter<int32_t>("Nbins"),
						       Parameters.getParameter<double>("xmin"),
						       Parameters.getParameter<double>("xmax"),
						       100, //that parameter should not be there !?
						       Parameters.getParameter<double>("ymin"),
						       Parameters.getParameter<double>("ymax"),
						       "" );
    subdetMEs.SubDetTotDigiProf->setAxisTitle("Event Time in Seconds",1);
    if (subdetMEs.SubDetTotDigiProf->kind() == MonitorElement::DQM_KIND_TPROFILE) subdetMEs.SubDetTotDigiProf->getTH1()->SetBit(TH1::kCanRebin);
  }
  
  // Number of Digi vs Bx - Profile
  if(subdetswitchapvcycleprofon){
    edm::ParameterSet Parameters =  conf_.getParameter<edm::ParameterSet>("TProfDigiApvCycle");
    HistoName = "Digi_vs_ApvCycle__" + label;
    subdetMEs.SubDetDigiApvProf=dqmStore_->bookProfile(HistoName,HistoName,
						       Parameters.getParameter<int32_t>("Nbins"),
						       Parameters.getParameter<double>("xmin"),
						       Parameters.getParameter<double>("xmax"),
						       200, //that parameter should not be there !?
						       Parameters.getParameter<double>("ymin"),
						       Parameters.getParameter<double>("ymax"),
						       "" );
    subdetMEs.SubDetDigiApvProf->setAxisTitle("ApvCycle (Corrected Absolute Bx % 70)",1);
  }
  
  // Number of Digi vs Bx - TH2
  if(subdetswitchapvcycleth2on){
    edm::ParameterSet Parameters =  conf_.getParameter<edm::ParameterSet>("TH2DigiApvCycle");
    dqmStore_->setCurrentFolder("SiStrip/MechanicalView/"+label);
    HistoName = "Digi_vs_ApvCycle_2D__" + label;
    // Adjusting the scale for 2D histogram
    double h2ymax = 9999.0;
    double yfact = Parameters.getParameter<double>("yfactor");
    if(label.find("TIB") != std::string::npos) h2ymax = (6984.*256.)*yfact;
    else if (label.find("TID") != std::string::npos) h2ymax = (2208.*256.)*yfact;
    else if (label.find("TOB") != std::string::npos) h2ymax = (12906.*256.)*yfact;
    else if (label.find("TEC") != std::string::npos) h2ymax = (7552.*2.*256.)*yfact;
    subdetMEs.SubDetDigiApvTH2=dqmStore_->book2D(HistoName,HistoName,
						 Parameters.getParameter<int32_t>("Nbins"),
						 Parameters.getParameter<double>("xmin"),
						 Parameters.getParameter<double>("xmax"),
						 Parameters.getParameter<int32_t>("Nbinsy"), //it was 100 that parameter should not be there !?
						 Parameters.getParameter<double>("ymin"),
						 h2ymax);
    subdetMEs.SubDetDigiApvTH2->setAxisTitle("absolute Bx mod(70)",1);
  }
  SubDetMEsMap[label]=subdetMEs;
}
//
// -- Get DetSet vector for a given Detector
//
int SiStripMonitorDigi::getDigiSourceIndex(uint32_t id) {
  int location = -1;
  for (unsigned int ival = 0; ival <  digi_detset_handles.size(); ++ival){
    edm::DetSetVector<SiStripDigi>::const_iterator isearch = digi_detset_handles[ival]->find(id); 
    if(isearch != digi_detset_handles[ival]->end()) {
      location = ival;
      break;
    }
  }
  return location;
}
//define this as a plug-in
DEFINE_FWK_MODULE(SiStripMonitorDigi);
