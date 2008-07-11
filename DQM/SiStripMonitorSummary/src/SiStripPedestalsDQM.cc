#include "DQM/SiStripMonitorSummary/interface/SiStripPedestalsDQM.h"
#include "DQMServices/Core/interface/MonitorElement.h"

// -----
SiStripPedestalsDQM::SiStripPedestalsDQM(const edm::EventSetup & eSetup,
                                         edm::ParameterSet const& hPSet,
                                         edm::ParameterSet const& fPSet):SiStripBaseCondObjDQM(eSetup, hPSet, fPSet){

}
// -----



// -----
SiStripPedestalsDQM::~SiStripPedestalsDQM(){}
// -----


// -----
void SiStripPedestalsDQM::getActiveDetIds(const edm::EventSetup & eSetup){
  
  getConditionObject(eSetup);
  pedestalHandle_->getDetIds(activeDetIds);
  selectModules(activeDetIds);

}
// -----


// -----
void SiStripPedestalsDQM::fillModMEs(const std::vector<uint32_t> & selectedDetIds){
   
  ModMEs CondObj_ME;
  
  
  for(std::vector<uint32_t>::const_iterator detIter_ = selectedDetIds.begin();
                                            detIter_!= selectedDetIds.end();detIter_++){
      
    fillMEsForDet(CondObj_ME,*detIter_);
      
  }
}    
// -----




// -----
void SiStripPedestalsDQM::fillMEsForDet(ModMEs selModME_, uint32_t selDetId_){
  
  getModMEs(selModME_,selDetId_);
  
  SiStripPedestals::Range pedRange = pedestalHandle_->getRange(selDetId_);
  int nStrip =  reader->getNumberOfApvsAndStripLength(selDetId_).first*128;
  
  for( int istrip=0;istrip<nStrip;++istrip){
    try{      
      if( CondObj_fillId_ =="onlyProfile" || CondObj_fillId_ =="ProfileAndCumul"){
        selModME_.ProfileDistr->Fill(istrip+1,pedestalHandle_->getPed(istrip,pedRange));
      }
    } 
    catch(cms::Exception& e){
      edm::LogError("SiStripPedestalsDQM")          
	<< "[SiStripPedestalsDQM::fillMEsForDet] cms::Exception accessing pedestalHandle_->getPed(istrip,pedRange) for strip "  
	<< istrip 
	<< " and detid " 
	<< selDetId_  
	<< " :  " 
	<< e.what() ;
    }
  }// istrip
  
}    
// -----



// -----
void SiStripPedestalsDQM::fillSummaryMEs(const std::vector<uint32_t> & selectedDetIds){
   
  for(std::vector<uint32_t>::const_iterator detIter_ = selectedDetIds.begin();
                                            detIter_!= selectedDetIds.end();detIter_++){
    fillMEsForLayer(SummaryMEsMap_, *detIter_);
  } 
}           
  
// -----



// -----
void SiStripPedestalsDQM::fillMEsForLayer( std::map<uint32_t, ModMEs> selMEsMap_, uint32_t selDetId_){

  // ----
  int subdetectorId_ = ((selDetId_>>25)&0x7);
    
  if( subdetectorId_<3 || subdetectorId_>6 ){ 
    edm::LogError("SiStripPedestalsDQM")
       << "[SiStripPedestalsDQM::fillMEsForLayer] WRONG INPUT : no such subdetector type : "
       << subdetectorId_ << " no folder set!" 
       << std::endl;
    return;
  }
  // ----
     
  std::map<uint32_t, ModMEs>::iterator selMEsMapIter_ = selMEsMap_.find(getLayerNameAndId(selDetId_).second);
  ModMEs selME_;
  selME_ =selMEsMapIter_->second;
  getSummaryMEs(selME_,selDetId_);
    
  SiStripPedestals::Range pedRange = pedestalHandle_->getRange(selDetId_);
  
  int nStrip =  reader->getNumberOfApvsAndStripLength(selDetId_).first*128;
  
  SiStripHistoId hidmanager;
  
  // --> profile summary    
  
  std::string hSummaryOfProfile_description;
  hSummaryOfProfile_description  = hPSet_.getParameter<std::string>("SummaryOfProfile_description");
  
  std::string hSummaryOfProfile_name; 
  
  
  hSummaryOfProfile_name = hidmanager.createHistoLayer(hSummaryOfProfile_description, 
						       "layer", 
						       getLayerNameAndId(selDetId_).first, 
						       "") ;
  
  
  for( int istrip=0;istrip<nStrip;++istrip){
    
    try{ 
      if( CondObj_fillId_ =="onlyProfile" || CondObj_fillId_ =="ProfileAndCumul"){
	selME_.SummaryOfProfileDistr->Fill(istrip+1,pedestalHandle_->getPed(istrip,pedRange));
      }
    } 
    catch(cms::Exception& e){
      edm::LogError("SiStripPedestalsDQM")          
	<< "[SiStripPedestalsDQM::fillMEsForLayer] cms::Exception accessing pedestalHandle_->getPed(istrip,pedRange) for strip "  
	<< istrip 
	<< " and detid " 
	<< selDetId_  
	<< " :  " 
	<< e.what() ;
    }
  }// istrip	
  // --> cumul summary  
  
  std::string hSummary_description;
  hSummary_description  = hPSet_.getParameter<std::string>("Summary_description");
  
  std::string hSummary_name; 
  hSummary_name = hidmanager.createHistoLayer(hSummary_description, 
					      "layer", 
					      getLayerNameAndId(selDetId_).first, 
					      "") ;
  float meanPedestal=0;
  
  for( int istrip=0;istrip<nStrip;++istrip){
    
    try{
      meanPedestal = meanPedestal + pedestalHandle_->getPed(istrip,pedRange);
    }
    catch(cms::Exception& e){
      edm::LogError("SiStripNoisesDQM")          
	<< "[SiStripNoisesDQM::fillMEsForLayer] cms::Exception accessing pedestalHandle_->getPed(istrip,pedRange) for strip "  
	<< istrip 
	<< "and detid " 
	<< selDetId_  
	<< " :  " 
	<< e.what() ;      
    }
    
  }//istrip
  
  
  meanPedestal = meanPedestal/nStrip;
  
  
  // -----
  // get detIds belonging to same layer to fill X-axis with detId-number
  
  std::vector<uint32_t> sameLayerDetIds_;
  
  uint32_t subselDetId_ =  ((selDetId_>>25)&0x7);
  SiStripSubStructure substructure_;
  
  sameLayerDetIds_.clear();

  sameLayerDetIds_=GetSameLayerDetId(activeDetIds,subselDetId_);
  
  
  unsigned int iBin=0;
  for(unsigned int i=0;i<sameLayerDetIds_.size();i++){
    if(sameLayerDetIds_[i]==selDetId_){iBin=i+1;}
  }   
  
  selME_.SummaryDistr->Fill(iBin,meanPedestal);
  
  
  
}  
// -----

