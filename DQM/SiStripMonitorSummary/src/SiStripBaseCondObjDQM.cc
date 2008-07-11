#include "DQM/SiStripMonitorSummary/interface/SiStripBaseCondObjDQM.h"

// -----
SiStripBaseCondObjDQM::SiStripBaseCondObjDQM(const edm::EventSetup & eSetup,
                                             edm::ParameterSet const& hPSet,
                                             edm::ParameterSet const& fPSet ):
  eSetup_(eSetup),
  hPSet_(hPSet),
  fPSet_(fPSet),
  cacheID_memory(0),
  dqmStore_(edm::Service<DQMStore>().operator->()){

  reader = new SiStripDetInfoFileReader(edm::FileInPath(std::string("CalibTracker/SiStripCommon/data/SiStripDetInfo.dat") ).fullPath());
  
  Mod_On_                 = fPSet_.getParameter<bool>("Mod_On");
  SummaryOnLayerLevel_On_ = fPSet_.getParameter<bool>("SummaryOnLayerLevel_On");
  SummaryOnStringLevel_On_ = fPSet_.getParameter<bool>("SummaryOnStringLevel_On");
  
  CondObj_fillId_    = hPSet_.getParameter<std::string>("CondObj_fillId");
  CondObj_name_      = hPSet_.getParameter<std::string>("CondObj_name");

  FillSummaryAtLayerLevel = hPSet_.getParameter<bool>("FillSummaryAtLayerLevel");
  FillSummaryProfileAtLayerLevel = hPSet_.getParameter<bool>("FillSummaryProfileAtLayerLevel");
  FillCumulativeSummaryAtLayerLevel = hPSet_.getParameter<bool>("FillCumulativeSummaryAtLayerLevel");

  //Warning message from wrong input:
  if(SummaryOnLayerLevel_On_ && SummaryOnStringLevel_On_){
    edm::LogWarning("SiStripBaseCondObjDQM") 
       << "[SiStripBaseCondObjDQM::SiStripBaseCondObjDQMs] PLEASE CHECK : String and layer level options can not be activated together"
       << std::endl; 
  }
 

}
// -----



// -----
void SiStripBaseCondObjDQM::analysis(const edm::EventSetup & eSetup_){
 
  cacheID_current=  getCache(eSetup_);
  
  if (cacheID_memory == cacheID_current) return;
  
  getConditionObject(eSetup_);

  getActiveDetIds(eSetup_);
  
  if(Mod_On_ )                                            { fillModMEs(activeDetIds); }
  if(SummaryOnLayerLevel_On_ || SummaryOnStringLevel_On_ ){ fillSummaryMEs(activeDetIds);}

}
// -----



// -----
void SiStripBaseCondObjDQM::analysisOnDemand(const edm::EventSetup & eSetup_, 
                                            std::string requestedSubDetector, 
                                            uint32_t requestedSide, 
					    uint32_t requestedLayer){
  
  getConditionObject(eSetup_);
  getActiveDetIds(eSetup_);
   
  std::vector<uint32_t> requestedDetIds_;
  requestedDetIds_.clear();
  
  SiStripSubStructure substructure_;
  
  if(requestedSubDetector=="TIB"){ 
      substructure_.getTIBDetectors( activeDetIds, requestedDetIds_, requestedLayer,0,0,0);
  }
  else if(requestedSubDetector=="TID"){ 
      substructure_.getTIDDetectors( activeDetIds, requestedDetIds_, requestedSide,requestedLayer,0,0);
  }
  else if(requestedSubDetector=="TOB"){  
      substructure_.getTOBDetectors( activeDetIds, requestedDetIds_, requestedLayer,0,0);
  }
  else if(requestedSubDetector=="TEC"){  
      substructure_.getTECDetectors( activeDetIds, requestedDetIds_, requestedSide,requestedLayer,0,0,0,0);
  } 

  analysisOnDemand(eSetup_,requestedDetIds_);
 
}
// -----


// -----
void SiStripBaseCondObjDQM::analysisOnDemand(const edm::EventSetup & eSetup_, uint32_t  detIdOnDemand){
 
  unsigned long long cacheID_current=  getCache(eSetup_);
  
  if (cacheID_memory == cacheID_current) return;
  
  getConditionObject(eSetup_);
  
  std::vector<uint32_t> vdetIdsOnDemand_;
  vdetIdsOnDemand_.push_back(detIdOnDemand); // fillModMEs needs a vector 

  fillModMEs(vdetIdsOnDemand_); 
  
}
// -----

// -----
void SiStripBaseCondObjDQM::analysisOnDemand(const edm::EventSetup & eSetup_, std::vector<uint32_t>  detIdsOnDemand){
 
  unsigned long long cacheID_current=  getCache(eSetup_);
  
  if (cacheID_memory == cacheID_current) return;
  
  getConditionObject(eSetup_);
  

  fillSummaryMEs(detIdsOnDemand); 
  
}
// -----

// -----
std::vector<uint32_t> SiStripBaseCondObjDQM::getCabledModules() {     
 
  std::vector<uint32_t> cabledDetIds_;  
  eSetup_.get<SiStripDetCablingRcd>().get(detCablingHandle_);
  detCablingHandle_->addActiveDetectorsRawIds(cabledDetIds_);

  return cabledDetIds_;  

}
// -----



// -----
void SiStripBaseCondObjDQM::selectModules(std::vector<uint32_t> & detIds_){
   
  ModulesToBeExcluded_     = fPSet_.getParameter< std::vector<unsigned int> >("ModulesToBeExcluded");
  ModulesToBeIncluded_     = fPSet_.getParameter< std::vector<unsigned int> >("ModulesToBeIncluded");
  SubDetectorsToBeExcluded_= fPSet_.getParameter< std::vector<std::string> >("SubDetectorsToBeExcluded");  
 
  if(fPSet_.getParameter<bool>("restrictModules") 
     && ModulesToBeExcluded_.size()==0 
     && ModulesToBeIncluded_.size()==0 ){
    edm::LogWarning("SiStripBaseCondObjDQM") 
       << "[SiStripBaseCondObjDQM::selectModules] PLEASE CHECK : no modules to be exclude/included in your cfg"
       << std::endl; 
  }

 
 
  
  // --> detIds to start with

  if( fPSet_.getParameter<bool>("restrictModules")){
      
    if( ModulesToBeIncluded_.size()>0 ){
     std::vector<uint32_t> tmp;
     tmp.clear();
     set_intersection( detIds_.begin(), detIds_.end(),  
                       ModulesToBeIncluded_.begin(), ModulesToBeIncluded_.end(),
		       inserter(tmp,tmp.begin()));
     swap(detIds_,tmp);
   }
    
  }

  

  // -----
  // *** exclude modules ***
  
  if( fPSet_.getParameter<bool>("restrictModules") ){
    
    std::sort(detIds_.begin(),detIds_.end());

    for( std::vector<uint32_t>::const_iterator modIter_  = ModulesToBeExcluded_.begin(); 
                                               modIter_ != ModulesToBeExcluded_.end(); modIter_++){
      
      std::vector<uint32_t>::iterator detIter_=std::lower_bound(detIds_.begin(),detIds_.end(),*modIter_);
      detIds_.erase(detIter_);
      detIter_--;
     
    }
  
  }
  // *** exclude modules ***
  // -----
   

  // -----
  // *** restrict to a particular subdetector ***
   
  if( *(SubDetectorsToBeExcluded_.begin()) !="none" ){
    
    std::vector<uint32_t> tmp;

    SiStripSubStructure substructure_;
    
    for( std::vector<std::string>::const_iterator modIter_  = SubDetectorsToBeExcluded_.begin(); 
                                                 modIter_ != SubDetectorsToBeExcluded_.end(); modIter_++){
      tmp.clear();

      if (*modIter_=="TIB")     { substructure_.getTIBDetectors(detIds_, tmp, 0,0,0,0);}
      else if (*modIter_=="TOB") { substructure_.getTOBDetectors(detIds_, tmp, 0,0,0);}
      else if (*modIter_=="TID") { substructure_.getTIDDetectors(detIds_, tmp, 0,0,0,0);}
      else if (*modIter_=="TEC") { substructure_.getTECDetectors(detIds_, tmp, 0,0,0,0,0,0);}
      else {
        edm::LogWarning("SiStripBaseCondObjDQM") 
       << "[SiStripBaseCondObjDQM::selectModules] PLEASE CHECK : no correct (name) subdetector to be excluded in your cfg"
       << std::endl; 
      }

      std::vector<uint32_t>::iterator iterBegin_=std::lower_bound(detIds_.begin(),
	                                                          detIds_.end(),
								  *min_element(tmp.begin(), tmp.end()));
								  
      std::vector<uint32_t>::iterator iterEnd_=std::lower_bound(detIds_.begin(),
	                                                        detIds_.end(),
								*max_element(tmp.begin(), tmp.end()));
								
      for(std::vector<uint32_t>::iterator detIter_ = iterEnd_;
                                         detIter_!= iterBegin_-1;detIter_--){
	  detIds_.erase(detIter_);
      } 

    } // loop SubDetectorsToBeExcluded_
  }

  
  // -----
  // *** fill only one Module per layer ***

  if(fPSet_.getParameter<std::string>("ModulesToBeFilled") == "onlyOneModulePerLayer"){
   
    std::vector<uint32_t> tmp;
    std::vector<uint32_t> layerDetIds;

    SiStripSubStructure substructure_;
    
    for(unsigned int i=1; i<5 ; i++){
      tmp.clear();
      substructure_.getTIBDetectors(detIds_, tmp, i,0,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
    }
    for(unsigned int i=1; i<7 ; i++){
      tmp.clear();
      substructure_.getTOBDetectors(detIds_, tmp, i,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
    }
    for(unsigned int i=1; i<4 ; i++){
      tmp.clear();
      substructure_.getTIDDetectors(detIds_, tmp, 1,i,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
      substructure_.getTIDDetectors(detIds_, tmp, 2,i,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
    }  
    for(unsigned int i=1; i<10 ; i++){
      tmp.clear();
      substructure_.getTECDetectors(detIds_, tmp, 1,i,0,0,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
      substructure_.getTECDetectors(detIds_, tmp, 2,i,0,0,0,0);
      if(tmp.size() !=0) { layerDetIds.push_back(*(tmp.begin()));}
    }
    
    detIds_.clear();
    detIds_=layerDetIds;

  }
  // -----

   
} //selectModules
// -----



// -----
void SiStripBaseCondObjDQM::getModMEs(ModMEs& CondObj_ME, const uint32_t& detId_){
  
  std::map< uint32_t, ModMEs >::const_iterator ModMEsMap_iter = ModMEsMap_.find(detId_);

  if (ModMEsMap_iter != ModMEsMap_.end()){ 
  
    CondObj_ME=ModMEsMap_iter->second;
  
    if( ( CondObj_fillId_ =="ProfileAndCumul" || CondObj_fillId_ =="onlyProfile") && CondObj_ME.ProfileDistr ) { 
      CondObj_ME.ProfileDistr ->Reset();    
    }
       
    if( (CondObj_fillId_ =="ProfileAndCumul" || CondObj_fillId_ =="onlyCumul" ) &&  CondObj_ME.CumulDistr ){
      CondObj_ME.CumulDistr ->Reset();
    }
    else {
      edm::LogWarning("SiStripBaseCondObjDQM") 
	<< "[SiStripBaseCondObjDQM::getModMEs] PLEASE CHECK : CondObj_fillId option mispelled";
    }
    return; 

  }
  
  // --> profile defined for all CondData
  if ( (CondObj_fillId_ =="ProfileAndCumul" || CondObj_fillId_ =="onlyProfile") ) {
    bookProfileMEs(CondObj_ME,detId_);
  }  
  
  // --> cumul currently only defined for noise and apvgain
  if( (CondObj_fillId_ =="ProfileAndCumul" || CondObj_fillId_ =="onlyCumul" )
      &&(CondObj_name_ == "noise" || CondObj_name_ == "apvgain")          ) bookCumulMEs(CondObj_ME,detId_);
  
 
  ModMEsMap_.insert( std::make_pair(detId_,CondObj_ME) );
  
}
// ---- 


// -----
void SiStripBaseCondObjDQM::getSummaryMEs(ModMEs& CondObj_ME, const uint32_t& detId_){

  std::map<uint32_t, ModMEs>::const_iterator SummaryMEsMap_iter;

  if(CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_ ){
    SummaryMEsMap_iter = SummaryMEsMap_.find(getStringNameAndId(detId_).second);
  }
  else {
    SummaryMEsMap_iter = SummaryMEsMap_.find(getLayerNameAndId(detId_).second);
  }
   
  if (SummaryMEsMap_iter != SummaryMEsMap_.end()){ return;}

  // --> currently only profile summary defined for all condition objects except quality
  if(  CondObj_name_ == "pedestal"      || 
      CondObj_name_ == "noise"         || 
      CondObj_name_ == "apvgain"       || 
      CondObj_name_ == "lorentzangle" ) {
    if(FillSummaryProfileAtLayerLevel)	
      if (CondObj_ME.SummaryOfProfileDistr) { bookSummaryProfileMEs(CondObj_ME,detId_);}  
    
  }
    
  // --> currently only genuine cumul LA
  if(  CondObj_name_ == "lorentzangle" ||  CondObj_name_ == "noise"  ) {
    if(FillCumulativeSummaryAtLayerLevel)
      if (CondObj_ME.SummaryOfCumulDistr) { bookSummaryCumulMEs(CondObj_ME,detId_); } 
  } 
                          
  // --> currently only summary as a function of detId for noise, pedestal and apvgain 
  if(      CondObj_name_ == "noise"          ||
          CondObj_name_ == "apvgain"        || 
	  CondObj_name_ == "pedestal"       || 
	  CondObj_name_ == "quality"           ) {
    if(FillSummaryAtLayerLevel)          
      if (CondObj_ME.SummaryDistr) { bookSummaryMEs(CondObj_ME,detId_); } 
    
  } 
                          
  if(CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_) {
    SummaryMEsMap_.insert( std::make_pair(getStringNameAndId(detId_).second,CondObj_ME) );
  }
  else {
    SummaryMEsMap_.insert( std::make_pair(getLayerNameAndId(detId_).second,CondObj_ME) );
  }

}
// ---- 



// -----
void SiStripBaseCondObjDQM::bookProfileMEs(SiStripBaseCondObjDQM::ModMEs& CondObj_ME, const uint32_t& detId_){
     
  int   hProfile_NchX  = 0;
  double hProfile_LowX  = 0;
  double hProfile_HighX = 0;
  
  std::string hProfile_description;
  hProfile_description   = hPSet_.getParameter<std::string>("Profile_description");
      
  std::string hProfile_xTitle, hProfile_yTitle;
  hProfile_xTitle          = hPSet_.getParameter<std::string>("Profile_xTitle");
  hProfile_yTitle          = hPSet_.getParameter<std::string>("Profile_yTitle");
            
  if( CondObj_name_!= "apvgain" ){
    
    int nStrip      = reader->getNumberOfApvsAndStripLength(detId_).first*128;
	
    hProfile_NchX           = nStrip;
    hProfile_LowX           = 0.5;
    hProfile_HighX          = nStrip+0.5;
  }
  else {
	
    int nApv      = reader->getNumberOfApvsAndStripLength(detId_).first;
	
    hProfile_NchX           = nApv;
    hProfile_LowX           = 0.5;
    hProfile_HighX          = nApv+0.5;
  }
      
  folder_organizer.setDetectorFolder(detId_); 
      
  std::string hProfile_Name; 
  hProfile_Name = hidmanager.createHistoId(hProfile_description, "det", detId_); ;
      
  std::string hProfile;
  hProfile = hProfile_Name ;
    
  CondObj_ME.ProfileDistr = dqmStore_->book1D(hProfile_Name, hProfile, hProfile_NchX, hProfile_LowX, hProfile_HighX);
  CondObj_ME.ProfileDistr->setAxisTitle(hProfile_xTitle,1);
  CondObj_ME.ProfileDistr->setAxisTitle(hProfile_yTitle,2);
  dqmStore_->tag(CondObj_ME.ProfileDistr, detId_);
  
}
// -----


      
// -----
void SiStripBaseCondObjDQM::bookCumulMEs(SiStripBaseCondObjDQM::ModMEs& CondObj_ME, const uint32_t& detId_){

  int    hCumul_NchX    = 0;
  double hCumul_LowX    = 0;
  double hCumul_HighX   = 0;
     
  std::string hCumul_description;
  hCumul_description   = hPSet_.getParameter<std::string>("Cumul_description");
      
  std::string hCumul_xTitle, hCumul_yTitle;
  hCumul_xTitle        = hPSet_.getParameter<std::string>("Cumul_xTitle");
  hCumul_yTitle        = hPSet_.getParameter<std::string>("Cumul_yTitle");
      
  hCumul_NchX          = hPSet_.getParameter<int>("Cumul_NchX");
  hCumul_LowX          = hPSet_.getParameter<double>("Cumul_LowX");
  hCumul_HighX         = hPSet_.getParameter<double>("Cumul_HighX");
      
  folder_organizer.setDetectorFolder(detId_); 
      
  std::string hCumul_name; 
  hCumul_name   = hidmanager.createHistoId(hCumul_description  , "det", detId_); ;
      
  std::string hCumul_title;
  hCumul_title   = hCumul_name ;
      
  CondObj_ME.CumulDistr = dqmStore_->book1D(hCumul_name, 
                                            hCumul_title, 
					    hCumul_NchX, 
					    hCumul_LowX, 
					    hCumul_HighX);
  CondObj_ME.CumulDistr->setAxisTitle(hCumul_xTitle,1);
  CondObj_ME.CumulDistr->setAxisTitle(hCumul_yTitle,2);
  dqmStore_->tag(CondObj_ME.CumulDistr, detId_);
      
} 
// ---- 



// -----
void SiStripBaseCondObjDQM::bookSummaryProfileMEs(SiStripBaseCondObjDQM::ModMEs& CondObj_ME, const uint32_t& detId_){
  
  std::vector<uint32_t> sameLayerDetIds_;

  int   hSummaryOfProfile_NchX    = 0;
  double hSummaryOfProfile_LowX    = 0;
  double hSummaryOfProfile_HighX   = 0;
     
  std::string hSummaryOfProfile_description;
  hSummaryOfProfile_description  = hPSet_.getParameter<std::string>("SummaryOfProfile_description");
      
  std::string hSummaryOfProfile_xTitle, hSummaryOfProfile_yTitle;
  hSummaryOfProfile_xTitle        = hPSet_.getParameter<std::string>("SummaryOfProfile_xTitle");
  hSummaryOfProfile_yTitle        = hPSet_.getParameter<std::string>("SummaryOfProfile_yTitle");
  
  int hSummaryOfProfile_NchY;
  double hSummaryOfProfile_LowY, hSummaryOfProfile_HighY;
  hSummaryOfProfile_NchY          = hPSet_.getParameter<int>("SummaryOfProfile_NchY");
  hSummaryOfProfile_LowY          = hPSet_.getParameter<double>("SummaryOfProfile_LowY");
  hSummaryOfProfile_HighY         = hPSet_.getParameter<double>("SummaryOfProfile_HighY");
  
  int nStrip, nApv, layerId_;    
  
  if(CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_) { layerId_= getStringNameAndId(detId_).second;}
  else                                                          { layerId_= getLayerNameAndId(detId_).second;}


     
  if( CondObj_name_ == "pedestal" || CondObj_name_ == "noise" ){ // plot in strip number
    
    if( (layerId_ > 610 && layerId_ < 620) || // TID & TEC have 768 strips at maximum
        (layerId_ > 620 && layerId_ < 630) ||
        (layerId_ > 410 && layerId_ < 414) ||
        (layerId_ > 420 && layerId_ < 424) ){ nStrip =768;} 
    else { nStrip      = reader->getNumberOfApvsAndStripLength(detId_).first*128;}
    
    hSummaryOfProfile_NchX           = nStrip;
    hSummaryOfProfile_LowX           = 0.5;
    hSummaryOfProfile_HighX          = nStrip+0.5;
  
  }  
  else if( (CondObj_name_ == "lorentzangle" && SummaryOnLayerLevel_On_) || CondObj_name_ == "quality"){ // plot in detId-number

    // -----
    // get detIds belonging to same layer to fill X-axis with detId-number
  					   
    uint32_t subDetId_ =  ((detId_>>25)&0x7);
    SiStripSubStructure substructure_;
  
    sameLayerDetIds_.clear();
  
    if(subDetId_==3){  //  TIB
      substructure_.getTIBDetectors(activeDetIds, sameLayerDetIds_,TIBDetId(detId_).layerNumber(),0,0,TIBDetId(detId_).stringNumber());  
    }
    else if(subDetId_==4){  // TID
      substructure_.getTIDDetectors(activeDetIds, sameLayerDetIds_,0,0,0,0);
    }
    else if(subDetId_==5){  // TOB
      substructure_.getTOBDetectors(activeDetIds, sameLayerDetIds_, TOBDetId(detId_).layerNumber(),0,0);
    }
    else if(subDetId_==6){  // TEC
      substructure_.getTECDetectors(activeDetIds, sameLayerDetIds_, 0,0,0,0,0,0);
    }

    hSummaryOfProfile_NchX           = sameLayerDetIds_.size(); 
    hSummaryOfProfile_LowX           = 0.5;
    hSummaryOfProfile_HighX          = sameLayerDetIds_.size()+0.5;
 
  } 
  else if( CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_){ // plot in detId-number

    // -----
    // get detIds belonging to same string to fill X-axis with detId-number
  					   
    uint32_t subDetId_ =  ((detId_>>25)&0x7);
    SiStripSubStructure substructure_;
    
    sameLayerDetIds_.clear(); 
    
    if(subDetId_==3){  //  TIB    
      if(TIBDetId(detId_).isInternalString()){
      	substructure_.getTIBDetectors(activeDetIds, sameLayerDetIds_, TIBDetId(detId_).layerNumber(),0,1,TIBDetId(detId_).stringNumber()); }
      else if(TIBDetId(detId_).isExternalString()){
	substructure_.getTIBDetectors(activeDetIds, sameLayerDetIds_, TIBDetId(detId_).layerNumber(),0,2,TIBDetId(detId_).stringNumber()); }
    }
    else if(subDetId_==4){  // TID
      substructure_.getTIDDetectors(activeDetIds, sameLayerDetIds_, 0,0,0,0);
    }
    else if(subDetId_==5){  // TOB
      substructure_.getTOBDetectors(activeDetIds, sameLayerDetIds_, TOBDetId(detId_).layerNumber(),0,TOBDetId(detId_).rodNumber());
    }
    else if(subDetId_==6){  // TEC
      substructure_.getTECDetectors(activeDetIds, sameLayerDetIds_, 0,0,0,0,0,0);
    }


    hSummaryOfProfile_NchX           = sameLayerDetIds_.size(); 
    hSummaryOfProfile_LowX           = 0.5;
    hSummaryOfProfile_HighX          = sameLayerDetIds_.size()+0.5;
 
  } 
  else if( CondObj_name_ == "apvgain"){
 
    if( (layerId_ > 610 && layerId_ < 620) || // TID & TEC have 6 apvs at maximum
        (layerId_ > 620 && layerId_ < 630) ||
        (layerId_ > 410 && layerId_ < 414) ||
        (layerId_ > 420 && layerId_ < 424) ){ nApv =6;} 
    else { nApv     = reader->getNumberOfApvsAndStripLength(detId_).first;}
    
    hSummaryOfProfile_NchX           = nApv;
    hSummaryOfProfile_LowX           = 0.5;
    hSummaryOfProfile_HighX          = nApv+0.5;
 
  }
  else {
    edm::LogWarning("SiStripBaseCondObjDQM") 
       << "[SiStripBaseCondObjDQM::bookSummaryProfileMEs] PLEASE CHECK : x-axis label in your cfg"
       << std::endl; 
  }
  
  uint32_t layer_=0;
      
  layer_ = folder_organizer.GetSubDetAndLayer(detId_).second;
      
  folder_organizer.setLayerFolder(detId_,layer_); 
      
  std::string hSummaryOfProfile_name; 
  
  // ---
  int subdetectorId_ = ((detId_>>25)&0x7);
  
 
  if( subdetectorId_<3 ||subdetectorId_>6 ){ 
    edm::LogError("SiStripBaseCondObjDQM")
       << "[SiStripBaseCondObjDQM::bookSummaryProfileMEs] WRONG INPUT : no such subdetector type : "
       << subdetectorId_ << " no folder set!" 
       << std::endl;
    return;
  }
  // ---
  
  if(CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_) { 
    hSummaryOfProfile_name = hidmanager.createHistoLayer(hSummaryOfProfile_description, "layer" , getStringNameAndId(detId_).first,"") ;
  }
  else {
    hSummaryOfProfile_name = hidmanager.createHistoLayer(hSummaryOfProfile_description, "layer" , getLayerNameAndId(detId_).first,"") ;
  }
      
  std::string hSummaryOfProfile_title;
  hSummaryOfProfile_title   = hSummaryOfProfile_name ;
           
  CondObj_ME.SummaryOfProfileDistr = dqmStore_->bookProfile(hSummaryOfProfile_name, 
                                                            hSummaryOfProfile_title, 
						            hSummaryOfProfile_NchX, 
						            hSummaryOfProfile_LowX, 
						            hSummaryOfProfile_HighX, 
						            hSummaryOfProfile_NchY, 
						            hSummaryOfProfile_LowY, 
						            hSummaryOfProfile_HighY);
  CondObj_ME.SummaryOfProfileDistr->setAxisTitle(hSummaryOfProfile_xTitle,1);
  CondObj_ME.SummaryOfProfileDistr->setAxisTitle(hSummaryOfProfile_yTitle,2);
  CondObj_ME.SummaryOfProfileDistr->setAxisRange(hSummaryOfProfile_LowY, hSummaryOfProfile_HighY,2);
 
  // -----
  // in order to get the right detId-number labelled in right bin of x-axis
  
  if( CondObj_name_ == "quality" ){
    
    unsigned int iBin=0;
    
    for(unsigned int i=0;i< sameLayerDetIds_.size(); i++){
    
      iBin++;
      char sameLayerDetIds_Name[1024];
      sprintf(sameLayerDetIds_Name,"%u",sameLayerDetIds_[i]);
      CondObj_ME.SummaryOfProfileDistr->setBinLabel(iBin, sameLayerDetIds_Name);
    
    }
  } 
  if( CondObj_name_ == "lorentzangle"){

    // Put the detIds for the -z side as following the geometrical order:
      reverse(sameLayerDetIds_.begin(), sameLayerDetIds_.begin()+sameLayerDetIds_.size()/2);

      unsigned int iBin=0;
       for(unsigned int i=0;i< sameLayerDetIds_.size(); i++){ 
	 iBin++;
	 if (!SummaryOnStringLevel_On_){
	   // remove the label for detIds:
	     CondObj_ME.SummaryOfProfileDistr->setBinLabel(iBin, "");
	 }  
     
	 if (SummaryOnStringLevel_On_){
      // Label with module position instead of detIds:
	   char sameLayerDetIds_Name[1024];
	   if(subdetectorId_==3){//re-abelling for TIB
	     if(TIBDetId(sameLayerDetIds_[i]).isZPlusSide()){
	       sprintf(sameLayerDetIds_Name,"%i",TIBDetId(sameLayerDetIds_[i]).module());}
	     else if(TIBDetId(sameLayerDetIds_[i]).isZMinusSide()){
	       sprintf(sameLayerDetIds_Name,"%i",-TIBDetId(sameLayerDetIds_[i]).module());}
	     CondObj_ME.SummaryOfProfileDistr->setBinLabel(iBin, sameLayerDetIds_Name);
	   }
	   else if(subdetectorId_==5){//re-abelling for TOB
	     if(TOBDetId(sameLayerDetIds_[i]).isZPlusSide())      { sprintf(sameLayerDetIds_Name,"%i",TOBDetId(sameLayerDetIds_[i]).module());}
	     else if(TOBDetId(sameLayerDetIds_[i]).isZMinusSide()) { sprintf(sameLayerDetIds_Name,"%i",-TOBDetId(sameLayerDetIds_[i]).module());}
	     CondObj_ME.SummaryOfProfileDistr->setBinLabel(iBin, sameLayerDetIds_Name);
	   }
	 }
       } 
 
       


  // -----
      
  dqmStore_->tag(CondObj_ME.SummaryOfProfileDistr, layer_);
      
  } // if "lorentzangle"

}
// ---- 



// -----
void SiStripBaseCondObjDQM::bookSummaryCumulMEs(SiStripBaseCondObjDQM::ModMEs& CondObj_ME, const uint32_t& detId_){
    
  int   hSummaryOfCumul_NchX    = 0;
  double hSummaryOfCumul_LowX    = 0;
  double hSummaryOfCumul_HighX   = 0;
	
  std::string hSummaryOfCumul_description;
  hSummaryOfCumul_description  = hPSet_.getParameter<std::string>("SummaryOfCumul_description");
	
  std::string hSummaryOfCumul_xTitle, hSummaryOfCumul_yTitle;
  hSummaryOfCumul_xTitle        = hPSet_.getParameter<std::string>("SummaryOfCumul_xTitle");
  hSummaryOfCumul_yTitle        = hPSet_.getParameter<std::string>("SummaryOfCumul_yTitle");
	
  hSummaryOfCumul_NchX          = hPSet_.getParameter<int>("SummaryOfCumul_NchX");
  hSummaryOfCumul_LowX          = hPSet_.getParameter<double>("SummaryOfCumul_LowX");
  hSummaryOfCumul_HighX         = hPSet_.getParameter<double>("SummaryOfCumul_HighX");
	
  uint32_t layer_=0;
	
  layer_ = folder_organizer.GetSubDetAndLayer(detId_).second;
	
  folder_organizer.setLayerFolder(detId_,layer_); 
	
  std::string hSummaryOfCumul_name; 
  
  // ---
  int subdetectorId_ = ((detId_>>25)&0x7);
  
  if( subdetectorId_<3 || subdetectorId_>6 ){ 
    edm::LogError("SiStripBaseCondObjDQM")
       << "[SiStripBaseCondObjDQM::bookSummaryCumulMEs] WRONG INPUT : no such subdetector type : "
       << subdetectorId_ << " no folder set!" 
       << std::endl;
    return;
  }
  // ---
  
  // LA Histos are plotted for each string:
  if(CondObj_name_ == "lorentzangle" && SummaryOnStringLevel_On_) { 
    hSummaryOfCumul_name = hidmanager.createHistoLayer(hSummaryOfCumul_description, "layer" , getStringNameAndId(detId_).first, "") ;
  }
  else {  
    hSummaryOfCumul_name = hidmanager.createHistoLayer(hSummaryOfCumul_description, "layer" , getLayerNameAndId(detId_).first, "") ;
  }

	
  std::string hSummaryOfCumul_title;
  hSummaryOfCumul_title   = hSummaryOfCumul_name ;
	
  CondObj_ME.SummaryOfCumulDistr = dqmStore_->book1D(hSummaryOfCumul_name, 
                                                     hSummaryOfCumul_title, 
						     hSummaryOfCumul_NchX, 
						     hSummaryOfCumul_LowX, 
						     hSummaryOfCumul_HighX);

  CondObj_ME.SummaryOfCumulDistr->setAxisTitle(hSummaryOfCumul_xTitle,1);
  CondObj_ME.SummaryOfCumulDistr->setAxisTitle(hSummaryOfCumul_yTitle,2);
	
  dqmStore_->tag(CondObj_ME.SummaryOfCumulDistr, layer_);
	
}
// -----


// -----
void SiStripBaseCondObjDQM::bookSummaryMEs(SiStripBaseCondObjDQM::ModMEs& CondObj_ME, const uint32_t& detId_){
  
  std::vector<uint32_t> sameLayerDetIds_;

  int   hSummary_NchX    = 0;
  double hSummary_LowX    = 0;
  double hSummary_HighX   = 0;
     
  std::string hSummary_description;
  hSummary_description  = hPSet_.getParameter<std::string>("Summary_description");
      
  std::string hSummary_xTitle, hSummary_yTitle;
  hSummary_xTitle        = hPSet_.getParameter<std::string>("Summary_xTitle");
  hSummary_yTitle        = hPSet_.getParameter<std::string>("Summary_yTitle");
  
  int hSummary_NchY;
  double hSummary_LowY, hSummary_HighY;
  hSummary_NchY          = hPSet_.getParameter<int>("Summary_NchY");
  hSummary_LowY          = hPSet_.getParameter<double>("Summary_LowY");
  hSummary_HighY         = hPSet_.getParameter<double>("Summary_HighY");
       

  // -----
  // get detIds belonging to same layer to fill X-axis with detId-number
  					   
  uint32_t subDetId_ =  ((detId_>>25)&0x7);
  SiStripSubStructure substructure_;
  
  sameLayerDetIds_.clear();
   
  if(subDetId_==3){  //  TIB
    substructure_.getTIBDetectors(activeDetIds, sameLayerDetIds_, TIBDetId(detId_).layerNumber(),0,0,0);  
  }
  else if(subDetId_==4){  // TID
    substructure_.getTIDDetectors(activeDetIds, sameLayerDetIds_, TIDDetId(detId_).side(),TIDDetId(detId_).diskNumber(),0,0);
  }
  else if(subDetId_==5){  // TOB
    substructure_.getTOBDetectors(activeDetIds, sameLayerDetIds_, TOBDetId(detId_).layerNumber(),0,0);
  }
  else if(subDetId_==6){  // TEC
    substructure_.getTECDetectors(activeDetIds, sameLayerDetIds_, TECDetId(detId_).side(), TECDetId(detId_).wheelNumber(),0,0,0,0);
  }

  hSummary_NchX           = sameLayerDetIds_.size(); 
  hSummary_LowX           = 0.5;
  hSummary_HighX          = sameLayerDetIds_.size()+0.5;
  
  uint32_t layer_=0;
      
  layer_ = folder_organizer.GetSubDetAndLayer(detId_).second;
      
  folder_organizer.setLayerFolder(detId_,layer_); 
      
  std::string hSummary_name; 
  
  // ---
  int subdetectorId_ = ((detId_>>25)&0x7);
  
 
  if( subdetectorId_<3 ||subdetectorId_>6 ){ 
    edm::LogError("SiStripBaseCondObjDQM")
       << "[SiStripBaseCondObjDQM::bookSummaryMEs] WRONG INPUT : no such subdetector type : "
       << subdetectorId_ << " no folder set!" 
       << std::endl;
    return;
  }
  // ---
  
  hSummary_name = hidmanager.createHistoLayer(hSummary_description, 
                                                       "layer" , 
						        getLayerNameAndId(detId_).first,
							"") ;
      
  std::string hSummary_title;
  hSummary_title   = hSummary_name ;
           
  CondObj_ME.SummaryDistr = dqmStore_->bookProfile(hSummary_name, 
                                                   hSummary_title, 
						   hSummary_NchX, 
						   hSummary_LowX, 
						   hSummary_HighX, 
						   hSummary_NchY, 
						   hSummary_LowY, 
						   hSummary_HighY);
  CondObj_ME.SummaryDistr->setAxisTitle(hSummary_xTitle,1);
  CondObj_ME.SummaryDistr->setAxisTitle(hSummary_yTitle,2);
  CondObj_ME.SummaryDistr->setAxisRange(hSummary_LowY, hSummary_HighY,2);
 
  // -----
  // in order to get the right detId-number labelled in right bin of x-axis
  unsigned int iBin=0;
    
  for(unsigned int i=0;i< sameLayerDetIds_.size(); i++){
    
    iBin++;
    char sameLayerDetIds_Name[1024];
    sprintf(sameLayerDetIds_Name,"%u",sameLayerDetIds_[i]);
    CondObj_ME.SummaryDistr->setBinLabel(iBin, " ");
    
  }
  // -----
      
  dqmStore_->tag(CondObj_ME.SummaryDistr, layer_);
      
} 




// -----
std::pair<std::string,uint32_t> SiStripBaseCondObjDQM::getLayerNameAndId(const uint32_t& detId_){

  int subdetectorId_ = ((detId_>>25)&0x7);
  int layerId_=0;
  std::string layerName_;
  
  char tempLayerName_;
  char tempLayerNumber_[20];
  
  if(      subdetectorId_==3 ){ //TIB
    
    for( unsigned int i=1; i < 5 ;i++){
      
      sprintf(tempLayerNumber_,"%u",i);
      
      if(TIBDetId(detId_).layer()==i){ 
        sprintf( &tempLayerName_, "%s%s","TIB__layer__", tempLayerNumber_); 
	layerId_ = 300 + i;
      }
      layerName_ = &tempLayerName_;
      
    }
    
  }
  else if( subdetectorId_==4 ){ //TIDD
    
    if(TIDDetId(detId_).side()==1){ // TIDD side 1
      
      for( unsigned int i=1; i < 4 ;i++){
	
	sprintf(tempLayerNumber_,"%u",i);
	
	if(TIDDetId(detId_).wheel()==i){ 
	  sprintf( &tempLayerName_, "%s%s","TID__side__1__wheel__", tempLayerNumber_); 
	  layerId_ = 410 + i;
	}
	layerName_ = &tempLayerName_;
	
      }
      
      
    }
    else if(TIDDetId(detId_).side()==2){// TIDD side 2
      
      for( unsigned int i=1; i < 4 ;i++){
	
	sprintf(tempLayerNumber_,"%u",i);
	
	if(TIDDetId(detId_).wheel()==i){ 
	  sprintf( &tempLayerName_, "%s%s","TID__side__2__wheel__", tempLayerNumber_); 
	  layerId_ = 420 + i;
	}
	layerName_ = &tempLayerName_;
	
      }
      
    }
  }
  else if( subdetectorId_==5 ){ // TOB
    
    for( unsigned int i=1; i < 7 ;i++){
      
      sprintf(tempLayerNumber_,"%u",i);
      
      if(TOBDetId(detId_).layer()==i){ 
        sprintf( &tempLayerName_, "%s%s","TOB__layer__", tempLayerNumber_); 
	layerId_ = 500 + i;
      }
      layerName_ = &tempLayerName_;
      
    }
    
    
  }
  else if( subdetectorId_==6 ){ // TEC
    
    
    if(TECDetId(detId_).side()==1){ // TEC side 1
      
      for( unsigned int i=1; i < 10 ;i++){
	
	sprintf(tempLayerNumber_,"%u",i);
	
	if(TECDetId(detId_).wheel()==i){ 
	  sprintf( &tempLayerName_, "%s%s","TEC__side__1__wheel__", tempLayerNumber_); 
	  layerId_ = 610 + i;
	}
	layerName_ = &tempLayerName_;
	
      }
    }
    else if(TECDetId(detId_).side()==2){ // TEC side 2
      
      for( unsigned int i=1; i < 10 ;i++){
	
	sprintf(tempLayerNumber_,"%u",i);
	
	if(TECDetId(detId_).wheel()==i){ 
	  sprintf( &tempLayerName_, "%s%s","TEC__side__2__wheel__", tempLayerNumber_); 
	  layerId_ = 620 + i;
	}
	layerName_ = &tempLayerName_;
	
      }
     }
  }
  
  return std::make_pair(layerName_,layerId_);
}


//---------------


std::pair<std::string,uint32_t> SiStripBaseCondObjDQM::getStringNameAndId(const uint32_t& detId_){

  int subdetectorId_ = ((detId_>>25)&0x7);
  int layerStringId_=0;
  std::string layerStringName_;
  
  char tempLayerStringName_;
  char tempLayerStringNumber_[100];
  
  if(      subdetectorId_==3 ){ //TIB
    if(TIBDetId(detId_).layer()==1 && TIBDetId(detId_).isInternalString()){ //1st layer int
      for( unsigned int i=1; i < 27 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L1_Int_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 30110+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else  if(TIBDetId(detId_).layer()==1 && TIBDetId(detId_).isExternalString()){ //1st layer ext
      for( unsigned int i=1; i < 31 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L1_Ext_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 301200+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==2 && TIBDetId(detId_).isInternalString()){ //2nd layer int
      for( unsigned int i=1; i < 35 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L2_Int_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 302100+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==2 && TIBDetId(detId_).isExternalString()){ //2nd layer ext
      for( unsigned int i=1; i < 39 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L2_Ext_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 302200+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==3 && TIBDetId(detId_).isInternalString()){ //3rd layer int
      for( unsigned int i=1; i < 45 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L3_Int_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 303100+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==3 && TIBDetId(detId_).isExternalString()){ //3rd layer ext
      for( unsigned int i=1; i < 47 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L3_Ext_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 303200+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==4 && TIBDetId(detId_).isInternalString()){ //4th layer int
      for( unsigned int i=1; i < 53 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L4_Int_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 304100+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TIBDetId(detId_).layer()==4 && TIBDetId(detId_).isExternalString()){ //4th layer ext
      for( unsigned int i=1; i < 57 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TIBDetId(detId_).stringNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TIB_L4_Ext_Str_", tempLayerStringNumber_);  
	  layerStringId_ = 304200+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
  } //TIB


  else if( subdetectorId_==5 ){ // TOB
    if(TOBDetId(detId_).layer()==1){ //1st layer
      for( unsigned int i=1; i < 43 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L1_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50100+i;
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TOBDetId(detId_).layer()==2){ //2nd layer
      for( unsigned int i=1; i < 49 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L2_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50200+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TOBDetId(detId_).layer()==3){ //3rd layer
      for( unsigned int i=1; i < 55 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L3_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50300+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TOBDetId(detId_).layer()==4){ //4th layer
      for( unsigned int i=1; i < 61 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L4_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50400+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TOBDetId(detId_).layer()==5){ //5th layer
      for( unsigned int i=1; i < 67 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L5_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50500+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
    else if(TOBDetId(detId_).layer()==6){ //6st layer
      for( unsigned int i=1; i < 75 ;i++){
	sprintf(tempLayerStringNumber_,"%u",i); 
	if(TOBDetId(detId_).rodNumber()==i){  
	  sprintf(&tempLayerStringName_, "%s%s","TOB_L6_Rod_", tempLayerStringNumber_);  
	  layerStringId_ = 50600+i; 
	  layerStringName_ = &tempLayerStringName_; 
	}
      }      
    }
  }//TOB

  // std::cout<<"TEST__layerStringName\t"<<layerStringName_<<"\t detId\t"<<detId_<<std::endl;
  return std::make_pair(layerStringName_,layerStringId_);
}
    
//========================
std::vector<uint32_t> SiStripBaseCondObjDQM::GetSameLayerDetId(std::vector<uint32_t> activeDetIds,uint32_t selDetId ){
 
  std::vector<uint32_t> sameLayerDetIds;
  sameLayerDetIds.clear();

  SiStripSubStructure substructure_;
  
  uint32_t subselDetId_ =  ((selDetId>>25)&0x7);

  if(subselDetId_==3){  //  TIB
    substructure_.getTIBDetectors(activeDetIds, sameLayerDetIds, TIBDetId(selDetId).layer(),0,0,0);  
  }
  else if(subselDetId_==4){  // TID
    substructure_.getTIDDetectors(activeDetIds, sameLayerDetIds, TIDDetId(selDetId).side(),TIDDetId(selDetId).wheel(),0,0);
  }
  else if(subselDetId_==5){  // TOB
    substructure_.getTOBDetectors(activeDetIds, sameLayerDetIds, TOBDetId(selDetId).layer(),0,0);
  }
  else if(subselDetId_==6){  // TEC
    substructure_.getTECDetectors(activeDetIds, sameLayerDetIds, TECDetId(selDetId).side(),TECDetId(selDetId).wheel(),0,0,0,0);
  }

  return sameLayerDetIds;
  
}


