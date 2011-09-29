#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripAPVRestorer.h"

#include <cmath>
#include <iostream>
#include <algorithm>


SiStripAPVRestorer::SiStripAPVRestorer(const edm::ParameterSet& conf):
  quality_cache_id(-1), noise_cache_id(-1), pedestal_cache_id(-1),
  ForceNoRestore_(conf.getParameter<bool>("ForceNoRestore")),
  SelfSelectRestoreAlgo_(conf.getParameter<bool>("SelfSelectRestoreAlgo")),
  InspectAlgo_(conf.getParameter<std::string>("APVInspectMode")),
  RestoreAlgo_(conf.getParameter<std::string>("APVRestoreMode")),
  useRealMeanCM_(conf.getParameter<bool>("useRealMeanCM")),
  fraction_(conf.getParameter<double>("Fraction")),
  deviation_(conf.getParameter<uint32_t>("Deviation")),
  restoreThreshold_(conf.getParameter<double>("restoreThreshold")),
  DeltaCMThreshold_(conf.getParameter<uint32_t>("DeltaCMThreshold")),
  nSigmaNoiseDerTh_(conf.getParameter<uint32_t>("nSigmaNoiseDerTh")),
  consecThreshold_(conf.getParameter<uint32_t>("consecThreshold")),
  hitStripThreshold_(conf.getParameter<uint32_t>("hitStripThreshold")),
  nSmooth_(conf.getParameter<uint32_t>("nSmooth")),
  minStripsToFit_(conf.getParameter<uint32_t>("minStripsToFit")),
  distortionThreshold_(conf.getParameter<uint32_t>("distortionThreshold")),
  CutToAvoidSignal_(conf.getParameter<double>("CutToAvoidSignal")),
  nSaturatedStrip_(conf.getParameter<uint32_t>("nSaturatedStrip")),
  ApplyBaselineCleaner_(conf.getParameter<bool>("ApplyBaselineCleaner")),
  MeanCM_(conf.getParameter<int32_t>("MeanCM"))  
  
{
  apvFlags_.clear();
  median_.clear();
  SmoothedMaps_.clear();
  BaselineMap_.erase(BaselineMap_.begin(), BaselineMap_.end());
}


void SiStripAPVRestorer::init(const edm::EventSetup& es){
  uint32_t n_cache_id = es.get<SiStripNoisesRcd>().cacheIdentifier();
  uint32_t q_cache_id = es.get<SiStripQualityRcd>().cacheIdentifier();
  uint32_t p_cache_id = es.get<SiStripPedestalsRcd>().cacheIdentifier();
  
  if(n_cache_id != noise_cache_id) {
    es.get<SiStripNoisesRcd>().get( noiseHandle );
    noise_cache_id = n_cache_id;
  } else {
    noise_cache_id = n_cache_id;
  }
  if(q_cache_id != quality_cache_id) {
    es.get<SiStripQualityRcd>().get( qualityHandle );
    quality_cache_id = q_cache_id;
  }else {
    quality_cache_id = q_cache_id;
  }
  
  if(p_cache_id != pedestal_cache_id) {
		es.get<SiStripPedestalsRcd>().get( pedestalHandle );
		pedestal_cache_id = p_cache_id;
  }else {
    pedestal_cache_id = p_cache_id;
  }
  
}

 
int16_t SiStripAPVRestorer::InspectAndRestore( const uint32_t& detId, const uint16_t& firstAPV, std::vector<int16_t>& rawDigisPedSubtracted,  std::vector<int16_t>& processedRawDigi, const std::vector< std::pair<short,float> >& vmedians ){
  int16_t nAPVFlagged = this->inspect(detId, firstAPV, rawDigisPedSubtracted, vmedians);
  this->restore(firstAPV, processedRawDigi);
  return nAPVFlagged;
}


int16_t SiStripAPVRestorer::inspect( const uint32_t& detId, const uint16_t& firstAPV, std::vector<int16_t>& digis, const std::vector< std::pair<short,float> >& vmedians) {
  
  detId_ = detId;
  
  apvFlags_.clear();
  apvFlags_.insert(apvFlags_.begin(), 6, "");
  median_.clear();
  median_.insert(median_.begin(), 6, -999);
  badAPVs_.clear();
  badAPVs_.insert(badAPVs_.begin(), 6, false);
  SmoothedMaps_.erase(SmoothedMaps_.begin(),  SmoothedMaps_.end());
  BaselineMap_.erase(BaselineMap_.begin(), BaselineMap_.end()); 
    
  for(size_t i=0; i< vmedians.size(); ++i){
         short APV =  vmedians[i].first;
         median_[APV]= vmedians[i].second;
         badAPVs_[APV] = qualityHandle->IsApvBad(detId_, APV);
  }
	
  if(InspectAlgo_=="BaselineFollower") return this->BaselineFollowerInspect(firstAPV, digis); 
  if(InspectAlgo_=="AbnormalBaseline") return this->AbnormalBaselineInspect(firstAPV, digis);
  if(InspectAlgo_=="Null") return this->NullInspect(firstAPV, digis);
  if(InspectAlgo_=="BaselineAndSaturation") return this->BaselineAndSaturationInspect(firstAPV, digis);
  throw cms::Exception("Unregistered Inspect Algorithm") << "SiStripAPVRestorer possibilities: (Null), (AbnormalBaseline),(BaselineFollower)";
  
}


void SiStripAPVRestorer::restore(const uint16_t& firstAPV, std::vector<int16_t>& digis ) {
	
  if(ForceNoRestore_) return;
  
  for( uint16_t APV=firstAPV; APV< digis.size()/128 + firstAPV; ++APV){
    std::string	algoToUse = *( apvFlags_.begin() + APV );
    
    if ( algoToUse != ""){
      if(!SelfSelectRestoreAlgo_) algoToUse = RestoreAlgo_;
     
   
  
      if(algoToUse=="Flat"){
	this->FlatRestore(APV, firstAPV, digis);
      }else if(algoToUse=="BaselineFollower"){
	this->BaselineFollowerRestore(APV, firstAPV, median_[APV], digis);
      }else{
	throw cms::Exception("Unregistered Restore Algorithm") << "SiStripAPVRestorer possibilities: (Flat), (BaselineFollower)";
      }
      
      
    }
  }
  
}


//Inspect method implementation ==========================================
//========================================================================
template<typename T>
inline
int16_t SiStripAPVRestorer::BaselineFollowerInspect(const uint16_t& firstAPV, std::vector<T>& digis){
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
  std::vector<T> singleAPVdigi;
  singleAPVdigi.clear();
 
  int16_t nAPVflagged = 0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  for(uint16_t APV=firstAPV ; APV< digis.size()/128 + firstAPV; ++APV){

    DigiMap smoothedmap;
    
    if(!badAPVs_[APV]){
      float MeanAPVCM = MeanCM_;
      if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
    
      singleAPVdigi.clear(); 
      for(int16_t strip = (APV-firstAPV)*128; strip < (APV-firstAPV+1)*128; ++strip){
        singleAPVdigi.push_back(digis[strip]); 
      }
    
    
      float DeltaCM = median_[APV] - MeanAPVCM; 
      
      //std::cout << "Delta CM: " << DeltaCM << " CM: " << median_[APV] << " detId " << (uint32_t) detId_ << std::endl; 	
      if(DeltaCM < 0 && std::abs(DeltaCM) > DeltaCMThreshold_){
      
        bool isFlat = FlatRegionsFinder(singleAPVdigi,smoothedmap,APV);
        if(!isFlat){
	      apvFlags_[APV]= "BaselineFollower";    //specify any algo to make the restore
	      nAPVflagged++;
        }
      }	
      
    } 
    SmoothedMaps_.insert(SmoothedMaps_.end(), std::pair<uint16_t, DigiMap>(APV, smoothedmap));
   }
  
  return nAPVflagged;
}







//Restore method implementation =====================================
//===================================================================
inline
void SiStripAPVRestorer::BaselineFollowerRestore(const uint16_t& APVn, const uint16_t& firstAPV, const float& median, std::vector<int16_t>& digis){
  //typename std::vector<T>::iterator firstStrip(digis.begin() + APVn*128), lastStrip(firstStrip + 128), actualStrip;
  
  
  std::vector<int16_t> baseline;
  baseline.clear();
  baseline.insert(baseline.begin(),128, 0);
  	
	 
 
    
  //============================= Find Flat Regions & Interpolating the baseline & subtracting the baseline  =================	
  
  if(SmoothedMaps_.size()){
    std::map<uint16_t, DigiMap >::iterator itSmootedMap = SmoothedMaps_.find(APVn);
    //this->BaselineFollower(SmoothedMaps_[APVn], baseline, median);	
    this->BaselineFollower(itSmootedMap->second, baseline, median);
  } else {
    //median=0;
    DigiMap  smoothedpoints;
    std::vector<int16_t> singleAPVdigi;
    singleAPVdigi.clear(); 
    for(int16_t strip = (APVn-firstAPV)*128; strip < (APVn-firstAPV+1)*128; ++strip) singleAPVdigi.push_back(digis[strip]); 
    this->FlatRegionsFinder(singleAPVdigi,smoothedpoints, APVn);
    this->BaselineFollower(smoothedpoints, baseline, median);		
    
  }	
  
  //============================= subtracting the baseline =============================================
  
  for(int16_t itStrip= 0 ; itStrip< 128; ++itStrip){
    digis[(APVn-firstAPV)*128+itStrip] -= baseline[itStrip] - median;
  }
  
		
  //============================= storing baseline to the map =============================================	
  BaselineMap_.insert(BaselineMap_.end(),  std::pair< uint16_t, std::vector < int16_t> >(APVn, baseline));
  
}


inline
void SiStripAPVRestorer::FlatRestore(const uint16_t& APVn, const uint16_t& firstAPV, std::vector<int16_t>& digis ){
 
  std::vector<int16_t> baseline;
  baseline.clear();
  baseline.insert(baseline.begin(),128, 150);
  baseline[0]=0; baseline[127]=0;
  BaselineMap_.insert(BaselineMap_.end(),  std::pair< uint16_t, std::vector < int16_t> >(APVn, baseline));  
  
  // typename std::vector<T>::iterator strip(digis.begin() + APVn*128), lastStrip(strip + 128);
  
  for(int16_t itStrip= 0 ; itStrip< 128; ++itStrip){
    digis[(APVn-firstAPV)*128+itStrip] = baseline[itStrip];
  }
 
  
}



//Baseline calculation implementation ==============================================================================================
//==================================================================================================================================

bool inline SiStripAPVRestorer::FlatRegionsFinder(const std::vector<int16_t>& adcs, DigiMap& smoothedpoints, const uint16_t& APVn){
  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId_);
  
  DigiMap consecpoints;
  DigiMapIter itConsecpoints, itSmoothedpoints;
  consecpoints.erase(consecpoints.begin(), consecpoints.end());
  smoothedpoints.erase(smoothedpoints.begin(), smoothedpoints.end());
  
   
  //============================= Height above local minimum ===============================                    
  std::vector<float> adcsLocalMinSubtracted;
  adcsLocalMinSubtracted.clear();
  adcsLocalMinSubtracted.insert(adcsLocalMinSubtracted.begin(), 128,0);
  for(uint32_t istrip=0; istrip<128; ++istrip) {
    float localmin = 999.9;		
    for(uint16_t jstrip=std::max(0,(int)(istrip-nSmooth_/2)); jstrip<std::min(128,(int)(istrip+nSmooth_/2)); ++jstrip) {
      float nextvalue = adcs[jstrip];
      if(nextvalue < localmin) localmin=nextvalue;			
    }
    adcsLocalMinSubtracted[istrip] = adcs[istrip] - localmin;
  }
  
  
  //============================= Find regions with stable slopes ========================
  std::vector<uint16_t> nConsStrip;
  nConsStrip.clear();
  
  //Creating maps with all the neighborhood strip and putting in a nCosntStip vector how many we have
  uint16_t consecStrips=0;
  for(uint32_t istrip=0; istrip<128; ++istrip) {    
    int16_t adc = adcs[istrip]; 
 
   //if( adcsLocalMinSubtracted[istrip] < nSigmaNoiseDerTh_ * (float)noiseHandle->getNoise(istrip+APVn*128,detNoiseRange) && (adc - median) < hitStripThreshold_){
   if( adcsLocalMinSubtracted[istrip] < nSigmaNoiseDerTh_ * (float)noiseHandle->getNoiseFast(istrip+APVn*128,detNoiseRange)){
      consecpoints.insert(consecpoints.end(), std::pair<uint16_t, int16_t >(istrip, adc));
      ++consecStrips;
    }else if (consecStrips >0){
      nConsStrip.push_back(consecStrips);
      consecStrips = 0;
    }    
  }     		

  //to cope with the last flat region of the APV
  if(consecStrips >0) nConsStrip.push_back(consecStrips);

  //removing from the map the fist and last points in wide flat regions and erasing from the map too small regions
  itConsecpoints = consecpoints.begin();
  float MinSmoothValue=20000., MaxSmoothValue=0.;
  for(std::vector<uint16_t>::iterator itnConsStrip = nConsStrip.begin(); itnConsStrip < nConsStrip.end(); ++itnConsStrip){
    
    consecStrips = *itnConsStrip;
    if(consecStrips >=consecThreshold_){
      ++itConsecpoints;  //skipping first point
      uint16_t nFirstStrip = itConsecpoints->first;
      uint16_t nLastStrip;
      float smoothValue = 0.0;
      float stripCount =1;
      for(uint16_t n =0; n < consecStrips-2; ++n){
		smoothValue += itConsecpoints->second;
		if(stripCount == consecThreshold_){
		  smoothValue /= (float)stripCount;
	  	  nLastStrip = nFirstStrip + stripCount -1;				                    
	  	  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
		  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
		  if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
		  if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
		  nFirstStrip = nLastStrip+1;
		  smoothValue=0;
		  stripCount=0;
		}
		++stripCount;
		++itConsecpoints;
     }
     ++itConsecpoints;  //and putting the pointer to the new seies of point 
      
     if(stripCount>1) {
     //if(smoothValue>0){
		--stripCount;
		smoothValue /= (float)(stripCount);
		nLastStrip = nFirstStrip + stripCount -1;
		smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
		smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
		if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
		if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
     }
   } else{
      for(int n =0; n< consecStrips ; ++n) ++itConsecpoints;
   }
  }
  
  	
  if( (MaxSmoothValue-MinSmoothValue) > distortionThreshold_){
 	if(ApplyBaselineCleaner_) this->BaselineCleaner(adcs, smoothedpoints, APVn);
	return false;
  }
  return true;
}


void inline SiStripAPVRestorer::BaselineCleaner(const std::vector<int16_t>& adcs, DigiMap& smoothedpoints, const uint16_t& APVn){
	SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId_);

	// only run the cleaner if there are enough points to start with
	if(smoothedpoints.size() < 4) return;
	
   	DigiMapIter itSmoothedpoints, itSmoothedpointsNext, itSmoothedpointsBegin, itSmoothedpointsEnd;
	
	bool printout = false;
   // if(436311512 == detId_) printout = true;
     if(printout){
        itSmoothedpointsBegin = smoothedpoints.begin();
     	std::cout << "start cleaning ============================= " << detId_ << " size " << smoothedpoints.size() << std::endl; 
     	for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != smoothedpoints.end(); ++itSmoothedpoints){  
     		std::cout << "strip " << itSmoothedpoints->first << " adc " << itSmoothedpoints->second << std::endl;
     	}
     	std::cout << "===============================================" << std::endl; 
     }
	
	itSmoothedpoints=smoothedpoints.begin();
	while ( itSmoothedpoints != --(smoothedpoints.end()) ) { //while we are not at the last point
	    if(smoothedpoints.size() <2) break;
		// get info about current and next points
		itSmoothedpointsNext = itSmoothedpoints;
		++itSmoothedpointsNext;
		float strip1 = itSmoothedpoints->first;
		float strip2 = itSmoothedpointsNext->first;
		float adc1 = itSmoothedpoints->second;
		float adc2 = itSmoothedpointsNext->second;
	  	float m = (adc2 -adc1)/(strip2 -strip1);
		if(printout) std::cout << smoothedpoints.size() << " " << strip1 << " " << adc1 << " " << strip2 << " " << adc2 << std::endl;		
		if (m>2) { // in case of large positive slope, remove next point and try again from same current point
			smoothedpoints.erase(itSmoothedpointsNext);
		} else if (m<-2) { // in case of large negative slope, remove current point and either...
			// move to next point if we have reached the beginning (post-increment to avoid invalidating pointer during erase) or...
			if(itSmoothedpoints==smoothedpoints.begin()) smoothedpoints.erase(itSmoothedpoints++); 
			// try again from the previous point if we have not reached the beginning
			else smoothedpoints.erase(itSmoothedpoints--); 
		} else { // in case of a flat enough slope, continue on to the next point
			itSmoothedpoints++;
		}
		
	}
	
	
	 if(printout){
	   	std::cout << "ending cleaning ============================= " << detId_  << " size " << smoothedpoints.size() << std::endl;
		itSmoothedpointsBegin = smoothedpoints.begin();
		for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != smoothedpoints.end(); ++itSmoothedpoints){  
		  std::cout << "strip " << itSmoothedpoints->first << " adc " << itSmoothedpoints->second << std::endl;
		}
		std::cout << "===============================================" <<  std::endl; 
		if(printout) printout = false;
	}		
	
	

	//inserting extra point is case of local minimum
	//--------------------------------------------------------------------------------------------------
	// these should be reset now for the point-insertion that follows
	
	if(smoothedpoints.size() >= 2){
    	itSmoothedpointsBegin = smoothedpoints.begin();
    	itSmoothedpointsEnd = --(smoothedpoints.end());
		for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != itSmoothedpointsEnd; ++itSmoothedpoints){  
    		itSmoothedpointsNext = itSmoothedpoints;
			++itSmoothedpointsNext;
      		float strip1 = itSmoothedpoints->first;
      		float strip2 = itSmoothedpointsNext->first;
      		float adc1 = itSmoothedpoints->second;
      		float adc2 = itSmoothedpointsNext->second;
	  		float m = (adc2 -adc1)/(strip2 -strip1);
    
        
        	if((strip2 - strip1) >3 && abs(adc1 -adc2) >4){
				float itStrip = 1;
        		float strip = itStrip + strip1;
 				while(strip < strip2){
				
					float adc = adcs[strip];
                	if( adc < (adc1 + m * itStrip - 2 * (float)noiseHandle->getNoiseFast(strip+APVn*128,detNoiseRange))){
						//std::cout << "applying correction strip: " << strip + APVn*128 << " adc " << adc << " detId: " << detId_ << std::endl;
						smoothedpoints.insert(itSmoothedpointsNext, std::pair<uint16_t, int16_t >(strip,adc));
						++itSmoothedpoints;
						++itSmoothedpointsNext;
						itSmoothedpointsEnd = --(smoothedpoints.end());
					} 
					++itStrip;
					++strip;
				}
			

	    	}
		}
	}
	
	
    itSmoothedpointsBegin = smoothedpoints.begin();
    itSmoothedpointsEnd = --(smoothedpoints.end());
    uint16_t firstStripFlat = itSmoothedpointsBegin->first;
    uint16_t lastStripFlat = itSmoothedpointsEnd->first;
    int16_t firstStripFlatADC= itSmoothedpointsBegin->second;
    int16_t lastStripFlatADC= itSmoothedpointsEnd->second;
    	
    itSmoothedpoints = itSmoothedpointsBegin;
    if(firstStripFlat >3){
		float strip = 0;
       	while(strip < firstStripFlat){
			float adc = adcs[strip];
            if( adc < ( firstStripFlatADC - 2 * (float)noiseHandle->getNoiseFast(strip+APVn*128,detNoiseRange))){
					smoothedpoints.insert(itSmoothedpoints, std::pair<uint16_t, int16_t >(strip,adc));
					++itSmoothedpoints;
			} 
			++strip;
		}
	}
	
	itSmoothedpoints = itSmoothedpointsEnd;
	if(lastStripFlat <125){
		float strip = lastStripFlat+1;
       	while(strip < 128){
			float adc = adcs[strip];
            if( adc < ( lastStripFlatADC - 2 * (float)noiseHandle->getNoiseFast(strip+APVn*128,detNoiseRange))){
            	smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(strip,adc));
			} 
			++strip;
		}
	}
	
	
    
}

void inline SiStripAPVRestorer::BaselineFollower(DigiMap& smoothedpoints, std::vector<int16_t>& baseline, const float& median){
  
  baseline.clear();
  DigiMapIter itSmoothedpoints;
  
  
  //if not enough points
  if(smoothedpoints.size() < minStripsToFit_){
     baseline.insert(baseline.begin(),128, median);
  } else {
     baseline.insert(baseline.begin(),128, 0);  
    
    DigiMapIter itSmoothedpointsBegin, itSmoothedpointsEnd;
    itSmoothedpointsBegin = smoothedpoints.begin();
    itSmoothedpointsEnd = --(smoothedpoints.end());
    
				
    uint16_t firstStripFlat = itSmoothedpointsBegin->first;
    uint16_t lastStripFlat = itSmoothedpointsEnd->first;
    int16_t firstStripFlatADC= itSmoothedpointsBegin->second;
    int16_t lastStripFlatADC= itSmoothedpointsEnd->second;
   
    //adding here the costant line at the extremities 
    baseline.erase(baseline.begin(), baseline.begin()+firstStripFlat);
    baseline.insert(baseline.begin(), firstStripFlat, firstStripFlatADC);
    
    baseline.erase(baseline.begin()+lastStripFlat, baseline.end());
    baseline.insert(baseline.end(), 128 - lastStripFlat, lastStripFlatADC);
    
    
    //IMPORTANT: the itSmoothedpointsEnd should be at least smaller than smoothedpoints.end() -1
    for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != itSmoothedpointsEnd; ++itSmoothedpoints){  
      DigiMapIter itSmoothedpointsNext = itSmoothedpoints;
      ++itSmoothedpointsNext;
      float strip1 = itSmoothedpoints->first;
      float strip2 = itSmoothedpointsNext->first;
      float adc1 = itSmoothedpoints->second;
      float adc2 = itSmoothedpointsNext->second;
     
      baseline[strip1] = adc1;
      baseline[strip2] = adc2;
      float m = (adc2 -adc1)/(strip2 -strip1);
      uint16_t itStrip = strip1 +1;
      float stripadc = adc1 + m; 
      while(itStrip < strip2){
		baseline[itStrip] = stripadc;
		++itStrip;
		stripadc+=m;
      }
      
    }
    
  }
}






//Other methods implementation ==============================================
//==========================================================================
/*
void SiStripAPVRestorer::fixAPVsCM(edm::DetSet<SiStripProcessedRawDigi>& cmdigis) {
  
  // cmdigis should be the same size as apvFlags_
  // otherwise something pathological has happened and we do nothing
  //if ( cmdigis.size() != apvFlags_.size() ) return;
  
  edm::DetSet<SiStripProcessedRawDigi>::iterator cm_iter = cmdigis.begin();
  std::vector<std::string>::const_iterator apvf_iter = apvFlags_.begin();
  
  // No way to change the adc value of a SiStripProcessedRawDigi
  // so we just extract the values, clear the DetSet, and
  // replace with the proper values.
  
  std::vector<float> cmvalues;
  for( ; cm_iter != cmdigis.end(); ++cm_iter  ) cmvalues.push_back( (*cm_iter).adc() );
  cmdigis.clear();
  
  std::vector<float>::const_iterator cmv_iter = cmvalues.begin();
  while( apvf_iter != apvFlags_.end() ){
    if( *apvf_iter != "") {
      cmdigis.push_back( SiStripProcessedRawDigi( -999.) );
    }
    else
      cmdigis.push_back( SiStripProcessedRawDigi( *cmv_iter ) );
    apvf_iter++;
    cmv_iter++;
  }
}
*/

void SiStripAPVRestorer::LoadMeanCMMap(const edm::Event& iEvent){
  if(useRealMeanCM_){  
	edm::Handle< edm::DetSetVector<SiStripRawDigi> > input;
    iEvent.getByLabel("siStripDigis","VirginRaw", input);
   this->CreateCMMapRealPed(*input);
  } else {
    edm::Handle< edm::DetSetVector<SiStripProcessedRawDigi> > inputCM;
    iEvent.getByLabel("MEANAPVCM",inputCM);
    this->CreateCMMapCMstored(*inputCM);
  }
}


void SiStripAPVRestorer::CreateCMMapRealPed(const edm::DetSetVector<SiStripRawDigi>& input){
  
  MeanCMmap_.erase(MeanCMmap_.begin(), MeanCMmap_.end());
  	
 //std::cout<< "===============================================" << std::endl;
 
 for ( edm::DetSetVector<SiStripRawDigi>::const_iterator 
	  rawDigis = input.begin(); rawDigis != input.end(); rawDigis++) {
         SiStripPedestals::Range detPedestalRange = pedestalHandle->getRange(rawDigis->id);
		 std::vector<float> MeanCMDetSet;
		 MeanCMDetSet.clear();
		
		for(uint16_t APV = 0; APV < rawDigis->size()/128; ++APV){
			uint16_t MinPed =0;
			for(uint16_t strip = APV*128; strip< (APV+1)*128; ++strip){
			  uint16_t ped =  (uint16_t)pedestalHandle->getPed(strip,detPedestalRange);
			  //std::cout << "Pedestal: " << ped << " strip: " << strip << " detId: " <<  rawDigis->id << std::endl;
			  if(ped < MinPed) MinPed = ped;
			}
			if(MinPed>128) MinPed=128;
			MeanCMDetSet.push_back(MinPed);
			//std::cout<< "Mean CM: "<< (uint32_t)rawDigis->id << ", " << MinPed << std::endl;
		}
		MeanCMmap_.insert(std::pair<uint32_t, std::vector<float> >(rawDigis->id,MeanCMDetSet));
		
	}

   //std::cout<< "***********************************************" << std::endl;
   //std::cout<< "***********************************************" << std::endl;
 
}

void SiStripAPVRestorer::CreateCMMapCMstored(const edm::DetSetVector<SiStripProcessedRawDigi>& Input){

  MeanCMmap_.erase(MeanCMmap_.begin(), MeanCMmap_.end());
  uint32_t detId;
  edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itInput;
  edm::DetSet<SiStripProcessedRawDigi>::const_iterator itCM;
  std::vector<float> MeanCMNValue;
  
  for(itInput = Input.begin(); itInput != Input.end(); ++itInput){
    detId = itInput->id;
    MeanCMNValue.clear();
    for(itCM = itInput->begin(); itCM != itInput->end(); ++itCM) MeanCMNValue.push_back(itCM->adc()); 			
    MeanCMmap_.insert(std::pair<uint32_t, std::vector<float> >(detId,MeanCMNValue));
  }
}

std::vector<bool>& SiStripAPVRestorer::GetAPVFlags(){
    //apvf.clear();
    apvFlagsBool_.clear();
    for(size_t i =0; i < apvFlags_.size(); ++i){
      //if(apvFlags_[i] != "") apvf.push_back(true);
      // else apvf.push_back(false); 
      if(apvFlags_[i] != "") apvFlagsBool_.push_back(true);
      else apvFlagsBool_.push_back(false);
    }
    return apvFlagsBool_;	
}


///Code still to be reviewed.
///Is should go backat the beginnng at the file.
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================
//==========================================================================


template<typename T>
inline
int16_t SiStripAPVRestorer::BaselineAndSaturationInspect(const uint16_t& firstAPV, std::vector<T>& digis){
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
   
  std::vector<T> singleAPVdigi;
  singleAPVdigi.clear();
  
  
  int16_t nAPVflagged = 0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    apvFlags_.push_back( "" );
    if(!badAPVs_[APV]){
     float MeanAPVCM = MeanCM_;
     if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
    
     singleAPVdigi.clear();
   
     uint16_t nSatStrip =0;
     for(int16_t strip = APV*128; strip < (APV+1)*128; ++strip){
       singleAPVdigi.push_back(digis[strip]);
       if(digis[strip] >=1023) ++nSatStrip;
     }
    
     float DeltaCM = median_[APV] -MeanAPVCM; 
    
    
     if(DeltaCM < 0 && std::abs(DeltaCM) > DeltaCMThreshold_&&nSatStrip>= nSaturatedStrip_){
       apvFlags_[APV] = RestoreAlgo_;    //specify any algo to make the restore
       nAPVflagged++;
     } 
    }	
  }
  
  return nAPVflagged;
}


template<typename T>
inline
int16_t SiStripAPVRestorer::AbnormalBaselineInspect( const uint16_t& firstAPV, std::vector<T>& digis){

  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
  typename std::vector<T>::iterator fs;
  
  int16_t nAPVflagged=0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  
  int devCount = 0, qualityCount = 0, minstrip = 0; 
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    apvFlags_.push_back( "" );
    if(!badAPVs_[APV]){
      float MeanAPVCM = MeanCM_;
      if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
      for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip){
        fs = digis.begin() + istrip;
        if ( !qualityHandle->IsStripBad(detQualityRange,istrip) ){
	       qualityCount++; 
	       if ( std::abs((int) *fs - MeanAPVCM) > (int)deviation_ ){ 
                devCount++;
	            minstrip = std::min((int) *fs, minstrip);
           }
         }
      }
    
      if( devCount > fraction_ * qualityCount ) {
        apvFlags_[APV] = RestoreAlgo_;      //specify any algo to make the restore
        nAPVflagged++;
      } 
    } 
  }
  
  return nAPVflagged;
  
}



template<typename T>
inline
int16_t SiStripAPVRestorer::NullInspect(const uint16_t& firstAPV, std::vector<T>& digis){

  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);

  typename std::vector<T>::iterator fs;

  int16_t nAPVflagged = 0;

  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
   apvFlags_.push_back( "" );
   if(!badAPVs_[APV]){ 
     int zeroCount = 0, qualityCount = 0; 
     for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip){
       fs = digis.begin() + istrip;
       if ( !qualityHandle->IsStripBad(detQualityRange,istrip) ){
        qualityCount++; 
        if ( (int) *fs < 1 ) zeroCount++;
       }
      }
    
      if( zeroCount > restoreThreshold_ * qualityCount ) {
        apvFlags_[APV] = RestoreAlgo_;     //specify any algo to make the restore
        nAPVflagged++;
      } 
   } 
   }
 
  return nAPVflagged;

}
