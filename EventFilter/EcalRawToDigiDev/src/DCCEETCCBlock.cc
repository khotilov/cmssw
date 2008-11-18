#include "EventFilter/EcalRawToDigiDev/interface/DCCEETCCBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h"


DCCEETCCBlock::DCCEETCCBlock ( DCCDataUnpacker  * u,EcalElectronicsMapper * m, DCCEventBlock * e, bool unpack) : 
DCCTCCBlock(u,m,e,unpack)
{
   blockLength_ = 0;
  
}

void DCCEETCCBlock::updateCollectors(){
  tps_ = unpacker_->ecalTpsCollection();
}




void DCCEETCCBlock::addTriggerPrimitivesToCollection(){
  
  //point to trigger data
  data_++;
  bool processTPG2(true) ;

  uint16_t * tccP_= reinterpret_cast< uint16_t * >(data_);


  // Unpack TPG1 input block
  if(ps_){ 
	for(uint i = 1; i<= NUMB_PSEUDOSTRIPS; i++, tccP_++){
        //fill pseudostrip container
	}
  }
 
		
  // Unpack TPG1 output block 
  for(uint i = 1; i<= NUMB_TTS_TPG1; i++, tccP_++){
	 
    if( i<= nTTs_){
	  pTP_ = mapper_->getTPPointer(tccId_,i);
	  if(pTP_) pTP_->setSample(0, *tccP_ );
	}else {
	  processTPG2 = false;
	  break;
	}
  }
  
	

  if(processTPG2){

    // Unpack TPG2 input block
    if(ps_){ 
      for(uint i = 1; i<= NUMB_PSEUDOSTRIPS; i++, tccP_++){
          //fill pseudostrip container
	  }
	}
	
	// Unpack TPG2 output block 
    for(uint i = 1; i<= NUMB_TTS_TPG2; i++, tccP_++){
	  uint tt = i+NUMB_TTS_TPG1;
	  
	  if( tt <= nTTs_){
	    pTP_ = mapper_->getTPPointer(tccId_,tt);
	    if(pTP_) pTP_->setSample(0, *tccP_ ); 
	  }
	  else break;
    }
    
  }
	

}


bool DCCEETCCBlock::checkTccIdAndNumbTTs(){
	
	
  bool tccFound(false);
  bool errorOnNumbOfTTs(false);
  int  activeDCC =  mapper_->getActiveSM();
  std::vector<uint> * m = mapper_->getTccs(activeDCC);
  std::vector<uint>::iterator it;
  for(it= m->begin();it!=m->end();it++){
    if((*it) == tccId_){ 
      tccFound=true;
     
	  /*
	    expNumbTTs_= 28; //separate from inner and outer tcc 
	
	    For inner TCCs you expect 28 TTs (4 in phi, 7 in eta).
	    For outer TCCs you expect 16 TTs (4 in phi, 4 in eta).
    
	    to implement : map tccid-> number of tts
	  */
	  if( nTTs_ != 28 || nTTs_ !=16 ){
	    if( ! DCCDataUnpacker::silentMode_ ){
          edm::LogWarning("EcalRawToDigiDevTCC") 
		   <<"\n Error on event "<<event_->l1A()<<" with bx "<<event_->bx()<<" in fed <<"<<mapper_->getActiveDCC()
           <<"\n TCC id "<<tccId_<<" has "<<nTTs_<<" Trigger Towers ( 28 or 16 are the expected values in EE)"
           <<"\n => Skipping to next fed block...";
          errorOnNumbOfTTs = true; 
          //todo : add to error collection   
        }
	  }
    }
  }
	
  if(!tccFound){

    if( ! DCCDataUnpacker::silentMode_ ){
      edm::LogWarning("EcalRawToDigiDevTCC") 
        <<"\n Error on event "<<event_->l1A()<<" with bx "<<event_->bx()<<" in fed <<"<<mapper_->getActiveDCC()
        <<"\n TCC id "<<tccId_<<" is not valid for this dcc "
        <<"\n => Skipping to next fed block...";
       //todo : add to error collection   
     }
  }
  

 return (tccFound || errorOnNumbOfTTs);

}


uint DCCEETCCBlock::getLength(){

  uint64_t * temp = data_;
  temp++;

  ps_       = ( *temp>>TCC_PS_B ) & B_MASK;   
   
  uint numbTps = (NUMB_TTS_TPG1 + NUMB_TTS_TPG2);
  if(ps_){ numbTps += NUMB_PSEUDOSTRIPS;}
  
  uint length = numbTps/4;
  if(numbTps%4) length++ ;
  
  
  return length;

}


