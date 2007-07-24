
#include "EventFilter/EcalRawToDigiDev/interface/DCCEBSRPBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataBlockPrototype.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEventBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEventBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h"



DCCEBSRPBlock::DCCEBSRPBlock(
  DCCDataUnpacker * u,EcalElectronicsMapper * m, DCCEventBlock * e, bool unpack
) : DCCSRPBlock(u,m,e,unpack)
{
  
  expNumbSrFlags_ = SRP_EB_NUMBFLAGS;
  
}


void DCCEBSRPBlock::updateCollectors(){
 // Set SR flag digis
  ebSrFlagsDigis_ = unpacker_->ebSrFlagsCollection(); 
}


void DCCEBSRPBlock::addSRFlagToCollection(){
  
  // Point to SR flags 
  data_++;
  uint16_t * my16Bitp_ = reinterpret_cast<uint16_t *> (data_);


  for( uint n=0; n<expNumbSrFlags_ ; n++ ){
    
    if(n>0&&n%4==0) my16Bitp_++;
   
    ushort  srFlag =  ( *my16Bitp_ >> ( (n-(n/4)*4) * 3 ) )  &  SRP_SRFLAG_MASK ;
    srFlags_[n] = srFlag;

    if(unpackInternalData_){  
     EBSrFlag * sr = (EBSrFlag*) mapper_->getSrFlagPointer(n+1);
     sr->setValue(srFlag);
     (*ebSrFlagsDigis_)->push_back(*sr);
    } 
  }
  

}



bool DCCEBSRPBlock::checkSrpIdAndNumbSRFlags(){

   //todo : check srp id based on sm...

  // Check number of SR flags
  if( nSRFlags_ != expNumbSrFlags_ ){
    edm::LogWarning("EcalRawToDigi@SUB=DCCSRPBlock::unpack")
      <<"\nUnable to unpack SRP block for event "<<event_->l1A()<<" in dcc <<"<<mapper_->getActiveDCC()
      <<"\nNumber of flags "<<nSRFlags_<<" is different from expected "<<expNumbSrFlags_;
    //Note : add to error collection ?		 
    return false;
  }
  return true;

} 

