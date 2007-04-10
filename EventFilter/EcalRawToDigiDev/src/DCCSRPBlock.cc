
#include "EventFilter/EcalRawToDigiDev/interface/DCCSRPBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataBlockPrototype.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEventBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEventBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/ECALUnpackerException.h"

#include "EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h"


DCCSRPBlock::DCCSRPBlock(
  DCCDataUnpacker * u,EcalElectronicsMapper * m, DCCEventBlock * e, bool unpack
) : DCCDataBlockPrototype(u,m,e,unpack)
{
 
  // Todo : include data integrity collections
  blockLength_    = SRP_BLOCKLENGTH;
  // Set SR flags to zero
  for(uint i=0; i<SRP_NUMBFLAGS; i++){ srFlags_[i]=0; }

}


void DCCSRPBlock::unpack(uint64_t ** data, uint * dwToEnd, uint numbFlags ){    
  
  expNumbSrFlags_ = numbFlags;
  error_          = false;  
  datap_          = data;
  data_           = *data;
  dwToEnd_        = dwToEnd;
  
  // Check SRP Length
  if( (*dwToEnd_) < blockLength_ ){
    std::ostringstream output;
    output<<"EcalRawToDigi@SUB=DCCSRPBlock::unpack"
      <<"\n Unable to unpack SRP block for event "<<event_->l1A()<<" in dcc <<"<<mapper_->getActiveDCC()
      <<"\n Only "<<((*dwToEnd_)*8)<<" bytes are available while "<<(blockLength_*8)<<" are needed!"<<std::endl;
    //Note : add to error collection 
    throw ECALUnpackerException(output.str());
  }
  
  
  
  // Point to begin of block
  data_++;
  
  srpId_          = ( *data_ ) & SRP_ID_MASK; 
  bx_             = ( *data_>>SRP_BX_B     ) & SRP_BX_MASK;
  l1_             = ( *data_>>SRP_L1_B     ) & SRP_L1_MASK;
  nSRFlags_       = ( *data_>>SRP_NFLAGS_B ) & SRP_NFLAGS_MASK;
 
  checkSrpIdAndNumbSRFlags(); 
	 
  // Check synchronization
  if(sync_){
    uint dccL1 = ( event_->l1A() ) & SRP_BX_MASK;
    uint dccBx = ( event_->bx()  ) & SRP_L1_MASK;
    if( dccBx != bx_ || dccL1 != l1_ ){
      std::ostringstream output;
      output<<"EcalRawToDigi@SUB=DCCSRPBlock::unpack"
        <<"\nSynchronization error for SRP block in event "<<event_->l1A()<<" with bx "<<event_->bx()<<" in dcc <<"<<mapper_->getActiveDCC()
        <<"\n SRP local l1A is  "<<l1_<<" and local bx is "<<bx_;
       //Note : add to error collection ?		 
       throw ECALUnpackerException(output.str());
    }
  }  
  //display(cout); 
  addSRFlagToCollection();
  
  
  updateEventPointers();
        
}



void DCCSRPBlock::display(std::ostream& o){

  o<<"\n Unpacked Info for SRP Block"
  <<"\n DW1 ============================="
  <<"\n SRP Id "<<srpId_
  <<"\n Numb Flags "<<nSRFlags_ 
  <<"\n Bx "<<bx_
  <<"\n L1 "<<l1_;
 
  for(uint i=0; i<SRP_NUMBFLAGS; i++){ 
    o<<"\n SR flag "<<(i+1)<<" = "<<(srFlags_[i]); 
  } 
} 




