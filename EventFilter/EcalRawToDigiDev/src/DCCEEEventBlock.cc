
#include "EventFilter/EcalRawToDigiDev/interface/DCCEEEventBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h"
#include "EventFilter/EcalRawToDigiDev/interface/EcalDCCHeaderRuntypeDecoder.h"
#include "EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCFEBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCSCBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEETCCBlock.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCEESRPBlock.h"
#include <sys/time.h>

#include <iomanip>
#include <sstream>


DCCEEEventBlock::DCCEEEventBlock( DCCDataUnpacker * u,EcalElectronicsMapper * m, bool hU, bool srpU, bool tccU, bool feU, bool memU) : 
  DCCEventBlock(u,m,hU,srpU,tccU,feU,memU)
{
 
  //Builds a tower unpacker block
  towerBlock_ = new DCCSCBlock(u,m,this,feUnpacking_); 
  
  //Builds a srp unpacker block
  srpBlock_   = new DCCEESRPBlock(u,m,this,srpUnpacking_);
  
  //Builds a tcc unpacker block
  tccBlock_   = new DCCEETCCBlock(u,m,this,tccUnpacking_);
  
 
}

void DCCEEEventBlock::unpackTCCBlocks(){

  std::vector<short>::iterator it;
  for(it=tccChStatus_.begin();it!=tccChStatus_.end();it++){
    if( (*it) != CH_TIMEOUT &&  (*it) != CH_DISABLED){
      tccBlock_->unpack(&data_,&dwToEnd_);
    }
  }
}

