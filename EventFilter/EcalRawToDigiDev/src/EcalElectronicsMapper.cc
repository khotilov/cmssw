#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h>
#include <Geometry/EcalMapping/interface/EcalElectronicsMapping.h>
#include <DataFormats/EcalDigi/interface/EBSrFlag.h>
#include <DataFormats/EcalDigi/interface/EESrFlag.h>
#include <EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h>

EcalElectronicsMapper::EcalElectronicsMapper( uint numbXtalTSamples, uint numbTriggerTSamples)
: pathToMapFile_(""),
numbXtalTSamples_(numbXtalTSamples),
numbTriggerTSamples_(numbTriggerTSamples),
mappingBuilder_(0)

{
	
  
  // Reset Arrays
  for(uint sm=0; sm < NUMB_SM; sm++){
    for(uint fe=0; fe< NUMB_FE; fe++){
	  
      for(uint strip=0; strip<NUMB_STRIP;strip++){
        for(uint xtal=0; xtal<NUMB_XTAL;xtal++){
		    
   	  //Reset DFrames and xtalDetIds
	  xtalDetIds_[sm][fe][strip][xtal]=0;
	}
      }
      //Reset SC Det Ids
      scDetIds_[sm][fe]=0;
      scEleIds_[sm][fe]=0;
      srFlags_[sm][fe]=0;
    }
  }
  
  
  //Reset TT det Ids
  for( uint tccid=0; tccid < NUMB_TCC; tccid++){
    for(uint tpg =0; tpg<NUMB_FE;tpg++){
      ttDetIds_[tccid][tpg]=0;
      ttTPIds_[tccid][tpg]=0;
    }
  }

  

  //Fill map sm id to tcc ids
  std::vector<uint> * ids;
  ids = new std::vector<uint>;
  ids->push_back(1); ids->push_back(18);ids->push_back(19);ids->push_back(36);
  mapSmIdToTccIds_[1]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(2); ids->push_back(3);ids->push_back(20);ids->push_back(21);
  mapSmIdToTccIds_[2]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(4); ids->push_back(5);ids->push_back(22);ids->push_back(23);
  mapSmIdToTccIds_[3]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(6); ids->push_back(7);ids->push_back(24);ids->push_back(25);
  mapSmIdToTccIds_[4]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(8); ids->push_back(9);ids->push_back(26);ids->push_back(27);
  mapSmIdToTccIds_[5]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(10); ids->push_back(11);ids->push_back(28);ids->push_back(29);
  mapSmIdToTccIds_[6]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(12); ids->push_back(13);ids->push_back(30);ids->push_back(31);
  mapSmIdToTccIds_[7]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(14); ids->push_back(15);ids->push_back(32);ids->push_back(33);
  mapSmIdToTccIds_[8]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(16); ids->push_back(17);ids->push_back(34);ids->push_back(35);
  mapSmIdToTccIds_[9]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(73); ids->push_back(90);ids->push_back(91);ids->push_back(108);
  mapSmIdToTccIds_[46]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(74); ids->push_back(75);ids->push_back(92);ids->push_back(93);
  mapSmIdToTccIds_[47]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(76); ids->push_back(77);ids->push_back(94);ids->push_back(95);
  mapSmIdToTccIds_[48]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(78); ids->push_back(79);ids->push_back(96);ids->push_back(97);
  mapSmIdToTccIds_[49]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(80); ids->push_back(81);ids->push_back(98);ids->push_back(99);  
  mapSmIdToTccIds_[50]= ids;
		 
  ids = new std::vector<uint>;
  ids->push_back(82); ids->push_back(83);ids->push_back(100);ids->push_back(101);
  mapSmIdToTccIds_[51]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(84); ids->push_back(85);ids->push_back(102);ids->push_back(103);
  mapSmIdToTccIds_[52]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(86); ids->push_back(87);ids->push_back(104);ids->push_back(105);
  mapSmIdToTccIds_[53]= ids;
		
  ids = new std::vector<uint>;
  ids->push_back(88); ids->push_back(89);ids->push_back(106);ids->push_back(107);
  mapSmIdToTccIds_[54]= ids;
	
  
  //Compute data block sizes
  unfilteredFEBlockLength_= computeUnfilteredFEBlockLength();  
  ebTccBlockLength_       = computeEBTCCBlockLength();
  eeTccBlockLength_       = computeEETCCBlockLength();
    

}


EcalElectronicsMapper::~EcalElectronicsMapper(){


  if( mappingBuilder_ ){ delete mappingBuilder_; }
 

  //DETETE ARRAYS
  for(uint sm=0; sm < NUMB_SM; sm++){
    for(uint fe=0; fe< NUMB_FE; fe++){
      for(uint strip=0; strip<NUMB_STRIP;strip++){
        for(uint xtal=0; xtal<NUMB_XTAL;xtal++){
	  if(xtalDetIds_[sm][fe][strip][xtal]){ 
            delete xtalDetIds_[sm][fe][strip][xtal]; 
          }
        }
      }

      if(scDetIds_[sm][fe]){ 
        delete scDetIds_[sm][fe];
        delete scEleIds_[sm][fe];
        delete srFlags_[sm][fe];
      }

    }
 
  }
  
  for( uint tccid=0; tccid < NUMB_TCC; tccid++){
    for(uint tpg =0; tpg<NUMB_FE;tpg++){
      if(ttDetIds_[tccid][tpg]){ 
        delete ttDetIds_[tccid][tpg];
        delete ttTPIds_[tccid][tpg];
      }
    }
  }


  pathToMapFile_.clear();
  
  
  std::map<uint, std::vector<uint> *>::iterator it;
  for(it = mapSmIdToTccIds_.begin(); it != mapSmIdToTccIds_.end(); it++ ){ delete (*it).second; }
  
  mapSmIdToTccIds_.clear();
 
}

void EcalElectronicsMapper::setEcalElectronicsMapping(const EcalElectronicsMapping * m){
  mappingBuilder_= m;
  fillMaps();

}

bool EcalElectronicsMapper::setActiveDCC(uint dccId){
   
  bool ret(true);
	
  //Update active dcc and associated smId
  dccId_ = dccId;
   
  smId_  = getSMId(dccId_);
	
  if(!smId_) ret = false;
	
  return ret;
	
 } 
 
 
bool EcalElectronicsMapper::setDCCMapFilePath(std::string aPath_){

  
  //try to open a dccMapFile in the given path
  std::ifstream dccMapFile_(aPath_.c_str());

  //if not successful return false
  if(!dccMapFile_.is_open()) return false;

  //else close file and accept given path
  dccMapFile_.close();
  pathToMapFile_ = aPath_;

  return true;
}


// bool EcalElectronicsMapper::readDCCMapFile(){

//   //try to open a dccMapFile in the given path
//   std::ifstream dccMapFile_(pathToMapFile_.c_str());
  
//   //if not successful return false
//   if(!dccMapFile_.is_open()) return false;
  
//   char lineBuf_[100];
//   uint SMId_,DCCId_;
//   // loop while extraction from file is possible
//   dccMapFile_.getline(lineBuf_,10);       //read line from file
//   while (dccMapFile_.good()) {
//     sscanf(lineBuf_,"%u:%u",&SMId_,&DCCId_);
//     myDCCMap_[SMId_] = DCCId_;
//     dccMapFile_.getline(lineBuf_,10);       //read line from file
//   }
  
  
//   return true;
  
// }

// bool EcalElectronicsMapper::readDCCMapFile(std::string aPath_){
//   //test if path is good
//   edm::FileInPath eff(aPath_);
  
//   if(!setDCCMapFilePath(eff.fullPath())) return false;

//   //read DCC map file
//   readDCCMapFile();
//   return true;
// }


bool EcalElectronicsMapper::makeMapFromVectors( std::vector<int>& orderedFedUnpackList,
					       std::vector<int>& orderedDCCIdList )
{

  // in case as non standard set of DCCId:FedId pairs was provided
  if ( orderedFedUnpackList.size() == orderedDCCIdList.size() &&
       orderedFedUnpackList.size() > 0)
    {
      edm::LogInfo("EcalElectronicsMapper") << "DCCIdList/FedUnpackList lists given. Being loaded.";
      
      std::string correspondence("list of pairs DCCId:FedId :  ");
      char           onePair[50];
      for (int v=0;  v< ((int)orderedFedUnpackList.size());  v++)	{
	myDCCMap_[ orderedDCCIdList[v]  ] = orderedFedUnpackList[v] ;
	
	sprintf( onePair, "  %d:%d",  orderedDCCIdList[v], orderedFedUnpackList[v]);
	std::string                 tmp(onePair);
	correspondence += tmp;
      }
      edm::LogInfo("EcalElectronicsMapper") << correspondence;
      
    }
  else    
    {  // default set of DCCId:FedId for ECAL

      edm::LogInfo("EcalElectronicsMapper") << "No input DCCIdList/FedUnpackList lists given for ECAL unpacker"
					    << "(or given with different number of elements). "
					    << " Loading default association DCCIdList:FedUnpackList,"
					    << "i.e.  1:601 ... 53:653,  54:654.";
      
      for (uint v=1; v<=54; v++)	{
	myDCCMap_[ v ] = (v+600) ;   }
    }

  return true;
}


std::ostream &operator<< (std::ostream& o, const EcalElectronicsMapper &aMapper_) {
  //print class information
  o << "---------------------------------------------------------";

  if(aMapper_.pathToMapFile_.size() < 1){
    o << "No correct input for DCC map has been given yet...";
  }
  else{
    o << "DCC Map (Map file: " << aMapper_.pathToMapFile_ << " )" << "SM id\t\tDCCid ";

    //get DCC map and iterator
    std::map<uint ,uint > aMap;
    aMap=aMapper_.myDCCMap_;
    std::map<uint ,uint >::iterator iter;

    //print info contained in map
    for(iter = aMap.begin(); iter != aMap.end(); iter++)
      o << iter->first << "\t\t" << iter->second;
  }

  o << "---------------------------------------------------------";
  return o;
}

  

uint EcalElectronicsMapper::computeUnfilteredFEBlockLength(){

  return ((numbXtalTSamples_-2)/4+1)*25+1; 

}


uint EcalElectronicsMapper::computeEBTCCBlockLength(){

  uint nTT=68;
  uint tf;
	  
  //TCC block size: header (8 bytes) + 17 words with 4 trigger primitives (17*8bytes)
  if( (nTT*numbTriggerTSamples_)<4 || (nTT*numbTriggerTSamples_)%4 ) tf=1;  
  else tf=0;
    
  return 1 + ((nTT*numbTriggerTSamples_)/4) + tf ;

}

uint EcalElectronicsMapper::computeEETCCBlockLength(){
  //Todo : implement multiple tt samples for the endcap
  return 9;  

}


uint EcalElectronicsMapper::getDCCId(uint aSMId_) const{
  //get iterator for SM id
  std::map<uint ,uint>::const_iterator it = myDCCMap_.find(aSMId_);

  //check if SMid exists and return DCC id
  if(it!= myDCCMap_.end()) return it->second;
 
  //error return
<<<<<<< EcalElectronicsMapper.cc
  edm::LogWarning("EcalElectronicsMapper") << "DCC requested for SM id: " << aSMId_ << " not found";
=======
  if( ! DCCDataUnpacker::silentMode_ ){
    edm::LogError("EcalElectronicsMapper") << "DCC requested for SM id: " << aSMId_ << " not found";
  }
>>>>>>> 1.23
  return 0;
}


uint EcalElectronicsMapper::getSMId(uint aDCCId_) const {
  //get iterator map
  std::map<uint ,uint>::const_iterator it;

  //try to find SM id for given DCC id
  for(it = myDCCMap_.begin(); it != myDCCMap_.end(); it++)
    if(it->second == aDCCId_) 
      return it->first;

  //error return
<<<<<<< EcalElectronicsMapper.cc
  edm::LogWarning("EcalEcalElectronicsMapper") << "SM requested DCC id: " << aDCCId_ << " not found";
=======
  if( ! DCCDataUnpacker::silentMode_ ){
    edm::LogError("EcalEcalElectronicsMapper") << "SM requested DCC id: " << aDCCId_ << " not found";
  }
>>>>>>> 1.23
  return 0;
}




void EcalElectronicsMapper::fillMaps(){

 
  for( int smId=1 ; smId<= 54; smId++){

    
    // Fill EB arrays  
    if( smId > 9 && smId < 46 ){
	 
      for(int feChannel =1; feChannel<=68; feChannel++){
		   
        uint tccId = smId + TCCID_SMID_SHIFT_EB;
		  
        // Builds Ecal Trigger Tower Det Id 

        uint rawid = (mappingBuilder_->getTrigTowerDetId(tccId, feChannel)).rawId();
        EcalTrigTowerDetId * ttDetId = new EcalTrigTowerDetId(rawid);
        ttDetIds_[tccId-1][feChannel-1] = ttDetId;
	EcalElectronicsId * ttEleId = new EcalElectronicsId(smId, feChannel, 1, 1);
	ttEleIds_[tccId-1][feChannel-1] = ttEleId;
        EcalTriggerPrimitiveDigi * tp     = new EcalTriggerPrimitiveDigi(*ttDetId);
        tp->setSize(numbTriggerTSamples_);
        for(uint i=0;i<numbTriggerTSamples_;i++){
          tp->setSample( i, EcalTriggerPrimitiveSample(0) );
        }
        ttTPIds_[tccId-1][feChannel-1]  = tp;   
	
        // Buil SRP Flag
        srFlags_[smId-1][feChannel-1] = new EBSrFlag(*ttDetId,0);
 
        for(uint stripId =1; stripId<=5; stripId++){
		    
	  for(uint xtalId =1;xtalId<=5;xtalId++){
 		  
	      EcalElectronicsId eid(smId,feChannel,stripId,xtalId);
	      EBDetId * detId = new EBDetId( (mappingBuilder_->getDetId(eid)).rawId());
	      xtalDetIds_[smId-1][feChannel-1][stripId-1][xtalId-1] = detId;
			 
	   } // close loop over xtals
	}// close loop over strips
		  		
      }// close loop over fechannels
		
    }//close loop over sm ids in the EB
    // Fill EE arrays (Todo : waiting SC correction)
    
    else{
	 
	std::vector<uint> * pTCCIds = mapSmIdToTccIds_[smId];
	std::vector<uint>::iterator it;
		
	for(it= pTCCIds->begin(); it!= pTCCIds->end(); it++){
			
          uint tccId = *it;
          
          // creating arrays of pointers for trigger objects
	  for(uint towerInTCC =1; towerInTCC <= numChannelsInDcc_[smId-1]; towerInTCC++){
              
              // Builds Ecal Trigger Tower Det Id 
              EcalTrigTowerDetId ttDetId = mappingBuilder_->getTrigTowerDetId(tccId, towerInTCC);
              
              ttDetIds_[tccId-1][towerInTCC-1] = new EcalTrigTowerDetId(ttDetId.rawId());
              EcalTriggerPrimitiveDigi * tp   = new EcalTriggerPrimitiveDigi(ttDetId);
              tp->setSize(numbTriggerTSamples_);
              for(uint i=0;i<numbTriggerTSamples_;i++){
                  tp->setSample( i, EcalTriggerPrimitiveSample(0) );
              }
              
              ttTPIds_[tccId-1][towerInTCC-1]  = tp;
              
	  }
        }
        

        
        // creating arrays of pointers for digi objects
	for(uint feChannel = 1; feChannel <= numChannelsInDcc_[smId-1]; feChannel++){

            // to avoid gap in CCU_id's
            if((smId==SECTOR_EEM_CCU_JUMP   || smId== SECTOR_EEP_CCU_JUMP) &&
               (MIN_CCUID_JUMP <= feChannel && feChannel <=MAX_CCUID_JUMP )
               ) continue;
            
          EcalScDetId scDetId = mappingBuilder_->getEcalScDetId(smId,feChannel);
          scDetIds_[smId-1][feChannel-1] = new EcalScDetId(scDetId.rawId());
	  scEleIds_[smId-1][feChannel-1] = new EcalElectronicsId(smId,feChannel,1,1);
          srFlags_[smId-1][feChannel-1]  = new EESrFlag( EcalScDetId( scDetId.rawId() ) , 0 ); 

          std::vector<DetId> ecalDetIds = mappingBuilder_->dccTowerConstituents(smId,feChannel);
          std::vector<DetId>::iterator it;
				  
          //EEDetIds	
          for(it = ecalDetIds.begin(); it!= ecalDetIds.end(); it++){
		    
          EcalElectronicsId ids = mappingBuilder_->getElectronicsId((*it));
			 
          int stripId    = ids.stripId();
          int xtalId     = ids.xtalId();
		 
          EEDetId * detId = new EEDetId((*it).rawId());
          xtalDetIds_[smId-1][feChannel-1][stripId-1][xtalId-1] = detId;
          
	}// close loop over tower constituents
	  
	}// close loop over  FE Channels		
	
    }// closing loop over sm ids in EE
	
  }
  

}

// number of readout channels (TT in EB, SC in EE) in a DCC
const uint  EcalElectronicsMapper::numChannelsInDcc_[NUMB_SM] = {34,32,33,33,32,34,33,41,33,    // EE -
                                                                 68,68,68,68,68,68,68,68,68,68, // EB-
                                                                 68,68,68,68,68,68,68,68,
                                                                 68,68,68,68,68,68,68,68,68,68, // EB+
                                                                 68,68,68,68,68,68,68,68,
                                                                 34,32,33,33,32,34,33,41,33};   // EE+
