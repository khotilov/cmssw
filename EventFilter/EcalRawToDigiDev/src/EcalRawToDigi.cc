#include "EventFilter/EcalRawToDigiDev/interface/EcalRawToDigi.h"
#include "EventFilter/EcalRawToDigiDev/interface/EcalElectronicsMapper.h"
#include "EventFilter/EcalRawToDigiDev/interface/DCCDataUnpacker.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"

#include "DataFormats/EcalRawData/interface/EcalListOfFEDS.h"


EcalRawToDigiDev::EcalRawToDigiDev(edm::ParameterSet const& conf):

  
  //define the list of FED to be unpacked
  fedUnpackList_(conf.getUntrackedParameter<std::vector<int> >("FEDs", std::vector<int>())),

  //define the ordered FED list
  orderedFedUnpackList_(conf.getUntrackedParameter<std::vector<int> >("orderedFedList", std::vector<int>())),

  //define the ordered DCCId list
  orderedDCCIdList_(conf.getUntrackedParameter<std::vector<int> >("orderedDCCIdList", std::vector<int>())),

  //get number of Xtal Time Samples
  numbXtalTSamples_(conf.getUntrackedParameter<uint>("numbXtalTSamples",10)),

  //Get number of Trigger Time Samples
  numbTriggerTSamples_(conf.getUntrackedParameter<uint>("numbTriggerTSamples",1)),
  
  //See if header unpacking is enabled
  headerUnpacking_(conf.getUntrackedParameter<bool>("headerUnpacking",true)),
 
  //See if srp unpacking is enabled
  srpUnpacking_(conf.getUntrackedParameter<bool>("srpUnpacking",true)),
  
  //See if tcc unpacking is enabled
  tccUnpacking_(conf.getUntrackedParameter<bool>("tccUnpacking",true)),
  
  //See if fe unpacking is enabled
  feUnpacking_(conf.getUntrackedParameter<bool>("feUnpacking",true)),
  
  //See if fe unpacking is enabled
  memUnpacking_(conf.getUntrackedParameter<bool>("memUnpacking",true)), 

  //See if syncCheck is enabled
  syncCheck_(conf.getUntrackedParameter<bool>("syncCheck",true)), 
  
  put_(conf.getUntrackedParameter<bool>("eventPut",false)),
  
  dataLabel_(conf.getUntrackedParameter<std::string>("InputLabel","source")),

  REGIONAL_(conf.getUntrackedParameter<bool>("DoRegional",false)),

  fedsLabel_(conf.getUntrackedParameter<std::string>("FedLabel","listfeds")),

  myMap_(0),
  
  theUnpacker_(0)

{
  
  first_ = true;
  mmm_ = new EcalElectronicsMapping();

  
  if( numbXtalTSamples_ <6 || numbXtalTSamples_>64 || (numbXtalTSamples_-2)%4 ){
    std::ostringstream output;
    output      <<"\n Unsuported number of xtal time samples : "<<numbXtalTSamples_
		<<"\n Valid Number of xtal time samples are : 6,10,14,18,...,62"; 
    edm::LogError("EcalRawToDigiDev")<< output.str();
    // todo : throw an execption
  }
  
  if( numbTriggerTSamples_ !=1 && numbTriggerTSamples_ !=4 && numbTriggerTSamples_ !=8  ){
    std::ostringstream output;
    output      <<"\n Unsuported number of trigger time samples : "<<numbTriggerTSamples_
		<<"\n Valid number of trigger time samples are :  1, 4 or 8"; 
    edm::LogError("EcalRawToDigiDev")<< output.str();
    // todo : throw an execption
  }
  
  //NA : testing
  //nevts_=0;
  //RUNNING_TIME_=0;

  // if there are FEDs specified to unpack fill the vector of the fedUnpackList_
  // else fill with the entire ECAL fed range (600-670)
  if (fedUnpackList_.empty()) 
    for (int i=FEDNumbering::getEcalFEDIds().first; i<=FEDNumbering::getEcalFEDIds().second; i++)
      fedUnpackList_.push_back(i);

  //print the FEDs to unpack to the logger
  std::ostringstream loggerOutput_;
  if(fedUnpackList_.size()!=0){
    for (uint i=0; i<fedUnpackList_.size(); i++) 
      loggerOutput_ << fedUnpackList_[i] << " ";
    edm::LogInfo("EcalRawToDigiDev") << "EcalRawToDigi will unpack FEDs ( " << loggerOutput_.str() << ")";
  }
  
  edm::LogInfo("EcalRawToDigiDev")
    <<"\n ECAL RawToDigi configuration:"
    <<"\n Header  unpacking is "<<headerUnpacking_
    <<"\n SRP Bl. unpacking is "<<srpUnpacking_
    <<"\n TCC Bl. unpacking is "<<tccUnpacking_  
    <<"\n FE  Bl. unpacking is "<<feUnpacking_
    <<"\n MEM Bl. unpacking is "<<memUnpacking_<<"\n";
  
  // Producer products :
  produces<EBDigiCollection>("ebDigis"); 
  produces<EEDigiCollection>("eeDigis");
  produces<EBSrFlagCollection>();
  produces<EESrFlagCollection>();
  produces<EcalRawDataCollection>();
  produces<EcalPnDiodeDigiCollection>();
  produces<EcalTrigPrimDigiCollection>("EBTT"); //change name
  produces<EcalTrigPrimDigiCollection>("EETT"); //change name
  
  // Integrity for xtal data
  produces<EBDetIdCollection>("EcalIntegrityGainErrors");
  produces<EBDetIdCollection>("EcalIntegrityGainSwitchErrors");
  produces<EBDetIdCollection>("EcalIntegrityGainSwitchStayErrors");
  produces<EBDetIdCollection>("EcalIntegrityChIdErrors");

  // Integrity Errors
  produces<EcalTrigTowerDetIdCollection>("EcalIntegrityTTIdErrors");
  produces<EcalTrigTowerDetIdCollection>("EcalIntegrityBlockSizeErrors");
 
  // Mem channels' integrity
  produces<EcalElectronicsIdCollection>("EcalIntegrityMemTtIdErrors");
  produces<EcalElectronicsIdCollection>("EcalIntegrityMemBlockSizeErrors");
  produces<EcalElectronicsIdCollection>("EcalIntegrityMemChIdErrors");
  produces<EcalElectronicsIdCollection>("EcalIntegrityMemGainErrors");


 
  // Build a new Electronics mapper and parse default map file
  myMap_ = new EcalElectronicsMapper(numbXtalTSamples_,numbTriggerTSamples_);

  // in case of external  text file (deprecated by HLT environment) 
  //  bool readResult = myMap_->readDCCMapFile(conf.getUntrackedParameter<std::string>("DCCMapFile",""));

  // use two arrays from cfg to establish DCCId:FedId. If they are empy, than use hard coded correspondence 
  bool readResult = myMap_->makeMapFromVectors(orderedFedUnpackList_, orderedDCCIdList_);


  if(!readResult){
    edm::LogError("EcalRawToDigiDev")<<"\n unable to read file : "
      <<conf.getUntrackedParameter<std::string>("DCCMapFile","");
  }
  
  // Build a new ECAL DCC data unpacker
  theUnpacker_ = new DCCDataUnpacker(myMap_,headerUnpacking_,srpUnpacking_,tccUnpacking_,feUnpacking_,memUnpacking_,syncCheck_);
   
}

void EcalRawToDigiDev::produce(edm::Event& e, const edm::EventSetup& es) {

  //double TIME_START = clock();
  //nevts_++; //NUNO


  if (first_) {

   edm::ESHandle< EcalElectronicsMapping > ecalmapping;
   es.get< EcalMappingRcd >().get(ecalmapping);

   const EcalElectronicsMapping* TheMapping = ecalmapping.product();
   *mmm_ = *TheMapping;
   myMap_ -> setEcalElectronicsMapping(mmm_);   // because I can not pass a const to SetEcalElectronicsMapping

   first_ = false;

  }

  // Get list of FEDS :
  std::vector<int> FEDS_to_unpack;
  if (REGIONAL_) {
        edm::Handle<EcalListOfFEDS> listoffeds;
        e.getByLabel(fedsLabel_, listoffeds);
        FEDS_to_unpack = listoffeds -> GetList();
  }



  // Step A: Get Inputs    

  edm::Handle<FEDRawDataCollection> rawdata;  
  e.getByLabel(dataLabel_,rawdata);


  // Step B: encapsulate vectors in actual collections and set unpacker pointers

  // create the collection of Ecal Digis
  std::auto_ptr<EBDigiCollection> productDigisEB(new EBDigiCollection);
  productDigisEB->reserve(1700);
  theUnpacker_->setEBDigisCollection(&productDigisEB);
  
  // create the collection of Ecal Digis
  std::auto_ptr<EEDigiCollection> productDigisEE(new EEDigiCollection);
  theUnpacker_->setEEDigisCollection(&productDigisEE);
  
  // create the collection for headers
  std::auto_ptr<EcalRawDataCollection> productDccHeaders(new EcalRawDataCollection);
  theUnpacker_->setDccHeadersCollection(&productDccHeaders); 

  // create the collection for invalid gains
  std::auto_ptr< EBDetIdCollection> productInvalidGains(new EBDetIdCollection);
  theUnpacker_->setInvalidGainsCollection(&productInvalidGains); 

  // create the collection for invalid gain Switch
  std::auto_ptr< EBDetIdCollection> productInvalidGainsSwitch(new EBDetIdCollection);
  theUnpacker_->setInvalidGainsSwitchCollection(&productInvalidGainsSwitch);
   
  // create the collection for invalid gain switch stay
  std::auto_ptr< EBDetIdCollection> productInvalidGainsSwitchStay(new EBDetIdCollection);
  theUnpacker_->setInvalidGainsSwitchStayCollection(&productInvalidGainsSwitch);
   
  // create the collection for invalid chids
  std::auto_ptr< EBDetIdCollection> productInvalidChIds(new EBDetIdCollection);
  theUnpacker_->setInvalidChIdsCollection(&productInvalidChIds);
         
  // create the collection for EB srflags       
  std::auto_ptr<EBSrFlagCollection> productEBSrFlags(new EBSrFlagCollection);
  theUnpacker_->setEBSrFlagsCollection(&productEBSrFlags);
  
  // create the collection for EB srflags       
  std::auto_ptr<EESrFlagCollection> productEESrFlags(new EESrFlagCollection);
  theUnpacker_->setEESrFlagsCollection(&productEESrFlags);
  
  // create the collection for EB tpgs
  std::auto_ptr<EcalTrigPrimDigiCollection> productEBTps(new EcalTrigPrimDigiCollection);
  theUnpacker_->setEBTpsCollection(&productEBTps);
  
   // create the collection for EE tpgs
  std::auto_ptr<EcalTrigPrimDigiCollection> productEETps(new EcalTrigPrimDigiCollection);
  theUnpacker_->setEETpsCollection(&productEETps);

  // create the collection for invalid TTIds
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidTTIds(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidTTIdsCollection(&productInvalidTTIds);
  
  // create the collection for invalid BlockLengths
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidBlockLengths(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidBlockLengthsCollection(&productInvalidBlockLengths);

  // MEMs Collections
  // create the collection for the Pn Diode Digis
  std::auto_ptr<EcalPnDiodeDigiCollection> productPnDiodeDigis(new EcalPnDiodeDigiCollection);
  theUnpacker_->setPnDiodeDigisCollection(&productPnDiodeDigis);

  // create the collection for invalid Mem Tt id 
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidMemTtIds(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidMemTtIdsCollection(& productInvalidMemTtIds);
  
  // create the collection for invalid Mem Block Size 
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidMemBlockSizes(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidMemBlockSizesCollection(& productInvalidMemBlockSizes);
  
  // create the collection for invalid Mem Block Size 
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidMemChIds(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidMemChIdsCollection(& productInvalidMemChIds);
 
  // create the collection for invalid Mem Gain Errors 
  std::auto_ptr<EcalElectronicsIdCollection> productInvalidMemGains(new EcalElectronicsIdCollection);
  theUnpacker_->setInvalidMemGainsCollection(& productInvalidMemGains);
 
  //  double TIME_START = clock(); 
  

  // Step C: unpack all requested FEDs    
  for (std::vector<int>::const_iterator i=fedUnpackList_.begin(); i!=fedUnpackList_.end(); i++) {

    if (REGIONAL_) {
      std::vector<int>::const_iterator fed_it = find(FEDS_to_unpack.begin(), FEDS_to_unpack.end(), *i);
      if (fed_it == FEDS_to_unpack.end()) continue;
    }

  
    // get fed raw data and SM id
    const FEDRawData & fedData = rawdata->FEDData(*i);
    int length = fedData.size();

    edm::LogInfo("EcalRawToDigiDev") << "raw data lenght: " << length ;
    //if data size is not null interpret data
    if ( length >= EMPTYEVENTSIZE ){
      
      if(myMap_->setActiveDCC(*i)){

	int smId = myMap_->getActiveSM();
	edm::LogInfo("EcalRawToDigiDev") << "Getting FED = " << *i <<"(SM = "<<smId<<")"<<" data size is: " << length;

	uint64_t * pData = (uint64_t *)(fedData.data());
	theUnpacker_->unpack( pData, static_cast<uint>(length),smId,*i);
      }
    }
    
  }// loop on FEDs
  
  //if(nevts_>1){   //NUNO
  //  double TIME_END = clock(); //NUNO
  //  RUNNING_TIME_ += TIME_END-TIME_START; //NUNO
  // }
  
  
  // Add collections to the event 
  
  if(put_){
    
    if( headerUnpacking_){ 
      e.put(productDccHeaders); 
    }
    
    if(feUnpacking_){
      e.put(productDigisEB,"ebDigis");
      e.put(productDigisEE,"eeDigis");
      e.put(productInvalidGains,"EcalIntegrityGainErrors");
      e.put(productInvalidGainsSwitch, "EcalIntegrityGainSwitchErrors");
      e.put(productInvalidGainsSwitchStay, "EcalIntegrityGainSwitchStayErrors");
      e.put(productInvalidChIds, "EcalIntegrityChIdErrors");
      e.put(productPnDiodeDigis);
    }
    if(memUnpacking_){
      e.put(productInvalidMemTtIds,"EcalIntegrityMemTtIdErrors");
      e.put(productInvalidMemBlockSizes,"EcalIntegrityMemBlockSizeErrors");
      e.put(productInvalidMemChIds,"EcalIntegrityMemChIdErrors");
      e.put(productInvalidMemGains,"EcalIntegrityMemGainErrors");
    }
    if(srpUnpacking_){
      e.put(productEBSrFlags);
      e.put(productEESrFlags);
    }
    if(tccUnpacking_){
      e.put(productEBTps,"EBTT"); //note change the name
      e.put(productEETps,"EETT"); //note change the name
      e.put(productInvalidTTIds,"EcalIntegrityTTIdErrors");
      e.put(productInvalidBlockLengths,"EcalIntegrityBlockSizeErrors");
    }
  }
  
//if(nevts_>1){   //NUNO
//  double TIME_END = clock(); //NUNO 
//  RUNNING_TIME_ += TIME_END-TIME_START; //NUNO
//}
  
}

EcalRawToDigiDev::~EcalRawToDigiDev() { 
  
  //cout << "EcalDCCUnpackingModule  " << "N events        " << (nevts_-1)<<endl;
  //cout << "EcalDCCUnpackingModule  " << " --- SETUP time " << endl;
  //cout << "EcalDCCUnpackingModule  " << "Time (sys)      " << SETUP_TIME_ << endl;
  //cout << "EcalDCCUnpackingModule  " << "Time in sec.    " << SETUP_TIME_/ CLOCKS_PER_SEC  << endl;
  //cout << "EcalDCCUnpackingModule  " << " --- Per event  " << endl;
  
  //RUNNING_TIME_ = RUNNING_TIME_ / (nevts_-1);
  
  //cout << "EcalDCCUnpackingModule  "<< "Time (sys)      " << RUNNING_TIME_  << endl;
  //cout << "EcalDCCUnpackingModule  "<< "Time in sec.    " << RUNNING_TIME_ / CLOCKS_PER_SEC  << endl;
  
  
  
  if(myMap_      ) delete myMap_;
  if(theUnpacker_) delete theUnpacker_;
  
}
