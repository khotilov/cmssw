#include "DQM/SiPixelMonitorRawData/interface/SiPixelRawDataErrorModule.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQM/SiPixelCommon/interface/SiPixelHistogramId.h"
// Framework
#include "FWCore/ServiceRegistry/interface/Service.h"
// STL
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <sstream>

// Data Formats
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

const int SiPixelRawDataErrorModule::LINK_bits = 6;
const int SiPixelRawDataErrorModule::ROC_bits  = 5;
const int SiPixelRawDataErrorModule::DCOL_bits = 5;
const int SiPixelRawDataErrorModule::PXID_bits = 8;
const int SiPixelRawDataErrorModule::ADC_bits  = 8;
const int SiPixelRawDataErrorModule::DataBit_bits = 1;

const int SiPixelRawDataErrorModule::ADC_shift  = 0;
const int SiPixelRawDataErrorModule::PXID_shift = ADC_shift + ADC_bits;
const int SiPixelRawDataErrorModule::DCOL_shift = PXID_shift + PXID_bits;
const int SiPixelRawDataErrorModule::ROC_shift  = DCOL_shift + DCOL_bits;
const int SiPixelRawDataErrorModule::LINK_shift = ROC_shift + ROC_bits;
const int SiPixelRawDataErrorModule::DB0_shift = 0;
const int SiPixelRawDataErrorModule::DB1_shift = DB0_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB2_shift = DB1_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB3_shift = DB2_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB4_shift = DB3_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB5_shift = DB4_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB6_shift = DB5_shift + DataBit_bits;
const int SiPixelRawDataErrorModule::DB7_shift = DB6_shift + DataBit_bits;

const uint32_t SiPixelRawDataErrorModule::LINK_mask = ~(~uint32_t(0) << LINK_bits);
const uint32_t SiPixelRawDataErrorModule::ROC_mask  = ~(~uint32_t(0) << ROC_bits);
const uint32_t SiPixelRawDataErrorModule::DCOL_mask = ~(~uint32_t(0) << DCOL_bits);
const uint32_t SiPixelRawDataErrorModule::PXID_mask = ~(~uint32_t(0) << PXID_bits);
const uint32_t SiPixelRawDataErrorModule::ADC_mask  = ~(~uint32_t(0) << ADC_bits);
const uint32_t SiPixelRawDataErrorModule::DataBit_mask = ~(~uint32_t(0) << DataBit_bits);

const int SiPixelRawDataErrorModule::TRLRBGN_bits = 32;
const int SiPixelRawDataErrorModule::EVTLGT_bits  = 24;
const int SiPixelRawDataErrorModule::TRLREND_bits = 8;

const int SiPixelRawDataErrorModule::TRLRBGN_shift = 0;
const int SiPixelRawDataErrorModule::EVTLGT_shift  = TRLRBGN_shift + TRLRBGN_bits;
const int SiPixelRawDataErrorModule::TRLREND_shift = EVTLGT_shift + EVTLGT_bits;

const long long SiPixelRawDataErrorModule::TRLREND_mask = ~(~(long long)(0) << TRLREND_bits);
const long long SiPixelRawDataErrorModule::EVTLGT_mask  = ~(~(long long)(0) << EVTLGT_bits);
const long long SiPixelRawDataErrorModule::TRLRBGN_mask = ~(~(long long)(0) << TRLRBGN_bits);
//
// Constructors
//
SiPixelRawDataErrorModule::SiPixelRawDataErrorModule() : 
  id_(0),
  ncols_(416),
  nrows_(160) 
{ 
  _debug_ = false;
  for(int i=0; i!=40; i++) for(int j=0; j!=37; j++){
    FedChNErrArray[i][j]=0; 
    FedChLErrArray[i][j]=0; 
    if(j<15) FedETypeNErrArray[i][j]=0;
  }
}
//
SiPixelRawDataErrorModule::SiPixelRawDataErrorModule(const uint32_t& id) : 
  id_(id),
  ncols_(416),
  nrows_(160)
{ 
  _debug_ = false;
  for(int i=0; i!=40; i++) for(int j=0; j!=37; j++){
    FedChNErrArray[i][j]=0; 
    FedChLErrArray[i][j]=0; 
    if(j<15) FedETypeNErrArray[i][j]=0;
  }
}
//
SiPixelRawDataErrorModule::SiPixelRawDataErrorModule(const uint32_t& id, const int& ncols, const int& nrows) : 
  id_(id),
  ncols_(ncols),
  nrows_(nrows)
{ 
  _debug_ = false;
  for(int i=0; i!=40; i++) for(int j=0; j!=37; j++){
    FedChNErrArray[i][j]=0; 
    FedChLErrArray[i][j]=0; 
    if(j<15) FedETypeNErrArray[i][j]=0;
  }
} 
//
// Destructor
//
SiPixelRawDataErrorModule::~SiPixelRawDataErrorModule() {}
//
// Book histograms for errors with detId
//
void SiPixelRawDataErrorModule::book(const edm::ParameterSet& iConfig, bool reducedSet, int type) {
//std::cout<<"Entering SiPixelRawDataErrorModule::book: "<<std::endl;
  std::string hid;
  // Get collection name and instantiate Histo Id builder
  edm::InputTag src = iConfig.getParameter<edm::InputTag>( "src" );
  bool barrel = DetId(id_).subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);
  bool endcap = DetId(id_).subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap);
  bool isHalfModule = false;
  if(barrel){
    isHalfModule = PixelBarrelName(DetId(id_)).isHalfModule(); 
  }

  // Get DQM interface
  DQMStore* theDMBE = edm::Service<DQMStore>().operator->();

  if(type==0){
    SiPixelHistogramId* theHistogramId = new SiPixelHistogramId( src.label() );
    // Types of errors
    hid = theHistogramId->setHistoId("errorType",id_);
    meErrorType_ = theDMBE->book1D(hid,"Type of errors",15,24.5,39.5);
    meErrorType_->setAxisTitle("Type of errors",1);
    // Number of errors
    hid = theHistogramId->setHistoId("NErrors",id_);
    meNErrors_ = theDMBE->book1D(hid,"Number of errors",10,0.,10.);
    meNErrors_->setAxisTitle("Number of errors",1);
    // For error type 30, the type of problem encoded in the TBM trailer
    // 0 = stack full, 1 = Pre-cal issued, 2 = clear trigger counter, 3 = sync trigger, 
    // 4 = sync trigger error, 5 = reset ROC, 6 = reset TBM, 7 = no token bit pass
    hid = theHistogramId->setHistoId("TBMMessage",id_);
    meTBMMessage_ = theDMBE->book1D(hid,"TBM trailer message",8,-0.5,7.5);
    meTBMMessage_->setAxisTitle("TBM message",1);
    // For error type 30, the type of problem encoded in the FSM bits, 0 = none
    // 1 = FSM errors, 2 = invalid # of ROCs, 3 = data stream too long, 4 = multiple
    hid = theHistogramId->setHistoId("TBMType",id_);
    meTBMType_ = theDMBE->book1D(hid,"State Machine message",5,-0.5,4.5);
    meTBMType_->setAxisTitle("FSM Type",1);
    if(!reducedSet){
      // For error type 31, the event number of the TBM header with the error
      hid = theHistogramId->setHistoId("EvtNbr",id_);
      meEvtNbr_ = theDMBE->bookInt(hid);
      // For error type 36, the invalid ROC number
      hid = theHistogramId->setHistoId("ROCId",id_);
      meROCId_ = theDMBE->bookInt(hid);
      // For error type 37, the invalid dcol values
      hid = theHistogramId->setHistoId("DCOLId",id_);
      meDCOLId_ = theDMBE->bookInt(hid);
      // For error type 37, the invalid ROC values
      hid = theHistogramId->setHistoId("PXId",id_);
      mePXId_ = theDMBE->bookInt(hid);
      // For error type 38, the ROC that is being read out of order
      hid = theHistogramId->setHistoId("ROCNmbr",id_);
      meROCNmbr_ = theDMBE->bookInt(hid);
    }
    delete theHistogramId;
  } else if(type==1 && barrel){
    uint32_t DBladder = PixelBarrelName(DetId(id_)).ladderName();
    char sladder[80]; sprintf(sladder,"Ladder_%02i",DBladder);
    hid = src.label() + "_" + sladder;
    if(isHalfModule) hid += "H";
    else hid += "F";
    // Types of errors
    meErrorTypeLad_ = theDMBE->book1D("errorType_"+hid,"Type of errors",15,24.5,39.5);
    meErrorTypeLad_->setAxisTitle("Type of errors",1);
    // Number of errors
    meNErrorsLad_ = theDMBE->book1D("NErrors_"+hid,"Number of errors",10,0.,10.);
    meNErrorsLad_->setAxisTitle("Number of errors",1);
    // For error type 30, the type of problem encoded in the TBM trailer
    // 0 = stack full, 1 = Pre-cal issued, 2 = clear trigger counter, 3 = sync trigger, 
    // 4 = sync trigger error, 5 = reset ROC, 6 = reset TBM, 7 = no token bit pass
    meTBMMessageLad_ = theDMBE->book1D("TBMMessage_"+hid,"TBM trailer message",8,-0.5,7.5);
    meTBMMessageLad_->setAxisTitle("TBM message",1);
    // For error type 30, the type of problem encoded in the FSM bits, 0 = none
    // 1 = FSM errors, 2 = invalid # of ROCs, 3 = data stream too long, 4 = multiple
    meTBMTypeLad_ = theDMBE->book1D("TBMType_"+hid,"State Machine message",5,-0.5,4.5);
    meTBMTypeLad_->setAxisTitle("FSM Type",1);
    if(!reducedSet){
      // For error type 31, the event number of the TBM header with the error
      meEvtNbrLad_ = theDMBE->bookInt("EvtNbr_"+hid);
      // For error type 36, the invalid ROC number
      meROCIdLad_ = theDMBE->bookInt("ROCId_"+hid);
      // For error type 37, the invalid dcol values
      meDCOLIdLad_ = theDMBE->bookInt("DCOLId_"+hid);
      // For error type 37, the invalid ROC values
      mePXIdLad_ = theDMBE->bookInt("PXId_"+hid);
      // For error type 38, the ROC that is being read out of order
      meROCNmbrLad_ = theDMBE->bookInt("ROCNmbr_"+hid);
    }
  } else if(type==4 && endcap){
    uint32_t blade= PixelEndcapName(DetId(id_)).bladeName();
    char sblade[80]; sprintf(sblade, "Blade_%02i",blade);
    hid = src.label() + "_" + sblade;
    // Types of errors
    meErrorTypeBlade_ = theDMBE->book1D("errorType_"+hid,"Type of errors",15,24.5,39.5);
    meErrorTypeBlade_->setAxisTitle("Type of errors",1);
    // Number of errors
    meNErrorsBlade_ = theDMBE->book1D("NErrors_"+hid,"Number of errors",10,0.,10.);
    meNErrorsBlade_->setAxisTitle("Number of errors",1);
    // For error type 30, the type of problem encoded in the TBM trailer
    // 0 = stack full, 1 = Pre-cal issued, 2 = clear trigger counter, 3 = sync trigger, 
    // 4 = sync trigger error, 5 = reset ROC, 6 = reset TBM, 7 = no token bit pass
    meTBMMessageBlade_ = theDMBE->book1D("TBMMessage_"+hid,"TBM trailer message",8,-0.5,7.5);
    meTBMMessageBlade_->setAxisTitle("TBM message",1);
    // For error type 30, the type of problem encoded in the FSM bits, 0 = none
    // 1 = FSM errors, 2 = invalid # of ROCs, 3 = data stream too long, 4 = multiple
    meTBMTypeBlade_ = theDMBE->book1D("TBMType_"+hid,"State Machine message",5,-0.5,4.5);
    meTBMTypeBlade_->setAxisTitle("FSM Type",1);
    if(!reducedSet){
      // For error type 31, the event number of the TBM header with the error
      meEvtNbrBlade_ = theDMBE->bookInt("EvtNbr_"+hid);
      // For error type 36, the invalid ROC number
      meROCIdBlade_ = theDMBE->bookInt("ROCId_"+hid);
      // For error type 37, the invalid dcol values
      meDCOLIdBlade_ = theDMBE->bookInt("DCOLId_"+hid);
      // For error type 37, the invalid ROC values
      mePXIdBlade_ = theDMBE->bookInt("PXId_"+hid);
      // For error type 38, the ROC that is being read out of order
      meROCNmbrBlade_ = theDMBE->bookInt("ROCNmbr_"+hid);
    }
  }
    
//std::cout<<"...leaving SiPixelRawDataErrorModule::book. "<<std::endl;
}
//
// Book histograms for errors within a FED
//
void SiPixelRawDataErrorModule::bookFED(const edm::ParameterSet& iConfig) {
//std::cout<<"Entering SiPixelRawDataErrorModule::bookFED: "<<std::endl;
  std::string hid;
  // Get collection name and instantiate Histo Id builder
  edm::InputTag src = iConfig.getParameter<edm::InputTag>( "src" );
  SiPixelHistogramId* theHistogramId = new SiPixelHistogramId( src.label() );
  // Get DQM interface
  DQMStore* theDMBE = edm::Service<DQMStore>().operator->();
  // Types of errors
  hid = theHistogramId->setHistoId("errorType",id_);
  meErrorType_ = theDMBE->book1D(hid,"Type of errors",15,24.5,39.5);
  meErrorType_->setAxisTitle("Type of errors",1);
  // Number of errors
  hid = theHistogramId->setHistoId("NErrors",id_);
  meNErrors_ = theDMBE->book1D(hid,"Number of errors",36,0.,36.);
  meNErrors_->setAxisTitle("Number of errors",1);
  // Type of FIFO full (errorType = 28).  FIFO 1 is 1-5 (where fullType = channel of FIFO 1), 
  // fullType = 6 signifies FIFO 2 nearly full, 7 signifies trigger FIFO nearly full, 8 
  // indicates an unexpected result
  hid = theHistogramId->setHistoId("fullType",id_);
  meFullType_ = theDMBE->book1D(hid,"Type of FIFO full",7,0.5,7.5);
  meFullType_->setAxisTitle("FIFO type",1);
  // For errorType = 29, channel with timeout error (0 indicates unexpected result)
  hid = theHistogramId->setHistoId("chanNmbr",id_);
  meChanNmbr_ = theDMBE->book1D(hid,"Timeout Channel",36,1.,37.);
  meChanNmbr_->setAxisTitle("Channel number",1);
  // For error type 30, the type of problem encoded in the TBM trailer
  // 0 = stack full, 1 = Pre-cal issued, 2 = clear trigger counter, 3 = sync trigger, 
  // 4 = sync trigger error, 5 = reset ROC, 6 = reset TBM, 7 = no token bit pass
  hid = theHistogramId->setHistoId("TBMMessage",id_);
  meTBMMessage_ = theDMBE->book1D(hid,"TBM trailer message",8,-0.5,7.5);
  meTBMMessage_->setAxisTitle("TBM message",1);
  // For error type 30, the type of problem encoded in the TBM error trailer 0 = none
  // 1 = data stream too long, 2 = FSM errors, 3 = invalid # of ROCs, 4 = multiple
  hid = theHistogramId->setHistoId("TBMType",id_);
  meTBMType_ = theDMBE->book1D(hid,"Type of TBM trailer",5,-0.5,4.5);
  meTBMType_->setAxisTitle("TBM Type",1);
  // For error type 31, the event number of the TBM header with the error
  hid = theHistogramId->setHistoId("EvtNbr",id_);
  meEvtNbr_ = theDMBE->bookInt(hid);
  // For errorType = 34, datastream size according to error word
  hid = theHistogramId->setHistoId("evtSize",id_);
  meEvtSize_ = theDMBE->bookInt(hid);
  // For errorType = 35, the bad channel number
  hid = theHistogramId->setHistoId("linkId",id_);
  meLinkId_ = theDMBE->book1D(hid,"Invalid Channel Number",36,1.,37.);
  meLinkId_->setAxisTitle("FED Channel number",1);
  // For error type 36, the invalid ROC number matched to FED Channel
  hid = theHistogramId->setHistoId("Type36Hitmap",id_);
  meInvROC_ = theDMBE->book2D(hid,"Invalid ROC by Channel Number",36,1.,37.,16,0.,16.);
  meInvROC_->setAxisTitle("FED Channel Number",1);
  meInvROC_->setAxisTitle("ROC Number",2);
  
  //for(int i=0; i!=40; i++)
  for(int j=0; j!=37; j++){
    //if(static_cast<int>(id_) == i){
      std::stringstream temp; temp << j;
      hid = "FedChNErrArray_" + temp.str();
      meFedChNErrArray_[j] = theDMBE->bookInt(hid);
      hid = "FedChLErrArray_" + temp.str();
      meFedChLErrArray_[j] = theDMBE->bookInt(hid);
      hid = "FedETypeNErrArray_" + temp.str();
      if(j<15) meFedETypeNErrArray_[j] = theDMBE->bookInt(hid);
    //}
  }
  delete theHistogramId;
//std::cout<<"...leaving SiPixelRawDataErrorModule::bookFED. "<<std::endl;
}
//
// Fill histograms
//
int SiPixelRawDataErrorModule::fill(const edm::DetSetVector<SiPixelRawDataError>& input, bool reducedSet, bool modon, bool ladon, bool bladeon) {
//std::cout<<"Entering SiPixelRawDataErrorModule::fill: "<<std::endl;
  bool barrel = DetId(id_).subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel);
  bool endcap = DetId(id_).subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap);
  bool isHalfModule = false;
  uint32_t DBladder = 0;
  if(barrel){
    isHalfModule = PixelBarrelName(DetId(id_)).isHalfModule(); 
    DBladder = PixelBarrelName(DetId(id_)).ladderName();
//    int FedId = DetId(id_).getFedId();                  // FED the error came from

  }
  
  unsigned int numberOfErrors = 0;  
  unsigned int numberOfSeriousErrors = 0;
  
  edm::DetSetVector<SiPixelRawDataError>::const_iterator isearch = input.find(id_); // search  errors of detid
  //edm::DetSetVector<SiPixelRawDataError>::const_iterator isearch = input.find(0xffffffff); // search  errors of detid
  
  if( isearch != input.end() ) {  // Not at empty iterator
    // Look at errors now
    edm::DetSet<SiPixelRawDataError>::const_iterator  di;
    for(di = isearch->data.begin(); di != isearch->data.end(); di++) {
      numberOfErrors++;
      int FedId = di->getFedId();                  // FED the error came from
      int chanNmbr = -1;
      int errorType = di->getType();               // type of error
      int TBMType;
      if(modon){
        (meErrorType_)->Fill((int)errorType);
	const int LINK_bits = 6;
	const int LINK_shift = 26;
	const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
        if(!reducedSet && (errorType == 32 || errorType == 33 || errorType == 34) ) {
	  long long errorWord = di->getWord64();     // for 64-bit error words
	  chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	  if(errorType == 34) {
	    int evtSize = (errorWord >> EVTLGT_shift) & EVTLGT_mask;
	    (meEvtSize_)->Fill((int)evtSize); 
	  }
        } else {
	  uint32_t errorWord = di->getWord32();      // for 32-bit error words
	  chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	  switch(errorType) {  // fill in the appropriate monitorables based on the information stored in the error word
	  case(30) : {
	    int T0 = (errorWord >> DB0_shift) & DataBit_mask;
	    int T1 = (errorWord >> DB1_shift) & DataBit_mask;
	    int T2 = (errorWord >> DB2_shift) & DataBit_mask;
	    int T3 = (errorWord >> DB3_shift) & DataBit_mask;
	    int T4 = (errorWord >> DB4_shift) & DataBit_mask;
	    int T5 = (errorWord >> DB5_shift) & DataBit_mask;
	    int T6 = (errorWord >> DB6_shift) & DataBit_mask;
	    int T7 = (errorWord >> DB7_shift) & DataBit_mask;
	    int TBMMessage;
	    if (T0==1) { TBMMessage=0; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T1==1) { TBMMessage=1; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T2==1) { TBMMessage=2; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T3==1) { TBMMessage=3; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T4==1) { TBMMessage=4; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T5==1) { TBMMessage=5; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T6==1) { TBMMessage=6; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T7==1) { TBMMessage=7; (meTBMMessage_)->Fill((int)TBMMessage); }
	    int StateMach_bits      = 4;
	    int StateMach_shift     = 8;
	    uint32_t StateMach_mask = ~(~uint32_t(0) << StateMach_bits);
	    int StateMach = (errorWord >> StateMach_shift) & StateMach_mask;
	    switch(StateMach) {
	    case(0) : {
	      TBMType = 0;
	      break; }
	    case(1) : {
	      TBMType = 1;
	      break; }
	    case(2) : case(4) : case(6) : {
	      TBMType = 2;
	      break; }
	    case(8) : {
	      TBMType = 3;
	      break; }
	    default : TBMType = 4;
	    };
	    (meTBMType_)->Fill((int)TBMType);
	    break; }
	  case(31) : {
            if(!reducedSet){
	      int evtNbr = (errorWord >> ADC_shift) & ADC_mask;
	      (meEvtNbr_)->Fill((int)evtNbr);
	    }
	    break; }
	  case(36) : {
            if(!reducedSet){
	      int ROCId = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCId_)->Fill((int)ROCId);
	    }
	    break; }
	  case(37) : {
            if(!reducedSet){
	      int DCOLId = (errorWord >> DCOL_shift) & DCOL_mask;
	      (meDCOLId_)->Fill((int)DCOLId);
	      int PXId = (errorWord >> PXID_shift) & PXID_mask;
	      (mePXId_)->Fill((int)PXId);
	    }
	    break; }
	  case(38) : {
            if(!reducedSet){
	      int ROCNmbr = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCNmbr_)->Fill((int)ROCNmbr);
	    }
	    break; }
	  default : break;
	  };
	}  
	//std::cout<<"INFO: "<<FedId<<" , "<<chanNmbr<<" , "<<errorType<<std::endl;
	if(!(errorType==30)){
	  numberOfSeriousErrors++;
	  FedChNErrArray[FedId][chanNmbr]++;
	  FedChLErrArray[FedId][chanNmbr] = errorType;
	  FedETypeNErrArray[FedId][errorType-25]++;
	}
      }
      
      if(ladon && barrel){
        (meErrorTypeLad_)->Fill((int)errorType);
        if(!reducedSet && (errorType == 32 || errorType == 33 || errorType == 34) ) {
	  long long errorWord = di->getWord64();     // for 64-bit error words
	  if(errorType == 34) {
	    int evtSize = (errorWord >> EVTLGT_shift) & EVTLGT_mask;
	    (meEvtSizeLad_)->Fill((int)evtSize); 
	  }
        } else {
	  uint32_t errorWord = di->getWord32();      // for 32-bit error words
	  switch(errorType) {  // fill in the appropriate monitorables based on the information stored in the error word
	  case(30) : {
	    int T0 = (errorWord >> DB0_shift) & DataBit_mask;
	    int T1 = (errorWord >> DB1_shift) & DataBit_mask;
	    int T2 = (errorWord >> DB2_shift) & DataBit_mask;
	    int T3 = (errorWord >> DB3_shift) & DataBit_mask;
	    int T4 = (errorWord >> DB4_shift) & DataBit_mask;
	    int T5 = (errorWord >> DB5_shift) & DataBit_mask;
	    int T6 = (errorWord >> DB6_shift) & DataBit_mask;
	    int T7 = (errorWord >> DB7_shift) & DataBit_mask;
	    int TBMMessage;
	    if (T0==1) { TBMMessage=0; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T1==1) { TBMMessage=1; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T2==1) { TBMMessage=2; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T3==1) { TBMMessage=3; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T4==1) { TBMMessage=4; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T5==1) { TBMMessage=5; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T6==1) { TBMMessage=6; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    if (T7==1) { TBMMessage=7; (meTBMMessageLad_)->Fill((int)TBMMessage); }
	    int StateMach_bits      = 4;
	    int StateMach_shift     = 8;
	    uint32_t StateMach_mask = ~(~uint32_t(0) << StateMach_bits);
	    int StateMach = (errorWord >> StateMach_shift) & StateMach_mask;
	    int TBMType;
	    switch(StateMach) {
	    case(0) : {
	      TBMType = 0;
	      break; }
	    case(1) : {
	      TBMType = 1;
	      break; }
	    case(2) : case(4) : case(6) : {
	      TBMType = 2;
	      break; }
	    case(8) : {
	      TBMType = 3;
	      break; }
	    default : TBMType = 4;
	    };
	    (meTBMTypeLad_)->Fill((int)TBMType);
	    break; }
	  case(31) : {
            if(!reducedSet){
	      int evtNbr = (errorWord >> ADC_shift) & ADC_mask;
	      (meEvtNbrLad_)->Fill((int)evtNbr);
	    }
	    break; }
	  case(36) : {
            if(!reducedSet){
	      int ROCId = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCIdLad_)->Fill((int)ROCId);
	    }
	    break; }
	  case(37) : {
            if(!reducedSet){
	      int DCOLId = (errorWord >> DCOL_shift) & DCOL_mask;
	      (meDCOLIdLad_)->Fill((int)DCOLId);
	      int PXId = (errorWord >> PXID_shift) & PXID_mask;
	      (mePXIdLad_)->Fill((int)PXId);
	    }
	    break; }
	  case(38) : {
            if(!reducedSet){
	      int ROCNmbr = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCNmbrLad_)->Fill((int)ROCNmbr);
	    }
	    break; }
	  default : break;
	  };
	}
      }

      if(bladeon && endcap){
        (meErrorTypeBlade_)->Fill((int)errorType);
        if(!reducedSet && (errorType == 32 || errorType == 33 || errorType == 34) ) {
	  long long errorWord = di->getWord64();     // for 64-bit error words
	  if(errorType == 34) {
	    int evtSize = (errorWord >> EVTLGT_shift) & EVTLGT_mask;
	    (meEvtSizeBlade_)->Fill((int)evtSize); 
	  }
        } else {
	  uint32_t errorWord = di->getWord32();      // for 32-bit error words
	  switch(errorType) {  // fill in the appropriate monitorables based on the information stored in the error word
	  case(30) : {
	    int T0 = (errorWord >> DB0_shift) & DataBit_mask;
	    int T1 = (errorWord >> DB1_shift) & DataBit_mask;
	    int T2 = (errorWord >> DB2_shift) & DataBit_mask;
	    int T3 = (errorWord >> DB3_shift) & DataBit_mask;
	    int T4 = (errorWord >> DB4_shift) & DataBit_mask;
	    int T5 = (errorWord >> DB5_shift) & DataBit_mask;
	    int T6 = (errorWord >> DB6_shift) & DataBit_mask;
	    int T7 = (errorWord >> DB7_shift) & DataBit_mask;
	    int TBMMessage;
	    if (T0==1) { TBMMessage=0; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T1==1) { TBMMessage=1; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T2==1) { TBMMessage=2; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T3==1) { TBMMessage=3; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T4==1) { TBMMessage=4; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T5==1) { TBMMessage=5; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T6==1) { TBMMessage=6; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    if (T7==1) { TBMMessage=7; (meTBMMessageBlade_)->Fill((int)TBMMessage); }
	    int StateMach_bits      = 4;
	    int StateMach_shift     = 8;
	    uint32_t StateMach_mask = ~(~uint32_t(0) << StateMach_bits);
	    int StateMach = (errorWord >> StateMach_shift) & StateMach_mask;
	    int TBMType;
	    switch(StateMach) {
	    case(0) : {
	      TBMType = 0;
	      break; }
	    case(1) : {
	      TBMType = 1;
	      break; }
	    case(2) : case(4) : case(6) : {
	      TBMType = 2;
	      break; }
	    case(8) : {
	      TBMType = 3;
	      break; }
	    default : TBMType = 4;
	    };
	    (meTBMTypeBlade_)->Fill((int)TBMType);
	    break; }
	  case(31) : {
            if(!reducedSet){
	      int evtNbr = (errorWord >> ADC_shift) & ADC_mask;
	      (meEvtNbrBlade_)->Fill((int)evtNbr);
	    }
	    break; }
	  case(36) : {
            if(!reducedSet){
	      int ROCId = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCIdBlade_)->Fill((int)ROCId);
	    }
	    break; }
	  case(37) : {
            if(!reducedSet){
	      int DCOLId = (errorWord >> DCOL_shift) & DCOL_mask;
	      (meDCOLIdBlade_)->Fill((int)DCOLId);
	      int PXId = (errorWord >> PXID_shift) & PXID_mask;
	      (mePXIdBlade_)->Fill((int)PXId);
	    }
	    break; }
	  case(38) : {
            if(!reducedSet){
	      int ROCNmbr = (errorWord >> ROC_shift) & ROC_mask;
	      (meROCNmbrBlade_)->Fill((int)ROCNmbr);
	    }
	    break; }
	  default : break;
	  };
	}
      }
    }//for loop
    
    if(modon && numberOfErrors>0) (meNErrors_)->Fill((float)numberOfErrors);
    if(ladon && barrel && numberOfErrors>0) (meNErrorsLad_)->Fill((float)numberOfErrors);
    if(bladeon && endcap && numberOfErrors>0) (meNErrorsBlade_)->Fill((float)numberOfErrors);
    //std::cout<<"number of errors="<<numberOfErrors<<std::endl;
    
    std::string hid;
    // Get DQM interface
    DQMStore* theDMBE = edm::Service<DQMStore>().operator->();
    for(int i=0; i!=40; i++){
      //std::cout<<"FED: "<<i<<" , "<<static_cast<int>(id_)<<std::endl;
      //if(static_cast<int>(id_) == i){
	for(int j=0; j!=37; j++){
          if(FedChNErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedChNErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
            meFedChNErrArray_[j] = theDMBE->get(hid); 
	    if(meFedChNErrArray_[j]) meFedChNErrArray_[j]->Fill(FedChNErrArray[i][j]); 
          }
          if(FedChLErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedChLErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
            meFedChLErrArray_[j] = theDMBE->get(hid); 
	    if(meFedChLErrArray_[j]) meFedChLErrArray_[j]->Fill(FedChLErrArray[i][j]); 
          }
          if(j<=14 && FedETypeNErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedETypeNErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
	    meFedETypeNErrArray_[j] = theDMBE->get(hid); 
          if(meFedETypeNErrArray_[j]) meFedETypeNErrArray_[j]->Fill(FedETypeNErrArray[i][j]); 
          }
        }//end for
      //}//end if
    }//end if not an empty iterator
  }
  
  //std::cout<<"number of detector units="<<numberOfDetUnits<<std::endl;
//std::cout<<"...leaving SiPixelRawDataErrorModule::fill. "<<std::endl;
  return numberOfSeriousErrors;
}

int SiPixelRawDataErrorModule::fillFED(const edm::DetSetVector<SiPixelRawDataError>& input) {
  //std::cout<<"Entering   SiPixelRawDataErrorModule::fillFED: "<<static_cast<int>(id_)<<std::endl;
  unsigned int numberOfErrors = 0;
  unsigned int numberOfSeriousErrors = 0;
  
  edm::DetSetVector<SiPixelRawDataError>::const_iterator isearch = input.find(0xffffffff); // search  errors of detid
  if( isearch != input.end() ) {  // Not an empty iterator
    //std::cout<<"FED has some errors!"<<std::endl;    
    
    // Look at FED errors now
    int TBMType;
    edm::DetSet<SiPixelRawDataError>::const_iterator  di;
    for(di = isearch->data.begin(); di != isearch->data.end(); di++) {
      //std::cout<<"Looping over errors now!"<<std::endl;
      int FedId = di->getFedId();                  // FED the error came from
      int chanNmbr = -1;
      int errorType = 0;
      if(FedId==static_cast<int>(id_)) {
	numberOfErrors++;
	errorType = di->getType();               // type of error
	(meErrorType_)->Fill((int)errorType);
	int TBMMessage = -1;
	if((errorType == 32)||(errorType == 33)||(errorType == 34)) {
	  long long errorWord = di->getWord64();     // for 64-bit error words
	  switch(errorType) {  // fill in the appropriate monitorables based on the information stored in the error word
	  case(32) : {
	    chanNmbr = 0;
	    break; }
	  case(33) : {
	    chanNmbr = 0;
	    break; }
	  case(34) : {
	    chanNmbr = 0;
	    int evtSize = (errorWord >> EVTLGT_shift) & EVTLGT_mask;
	    (meEvtSize_)->Fill((int)evtSize); 
	    break; }
	  default : break;
	  };
	} else {
	  uint32_t errorWord = di->getWord32();      // for 32-bit error words
	  switch(errorType) {  // fill in the appropriate monitorables based on the information stored in the error word
	  case(25) : {
	    chanNmbr = 0;
	    break; }
	  case(28) : {
	    int NFa = (errorWord >> DB0_shift) & DataBit_mask;
	    int NFb = (errorWord >> DB1_shift) & DataBit_mask;
	    int NFc = (errorWord >> DB2_shift) & DataBit_mask;
	    int NFd = (errorWord >> DB3_shift) & DataBit_mask;
	    int NFe = (errorWord >> DB4_shift) & DataBit_mask;
	    int NF2 = (errorWord >> DB6_shift) & DataBit_mask;
	    int L1A = (errorWord >> DB7_shift) & DataBit_mask;
	    int fullType = 0;
	    if (NFa==1) fullType = 1; (meFullType_)->Fill((int)fullType);
	    if (NFb==1) fullType = 2; (meFullType_)->Fill((int)fullType);
	    if (NFc==1) fullType = 3; (meFullType_)->Fill((int)fullType);
	    if (NFd==1) fullType = 4; (meFullType_)->Fill((int)fullType);
	    if (NFe==1) fullType = 5; (meFullType_)->Fill((int)fullType);
	    if (NF2==1) fullType = 6; (meFullType_)->Fill((int)fullType);
	    if (L1A==1) fullType = 7; (meFullType_)->Fill((int)fullType);
	    chanNmbr = 0;
	    break; }
	  case(29) : {
	    int CH1 = (errorWord >> DB0_shift) & DataBit_mask;
	    int CH2 = (errorWord >> DB1_shift) & DataBit_mask;
	    int CH3 = (errorWord >> DB2_shift) & DataBit_mask;
	    int CH4 = (errorWord >> DB3_shift) & DataBit_mask;
	    int CH5 = (errorWord >> DB4_shift) & DataBit_mask;
	    int BLOCK_bits      = 3;
	    int BLOCK_shift     = 8;
	    uint32_t BLOCK_mask = ~(~uint32_t(0) << BLOCK_bits);
	    int BLOCK = (errorWord >> BLOCK_shift) & BLOCK_mask;
	    int localCH = 1*CH1+2*CH2+3*CH3+4*CH4+5*CH5;
	    if (BLOCK%2==0) chanNmbr=(BLOCK/2)*9+localCH;
	    else chanNmbr = (BLOCK-1)/2+4+localCH;
	    if ((chanNmbr<1)||(chanNmbr>36)) chanNmbr=0;  // signifies unexpected result
	    (meChanNmbr_)->Fill((int)chanNmbr);
	    break; }
	  case(30) : {
	    int T0 = (errorWord >> DB0_shift) & DataBit_mask;
	    int T1 = (errorWord >> DB1_shift) & DataBit_mask;
	    int T2 = (errorWord >> DB2_shift) & DataBit_mask;
	    int T3 = (errorWord >> DB3_shift) & DataBit_mask;
	    int T4 = (errorWord >> DB4_shift) & DataBit_mask;
	    int T5 = (errorWord >> DB5_shift) & DataBit_mask;
	    int T6 = (errorWord >> DB6_shift) & DataBit_mask;
	    int T7 = (errorWord >> DB7_shift) & DataBit_mask;
	    if (T0==1) { TBMMessage=0; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T1==1) { TBMMessage=1; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T2==1) { TBMMessage=2; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T3==1) { TBMMessage=3; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T4==1) { TBMMessage=4; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T5==1) { TBMMessage=5; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T6==1) { TBMMessage=6; (meTBMMessage_)->Fill((int)TBMMessage); }
	    if (T7==1) { TBMMessage=7; (meTBMMessage_)->Fill((int)TBMMessage); }
	    int StateMach_bits      = 4;
	    int StateMach_shift     = 8;
	    uint32_t StateMach_mask = ~(~uint32_t(0) << StateMach_bits);
	    int StateMach = (errorWord >> StateMach_shift) & StateMach_mask;
	    switch(StateMach) {
	    case(0) : {
	      TBMType = 0;
	      break; }
	    case(1) : {
	      TBMType = 1;
	      break; }
	    case(2) : case(4) : case(6) : {
	      TBMType = 2;
	      break; }
	    case(8) : {
	      TBMType = 3;
	      break; }
	    default : TBMType = 4;
	    };
	    (meTBMType_)->Fill((int)TBMType);
	    const int LINK_bits = 6;
	    const int LINK_shift = 26;
	    const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
	    chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	    break; }
	  case(31) : {
	    int evtNbr = (errorWord >> ADC_shift) & ADC_mask;
	    (meEvtNbr_)->Fill((int)evtNbr);
	    const int LINK_bits = 6;
	    const int LINK_shift = 26;
	    const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
	    chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	    break; }
	  case(35) : {
	    int linkId = (errorWord >> LINK_shift) & LINK_mask;
	    (meLinkId_)->Fill((float)linkId);
	    break; }
	  case(36) : {
	    int ROCId = (errorWord >> ROC_shift) & ROC_mask;
	    int ChanId = (errorWord >> LINK_shift) & LINK_mask;
	    (meInvROC_)->Fill((float)ChanId,(float)ROCId);
	    const int LINK_bits = 6;
	    const int LINK_shift = 26;
	    const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
	    chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	    break; }
	  case(37) : {
	    const int LINK_bits = 6;
	    const int LINK_shift = 26;
	    const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
	    chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	    break; }
	  case(38) : {
	    const int LINK_bits = 6;
	    const int LINK_shift = 26;
	    const uint32_t LINK_mask = ~(~(uint32_t)(0) << LINK_bits);
	    chanNmbr = (errorWord >> LINK_shift) & LINK_mask;
	    break; }
	  case(39) : {
	    chanNmbr = 0;
	    break; }
	  default : break;
	  };
	}// end if errorType
	if(!(errorType==30)){
	  numberOfSeriousErrors++;
	  FedChNErrArray[FedId][chanNmbr]++;
	  FedChLErrArray[FedId][chanNmbr] = errorType;
	  FedETypeNErrArray[FedId][errorType-25]++;
	}
      }// end if FedId
    }// end for loop over all error words
    
    if(numberOfErrors>0) (meNErrors_)->Fill((float)numberOfErrors);
    //std::cout<<"number of errors="<<numberOfErrors<<std::endl;
    std::string hid;
    // Get DQM interface
    DQMStore* theDMBE = edm::Service<DQMStore>().operator->();
    for(int i=0; i!=40; i++){
      if(static_cast<int>(id_) == i){
	for(int j=0; j!=37; j++){
          if(FedChNErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedChNErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
            meFedChNErrArray_[j] = theDMBE->get(hid); 
	    if(meFedChNErrArray_[j]) meFedChNErrArray_[j]->Fill(FedChNErrArray[i][j]); 
          }
          if(FedChLErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedChLErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
            meFedChLErrArray_[j] = theDMBE->get(hid); 
	    if(meFedChLErrArray_[j]) meFedChLErrArray_[j]->Fill(FedChLErrArray[i][j]); 
          }
          if(j<=14 && FedETypeNErrArray[i][j]>0){
            static const char fmt[] = "Pixel/AdditionalPixelErrors/FED_%d/FedETypeNErrArray_%d";
            char buf[sizeof(fmt) + 2*32]; // 32 digits is enough for up to 2^105 + sign.
            sprintf(buf, fmt, i, j);
            hid = buf;
	    meFedETypeNErrArray_[j] = theDMBE->get(hid); 
          if(meFedETypeNErrArray_[j]) meFedETypeNErrArray_[j]->Fill(FedETypeNErrArray[i][j]); 
          }
        }//end for
      }//end if
    }//end for
  }// end if not an empty iterator
  
//  std::cout<<"number of detector units="<<numberOfDetUnits<<std::endl;
//std::cout<<"...leaving   SiPixelRawDataErrorModule::fillFED. "<<std::endl;
  return numberOfSeriousErrors;
}
