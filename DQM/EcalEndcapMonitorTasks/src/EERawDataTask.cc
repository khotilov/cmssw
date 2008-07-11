#include <iostream>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/FEDRawData/src/fed_header.h"

#include <DQM/EcalCommon/interface/Numbers.h>

#include "DQM/EcalEndcapMonitorTasks/interface/EERawDataTask.h"

using namespace cms;
using namespace edm;
using namespace std;

EERawDataTask::EERawDataTask(const ParameterSet& ps) {

  init_ = false;

  dqmStore_ = Service<DQMStore>().operator->();

  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  mergeRuns_ = ps.getUntrackedParameter<bool>("mergeRuns", false);

  FEDRawDataCollection_ = ps.getParameter<edm::InputTag>("FEDRawDataCollection");
  EcalRawDataCollection_ = ps.getParameter<edm::InputTag>("EcalRawDataCollection");

  meEECRCErrors_ = 0;
  meEERunNumberErrors_ = 0;
  meEEL1AErrors_ = 0;
  meEEOrbitNumberErrors_ = 0;
  meEEBunchCrossingErrors_ = 0;
  meEETriggerTypeErrors_ = 0;

}

EERawDataTask::~EERawDataTask() {
}

void EERawDataTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EERawDataTask");
    dqmStore_->rmdir(prefixME_ + "/EERawDataTask");
  }

  Numbers::initGeometry(c, false);

}

void EERawDataTask::beginRun(const Run& r, const EventSetup& c) {

  if ( ! mergeRuns_ ) this->reset();

}

void EERawDataTask::endRun(const Run& r, const EventSetup& c) {

}

void EERawDataTask::reset(void) {

  if ( meEECRCErrors_ ) meEECRCErrors_->Reset();
  if ( meEERunNumberErrors_ ) meEERunNumberErrors_->Reset();
  if ( meEEL1AErrors_ ) meEEL1AErrors_->Reset();
  if ( meEEOrbitNumberErrors_ ) meEEOrbitNumberErrors_->Reset();
  if ( meEEBunchCrossingErrors_ ) meEEBunchCrossingErrors_->Reset();
  if ( meEETriggerTypeErrors_ ) meEETriggerTypeErrors_->Reset();

}

void EERawDataTask::setup(void){

  init_ = true;

  char histo[200];

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EERawDataTask");

    sprintf(histo, "EERDT CRC errors");
    meEECRCErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19);
    for (int i = 0; i < 18; i++) {
      meEECRCErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }

    sprintf(histo, "EERDT run number errors");
    meEERunNumberErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19); 
    for (int i = 0; i < 18; i++) {
      meEERunNumberErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }

    sprintf(histo, "EERDT L1A errors");
    meEEL1AErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19); 
    for (int i = 0; i < 18; i++) {
      meEEL1AErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }

    sprintf(histo, "EERDT orbit number errors");
    meEEOrbitNumberErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19); 
    for (int i = 0; i < 18; i++) {
      meEEOrbitNumberErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }

    sprintf(histo, "EERDT bunch crossing errors");
    meEEBunchCrossingErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19); 
    for (int i = 0; i < 18; i++) {
      meEEBunchCrossingErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }

    sprintf(histo, "EERDT trigger type errors");
    meEETriggerTypeErrors_ = dqmStore_->book1D(histo, histo, 18, 1, 19); 
    for (int i = 0; i < 18; i++) {
      meEETriggerTypeErrors_->setBinLabel(i+1, Numbers::sEB(i+1).c_str(), 1);
    }
  }

}

void EERawDataTask::cleanup(void){

  if ( ! init_ ) return;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EERawDataTask");

    if ( meEECRCErrors_ ) dqmStore_->removeElement( meEECRCErrors_->getName() );
    meEECRCErrors_ = 0;

    if ( meEERunNumberErrors_ ) dqmStore_->removeElement( meEERunNumberErrors_->getName() );
    meEERunNumberErrors_ = 0;

    if ( meEEL1AErrors_ ) dqmStore_->removeElement( meEEL1AErrors_->getName() );
    meEEL1AErrors_ = 0;

    if ( meEEOrbitNumberErrors_ ) dqmStore_->removeElement( meEEOrbitNumberErrors_->getName() );
    meEEOrbitNumberErrors_ = 0;

    if ( meEEBunchCrossingErrors_ ) dqmStore_->removeElement( meEEBunchCrossingErrors_->getName() );
    meEEBunchCrossingErrors_ = 0;

    if ( meEETriggerTypeErrors_ ) dqmStore_->removeElement( meEETriggerTypeErrors_->getName() );
    meEETriggerTypeErrors_ = 0;

  }

  init_ = false;

}

void EERawDataTask::endJob(void) {

  LogInfo("EERawDataTask") << "analyzed " << ievt_ << " events";

  if ( enableCleanup_ ) this->cleanup();

}

void EERawDataTask::analyze(const Event& e, const EventSetup& c){

  if ( ! init_ ) this->setup();

  ievt_++;

  int GT_L1A=0, GT_OrbitNumber=0, GT_BunchCrossing=0, GT_TriggerType=0;

  edm::Handle<FEDRawDataCollection> allFedRawData;

  if ( e.getByLabel(FEDRawDataCollection_, allFedRawData) ) {

    // GT FED data
    const FEDRawData& gtFedData = allFedRawData->FEDData(812);
    
    int length = gtFedData.size()/sizeof(uint64_t);
    
    std::cout << "GT FED length = " << length << endl;
    
    if ( length > 0 ) {
      
      FEDHeader header(gtFedData.data());
      
      GT_L1A = header.lvl1ID();
      GT_BunchCrossing = header.bxID();
      GT_TriggerType = header.triggerType();
      
      //      uint64_t * pData = (uint64_t *)(gtFedData.data());
      /// FIXME: how to get the orbit from the GT?
      GT_OrbitNumber = 0;

    }

    // ECAL endcap FEDs
    int EEFirstFED[2];
    EEFirstFED[0] = 601; // EE-
    EEFirstFED[1] = 646; // EE+
    for(int zside=0; zside<2; zside++) {

      int firstFedOnSide=EEFirstFED[zside];

      for(int i=0; i<9; i++) {
      
	const FEDRawData& fedData = allFedRawData->FEDData(firstFedOnSide+i);
	
	int length = fedData.size()/sizeof(uint64_t);
	
	if ( length > 0 ) {
	  
	  uint64_t * pData = (uint64_t *)(fedData.data());
	  uint64_t * fedTrailer = pData + (length - 1);
	  bool crcError = (*fedTrailer >> 2 ) & 0x1; 
	  
	  if (crcError) meEECRCErrors_->Fill( i+1 );
	  
	}
	
      }

    }

  } else {
    LogWarning("EERawDataTask") << FEDRawDataCollection_ << " not available";
  }


  Handle<EcalRawDataCollection> dcchs;

  if ( e.getByLabel(EcalRawDataCollection_, dcchs) ) {
    
    if ( dcchs.isValid() ) {
      
      for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin(); dcchItr != dcchs->end(); ++dcchItr ) {

	EcalDCCHeaderBlock dcch = (*dcchItr);
	
	if ( Numbers::subDet( dcch ) != EcalEndcap ) continue;
	
	int ism = Numbers::iSM( dcch, EcalEndcap );
	float xism = ism+0.5;

	int ECALDCC_L1A = dcch.getLV1();
	int ECALDCC_OrbitNumber = dcch.getOrbit();
	int ECALDCC_BunchCrossing = dcch.getBX();
	int ECALDCC_TriggerType = dcch.getBasicTriggerType();

	if ( GT_L1A != ECALDCC_L1A ) meEEL1AErrors_->Fill( xism );

	if ( GT_OrbitNumber != ECALDCC_OrbitNumber ) meEEOrbitNumberErrors_->Fill ( xism );
	
	if ( GT_BunchCrossing != ECALDCC_BunchCrossing ) meEEBunchCrossingErrors_->Fill( xism );

	if ( GT_TriggerType != ECALDCC_TriggerType ) meEETriggerTypeErrors_->Fill ( xism );

      }
      
    }

  } else {
    LogWarning("EERawDataTask") << EcalRawDataCollection_ << " not available";
  }

}

