/*
 * \file EELedTask.cc
 *
 * $Date: 2009/10/26 17:33:51 $
 * $Revision: 1.54 $
 * \author G. Della Ricca
 *
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <DQM/EcalCommon/interface/Numbers.h>
#include <DQM/EcalCommon/interface/NumbersPn.h>

#include <DQM/EcalEndcapMonitorTasks/interface/EELedTask.h>

using namespace cms;
using namespace edm;
using namespace std;

EELedTask::EELedTask(const ParameterSet& ps){

  init_ = false;

  dqmStore_ = Service<DQMStore>().operator->();

  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  mergeRuns_ = ps.getUntrackedParameter<bool>("mergeRuns", false);

  EcalRawDataCollection_ = ps.getParameter<edm::InputTag>("EcalRawDataCollection");
  EEDigiCollection_ = ps.getParameter<edm::InputTag>("EEDigiCollection");
  EcalPnDiodeDigiCollection_ = ps.getParameter<edm::InputTag>("EcalPnDiodeDigiCollection");
  EcalUncalibratedRecHitCollection_ = ps.getParameter<edm::InputTag>("EcalUncalibratedRecHitCollection");

  // vector of enabled wavelengths (Default to all 2)
  ledWavelengths_.reserve(2);
  for ( unsigned int i = 1; i <= 2; i++ ) ledWavelengths_.push_back(i);
  ledWavelengths_ = ps.getUntrackedParameter<vector<int> >("ledWavelengths", ledWavelengths_);

  for (int i = 0; i < 18; i++) {
    meShapeMapL1_[i] = 0;
    meAmplMapL1_[i] = 0;
    meTimeMapL1_[i] = 0;
    meAmplPNMapL1_[i] = 0;
    mePnAmplMapG01L1_[i] = 0;
    mePnPedMapG01L1_[i] = 0;
    mePnAmplMapG16L1_[i] = 0;
    mePnPedMapG16L1_[i] = 0;

    meShapeMapL2_[i] = 0;
    meAmplMapL2_[i] = 0;
    meTimeMapL2_[i] = 0;
    meAmplPNMapL2_[i] = 0;
    mePnAmplMapG01L2_[i] = 0;
    mePnPedMapG01L2_[i] = 0;
    mePnAmplMapG16L2_[i] = 0;
    mePnPedMapG16L2_[i] = 0;
  }

}

EELedTask::~EELedTask(){

}

void EELedTask::beginJob(void){

  ievt_ = 0;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask");
    dqmStore_->rmdir(prefixME_ + "/EELedTask");
  }

}

void EELedTask::beginRun(const Run& r, const EventSetup& c) {

  Numbers::initGeometry(c, false);

  if ( ! mergeRuns_ ) this->reset();

}

void EELedTask::endRun(const Run& r, const EventSetup& c) {

  for (int i = 0; i < 18; i++) {
    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {
      if ( meShapeMapL1_[i] )  meShapeMapL1_[i]->Reset();
      if ( meAmplMapL1_[i] ) meAmplMapL1_[i]->Reset();
      if ( meTimeMapL1_[i] ) meTimeMapL1_[i]->Reset();
      if ( meAmplPNMapL1_[i] ) meAmplPNMapL1_[i]->Reset();
    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {
      if ( meShapeMapL2_[i] )  meShapeMapL2_[i]->Reset();
      if ( meAmplMapL2_[i] ) meAmplMapL2_[i]->Reset();
      if ( meTimeMapL2_[i] ) meTimeMapL2_[i]->Reset();
      if ( meAmplPNMapL2_[i] ) meAmplPNMapL2_[i]->Reset();
    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {
      if ( mePnAmplMapG01L1_[i] ) mePnAmplMapG01L1_[i]->Reset();
      if ( mePnPedMapG01L1_[i] ) mePnPedMapG01L1_[i]->Reset();

      if ( mePnAmplMapG16L1_[i] ) mePnAmplMapG16L1_[i]->Reset();
      if ( mePnPedMapG16L1_[i] ) mePnPedMapG16L1_[i]->Reset();
    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {
      if ( mePnAmplMapG01L2_[i] ) mePnAmplMapG01L2_[i]->Reset();
      if ( mePnPedMapG01L2_[i] ) mePnPedMapG01L2_[i]->Reset();

      if ( mePnAmplMapG16L2_[i] ) mePnAmplMapG16L2_[i]->Reset();
      if ( mePnPedMapG16L2_[i] ) mePnPedMapG16L2_[i]->Reset();
    }
  }

}

void EELedTask::reset(void) {

}

void EELedTask::setup(void){

  init_ = true;

  char histo[200];

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask");

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT shape %s L1", Numbers::sEE(i+1).c_str());
        meShapeMapL1_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
        meShapeMapL1_[i]->setAxisTitle("channel", 1);
        meShapeMapL1_[i]->setAxisTitle("sample", 2);
        meShapeMapL1_[i]->setAxisTitle("amplitude", 3);
        dqmStore_->tag(meShapeMapL1_[i], i+1);
        sprintf(histo, "EELDT amplitude %s L1", Numbers::sEE(i+1).c_str());
        meAmplMapL1_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
        meAmplMapL1_[i]->setAxisTitle("jx", 1);
        meAmplMapL1_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meAmplMapL1_[i], i+1);
        sprintf(histo, "EELDT timing %s L1", Numbers::sEE(i+1).c_str());
        meTimeMapL1_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
        meTimeMapL1_[i]->setAxisTitle("jx", 1);
        meTimeMapL1_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meTimeMapL1_[i], i+1);
        sprintf(histo, "EELDT amplitude over PN %s L1", Numbers::sEE(i+1).c_str());
        meAmplPNMapL1_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
        meAmplPNMapL1_[i]->setAxisTitle("jx", 1);
        meAmplPNMapL1_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meAmplPNMapL1_[i], i+1);
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT shape %s L2", Numbers::sEE(i+1).c_str());
        meShapeMapL2_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
        meShapeMapL2_[i]->setAxisTitle("channel", 1);
        meShapeMapL2_[i]->setAxisTitle("sample", 2);
        meShapeMapL2_[i]->setAxisTitle("amplitude", 3);
        dqmStore_->tag(meShapeMapL2_[i], i+1);
        sprintf(histo, "EELDT amplitude %s L2", Numbers::sEE(i+1).c_str());
        meAmplMapL2_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
        meAmplMapL2_[i]->setAxisTitle("jx", 1);
        meAmplMapL2_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meAmplMapL2_[i], i+1);
        sprintf(histo, "EELDT timing %s L2", Numbers::sEE(i+1).c_str());
        meTimeMapL2_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
        meTimeMapL2_[i]->setAxisTitle("jx", 1);
        meTimeMapL2_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meTimeMapL2_[i], i+1);
        sprintf(histo, "EELDT amplitude over PN %s L2", Numbers::sEE(i+1).c_str());
        meAmplPNMapL2_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
        meAmplPNMapL2_[i]->setAxisTitle("jx", 1);
        meAmplPNMapL2_[i]->setAxisTitle("jy", 2);
        dqmStore_->tag(meAmplPNMapL2_[i], i+1);
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN");

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN/Gain01");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT PNs amplitude %s G01 L1", Numbers::sEE(i+1).c_str());
        mePnAmplMapG01L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnAmplMapG01L1_[i]->setAxisTitle("channel", 1);
        mePnAmplMapG01L1_[i]->setAxisTitle("amplitude", 2);
        dqmStore_->tag(mePnAmplMapG01L1_[i], i+1);
        sprintf(histo, "EELDT PNs pedestal %s G01 L1", Numbers::sEE(i+1).c_str());
        mePnPedMapG01L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnPedMapG01L1_[i]->setAxisTitle("channel", 1);
        mePnPedMapG01L1_[i]->setAxisTitle("pedestal", 2);
        dqmStore_->tag(mePnPedMapG01L1_[i], i+1);
      }

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN/Gain16");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT PNs amplitude %s G16 L1", Numbers::sEE(i+1).c_str());
        mePnAmplMapG16L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnAmplMapG16L1_[i]->setAxisTitle("channel", 1);
        mePnAmplMapG16L1_[i]->setAxisTitle("amplitude", 2);
        dqmStore_->tag(mePnAmplMapG16L1_[i], i+1);
        sprintf(histo, "EELDT PNs pedestal %s G16 L1", Numbers::sEE(i+1).c_str());
        mePnPedMapG16L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnPedMapG16L1_[i]->setAxisTitle("channel", 1);
        mePnPedMapG16L1_[i]->setAxisTitle("pedestal", 2);
        dqmStore_->tag(mePnPedMapG16L1_[i], i+1);
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN");

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN/Gain01");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT PNs amplitude %s G01 L2", Numbers::sEE(i+1).c_str());
        mePnAmplMapG01L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnAmplMapG01L2_[i]->setAxisTitle("amplitude", 2);
        mePnAmplMapG01L2_[i]->setAxisTitle("channel", 1);
        dqmStore_->tag(mePnAmplMapG01L2_[i], i+1);
        sprintf(histo, "EELDT PNs pedestal %s G01 L2", Numbers::sEE(i+1).c_str());
        mePnPedMapG01L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnPedMapG01L2_[i]->setAxisTitle("channel", 1);
        mePnPedMapG01L2_[i]->setAxisTitle("pedestal", 2);
        dqmStore_->tag(mePnPedMapG01L2_[i], i+1);
      }

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN/Gain16");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EELDT PNs amplitude %s G16 L2", Numbers::sEE(i+1).c_str());
        mePnAmplMapG16L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnAmplMapG16L2_[i]->setAxisTitle("channel", 1);
        mePnAmplMapG16L2_[i]->setAxisTitle("amplitude", 2);
        dqmStore_->tag(mePnAmplMapG16L2_[i], i+1);
        sprintf(histo, "EELDT PNs pedestal %s G16 L2", Numbers::sEE(i+1).c_str());
        mePnPedMapG16L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
        mePnPedMapG16L2_[i]->setAxisTitle("channel", 1);
        mePnPedMapG16L2_[i]->setAxisTitle("pedestal", 2);
        dqmStore_->tag(mePnPedMapG16L2_[i], i+1);
      }

    }

  }

}

void EELedTask::cleanup(void){

  if ( ! init_ ) return;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask");

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1");
      for (int i = 0; i < 18; i++) {
        if ( meShapeMapL1_[i] )  dqmStore_->removeElement( meShapeMapL1_[i]->getName() );
        meShapeMapL1_[i] = 0;
        if ( meAmplMapL1_[i] ) dqmStore_->removeElement( meAmplMapL1_[i]->getName() );
        meAmplMapL1_[i] = 0;
        if ( meTimeMapL1_[i] ) dqmStore_->removeElement( meTimeMapL1_[i]->getName() );
        meTimeMapL1_[i] = 0;
        if ( meAmplPNMapL1_[i] ) dqmStore_->removeElement( meAmplPNMapL1_[i]->getName() );
        meAmplPNMapL1_[i] = 0;
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2");
      for (int i = 0; i < 18; i++) {
        if ( meShapeMapL2_[i] )  dqmStore_->removeElement( meShapeMapL2_[i]->getName() );
        meShapeMapL2_[i] = 0;
        if ( meAmplMapL2_[i] ) dqmStore_->removeElement( meAmplMapL2_[i]->getName() );
        meAmplMapL2_[i] = 0;
        if ( meTimeMapL2_[i] ) dqmStore_->removeElement( meTimeMapL2_[i]->getName() );
        meTimeMapL2_[i] = 0;
        if ( meAmplPNMapL2_[i] ) dqmStore_->removeElement( meAmplPNMapL2_[i]->getName() );
        meAmplPNMapL2_[i] = 0;
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 1) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN");

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN/Gain01");
      for (int i = 0; i < 18; i++) {
        if ( mePnAmplMapG01L1_[i] ) dqmStore_->removeElement( mePnAmplMapG01L1_[i]->getName() );
        mePnAmplMapG01L1_[i] = 0;
        if ( mePnPedMapG01L1_[i] ) dqmStore_->removeElement( mePnPedMapG01L1_[i]->getName() );
        mePnPedMapG01L1_[i] = 0;
      }

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led1/PN/Gain16");
      for (int i = 0; i < 18; i++) {
        if ( mePnAmplMapG16L1_[i] ) dqmStore_->removeElement( mePnAmplMapG16L1_[i]->getName() );
        mePnAmplMapG16L1_[i] = 0;
        if ( mePnPedMapG16L1_[i] ) dqmStore_->removeElement( mePnPedMapG16L1_[i]->getName() );
        mePnPedMapG16L1_[i] = 0;
      }

    }

    if ( find(ledWavelengths_.begin(), ledWavelengths_.end(), 2) != ledWavelengths_.end() ) {

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN");

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN/Gain01");
      for (int i = 0; i < 18; i++) {
        if ( mePnAmplMapG01L2_[i] ) dqmStore_->removeElement( mePnAmplMapG01L2_[i]->getName() );
        mePnAmplMapG01L2_[i] = 0;
        if ( mePnPedMapG01L2_[i] ) dqmStore_->removeElement( mePnPedMapG01L2_[i]->getName() );
        mePnPedMapG01L2_[i] = 0;
      }

      dqmStore_->setCurrentFolder(prefixME_ + "/EELedTask/Led2/PN/Gain16");
      for (int i = 0; i < 18; i++) {
        if ( mePnAmplMapG16L2_[i] ) dqmStore_->removeElement( mePnAmplMapG16L2_[i]->getName() );
        mePnAmplMapG16L2_[i] = 0;
        if ( mePnPedMapG16L2_[i] ) dqmStore_->removeElement( mePnPedMapG16L2_[i]->getName() );
        mePnPedMapG16L2_[i] = 0;
      }

    }

  }

  init_ = false;

}

void EELedTask::endJob(void){

  LogInfo("EELedTask") << "analyzed " << ievt_ << " events";

  if ( enableCleanup_ ) this->cleanup();

}

void EELedTask::analyze(const Event& e, const EventSetup& c){

  bool enable = false;
  int runType[18];
  for (int i=0; i<18; i++) runType[i] = -1;
  int rtHalf[18];
  for (int i=0; i<18; i++) rtHalf[i] = -1;
  int waveLength[18];
  for (int i=0; i<18; i++) waveLength[i] = -1;

  Handle<EcalRawDataCollection> dcchs;

  if ( e.getByLabel(EcalRawDataCollection_, dcchs) ) {

    for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin(); dcchItr != dcchs->end(); ++dcchItr ) {

      if ( Numbers::subDet( *dcchItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *dcchItr, EcalEndcap );

      runType[ism-1] = dcchItr->getRunType();
      rtHalf[ism-1] = dcchItr->getRtHalf();
      waveLength[ism-1] = dcchItr->getEventSettings().wavelength;

      if ( dcchItr->getRunType() == EcalDCCHeaderBlock::LED_STD ||
           dcchItr->getRunType() == EcalDCCHeaderBlock::LED_GAP ) enable = true;

    }

  } else {

    LogWarning("EELedTask") << EcalRawDataCollection_ << " not available";

  }

  if ( ! enable ) return;

  if ( ! init_ ) this->setup();

  ievt_++;

  Handle<EEDigiCollection> digis;

  if ( e.getByLabel(EEDigiCollection_, digis) ) {

    int need = digis->size();
    LogDebug("EELedTask") << "event " << ievt_ << " digi collection size " << need;

    for ( EEDigiCollection::const_iterator digiItr = digis->begin(); digiItr != digis->end(); ++digiItr ) {

      EEDetId id = digiItr->id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ! ( runType[ism-1] == EcalDCCHeaderBlock::LED_STD ||
               runType[ism-1] == EcalDCCHeaderBlock::LED_GAP ) ) continue;

      if ( runType[ism-1] == EcalDCCHeaderBlock::LED_GAP &&
           rtHalf[ism-1] != Numbers::RtHalf(id) ) continue;

      int ic = Numbers::icEE(ism, ix, iy);

      EEDataFrame dataframe = (*digiItr);

      for (int i = 0; i < 10; i++) {

        int adc = dataframe.sample(i).adc();
        float gain = 1.;

        MonitorElement* meShapeMap = 0;

        if ( dataframe.sample(i).gainId() == 1 ) gain = 1./12.;
        if ( dataframe.sample(i).gainId() == 2 ) gain = 1./ 6.;
        if ( dataframe.sample(i).gainId() == 3 ) gain = 1./ 1.;

        if ( Numbers::RtHalf(id) == 0 || Numbers::RtHalf(id) == 1 ) {

          if ( waveLength[ism-1] == 0 ) meShapeMap = meShapeMapL1_[ism-1];
          if ( waveLength[ism-1] == 2 ) meShapeMap = meShapeMapL2_[ism-1];

        } else {

          LogWarning("EELedTask") << " RtHalf = " << Numbers::RtHalf(id);

        }

//        float xval = float(adc) * gain;
        float xval = float(adc);

        if ( meShapeMap ) meShapeMap->Fill(ic - 0.5, i + 0.5, xval);

      }

    }

  } else {

    LogWarning("EELedTask") << EEDigiCollection_ << " not available";

  }

  float adcPN[80];
  for ( int i = 0; i < 80; i++ ) {
    adcPN[i] = 0.;
  }

  Handle<EcalPnDiodeDigiCollection> pns;

  if ( e.getByLabel(EcalPnDiodeDigiCollection_, pns) ) {

    int nep = pns->size();
    LogDebug("EELedTask") << "event " << ievt_ << " pns collection size " << nep;

    for ( EcalPnDiodeDigiCollection::const_iterator pnItr = pns->begin(); pnItr != pns->end(); ++pnItr ) {

      if ( Numbers::subDet( pnItr->id() ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( pnItr->id() );

      int num = pnItr->id().iPnId();

      if ( ! ( runType[ism-1] == EcalDCCHeaderBlock::LED_STD ||
               runType[ism-1] == EcalDCCHeaderBlock::LED_GAP ) ) continue;

      float xvalped = 0.;

      for (int i = 0; i < 4; i++) {

        int adc = pnItr->sample(i).adc();

        MonitorElement* mePNPed = 0;

        if ( pnItr->sample(i).gainId() == 0 ) {
          if ( waveLength[ism-1] == 0 ) mePNPed = mePnPedMapG01L1_[ism-1];
          if ( waveLength[ism-1] == 2 ) mePNPed = mePnPedMapG01L2_[ism-1];
        }
        if ( pnItr->sample(i).gainId() == 1 ) {
          if ( waveLength[ism-1] == 0 ) mePNPed = mePnPedMapG16L1_[ism-1];
          if ( waveLength[ism-1] == 2 ) mePNPed = mePnPedMapG16L2_[ism-1];
        }

        float xval = float(adc);

        if ( mePNPed ) mePNPed->Fill(num - 0.5, xval);

        xvalped = xvalped + xval;

      }

      xvalped = xvalped / 4;

      float xvalmax = 0.;

      MonitorElement* mePN = 0;

      for (int i = 0; i < 50; i++) {

        int adc = pnItr->sample(i).adc();

        float xval = float(adc);

        if ( xval >= xvalmax ) xvalmax = xval;

      }

      xvalmax = xvalmax - xvalped;

      if ( pnItr->sample(0).gainId() == 0 ) {
        if ( waveLength[ism-1] == 0 ) mePN = mePnAmplMapG01L1_[ism-1];
        if ( waveLength[ism-1] == 2 ) mePN = mePnAmplMapG01L2_[ism-1];
      }
      if ( pnItr->sample(0).gainId() == 1 ) {
        if ( waveLength[ism-1] == 0 ) mePN = mePnAmplMapG16L1_[ism-1];
        if ( waveLength[ism-1] == 2 ) mePN = mePnAmplMapG16L2_[ism-1];
      }

      if ( mePN ) mePN->Fill(num - 0.5, xvalmax);

      int ipn = NumbersPn::ipnEE( ism, num );

      adcPN[ipn] = xvalmax;

    }

  } else {

    LogWarning("EELedTask") << EcalPnDiodeDigiCollection_ << " not available";

  }

  Handle<EcalUncalibratedRecHitCollection> hits;

  if ( e.getByLabel(EcalUncalibratedRecHitCollection_, hits) ) {

    int neh = hits->size();
    LogDebug("EELedTask") << "event " << ievt_ << " hits collection size " << neh;

    for ( EcalUncalibratedRecHitCollection::const_iterator hitItr = hits->begin(); hitItr != hits->end(); ++hitItr ) {

      EEDetId id = hitItr->id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( ! ( runType[ism-1] == EcalDCCHeaderBlock::LED_STD ||
               runType[ism-1] == EcalDCCHeaderBlock::LED_GAP ) ) continue;

      if ( runType[ism-1] == EcalDCCHeaderBlock::LED_GAP &&
           rtHalf[ism-1] != Numbers::RtHalf(id) ) continue;

      MonitorElement* meAmplMap = 0;
      MonitorElement* meTimeMap = 0;
      MonitorElement* meAmplPNMap = 0;

      if ( Numbers::RtHalf(id) == 0 || Numbers::RtHalf(id) == 1 ) {

        if ( waveLength[ism-1] == 0 ) {
          meAmplMap = meAmplMapL1_[ism-1];
          meTimeMap = meTimeMapL1_[ism-1];
          meAmplPNMap = meAmplPNMapL1_[ism-1];
        }
        if ( waveLength[ism-1] == 2 ) {
          meAmplMap = meAmplMapL2_[ism-1];
          meTimeMap = meTimeMapL2_[ism-1];
          meAmplPNMap = meAmplPNMapL2_[ism-1];
        }

      } else {

        LogWarning("EELedTask") << " RtHalf = " << Numbers::RtHalf(id);

      }

      float xval = hitItr->amplitude();
      if ( xval <= 0. ) xval = 0.0;
      float yval = hitItr->jitter() + 6.0;
      if ( yval <= 0. ) yval = 0.0;
      float zval = hitItr->pedestal();
      if ( zval <= 0. ) zval = 0.0;

      if ( meAmplMap ) meAmplMap->Fill(xix, xiy, xval);

      if ( xval > 16. ) {
        if ( meTimeMap ) meTimeMap->Fill(xix, xiy, yval);
      }

      float wval = 0.;

      std::vector<int> PNsInLM = NumbersPn::getPNs( ism, ix, iy );
      int refPn = PNsInLM[0];

      if ( adcPN[refPn] != 0. ) wval = xval /  adcPN[refPn];

      if ( meAmplPNMap ) meAmplPNMap->Fill(xix, xiy, wval);

    }

  } else {

    LogWarning("EELedTask") << EcalUncalibratedRecHitCollection_ << " not available";

  }

}

