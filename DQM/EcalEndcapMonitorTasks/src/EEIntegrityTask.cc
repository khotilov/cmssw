/*
 * \file EEIntegrityTask.cc
 *
 * $Date: 2009/06/23 06:41:51 $
 * $Revision: 1.44 $
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

#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/EcalDetId/interface/EcalDetIdCollections.h"

#include <DQM/EcalCommon/interface/Numbers.h>

#include <DQM/EcalEndcapMonitorTasks/interface/EEIntegrityTask.h>

using namespace cms;
using namespace edm;
using namespace std;

EEIntegrityTask::EEIntegrityTask(const ParameterSet& ps){

  init_ = false;

  dqmStore_ = Service<DQMStore>().operator->();

  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  mergeRuns_ = ps.getUntrackedParameter<bool>("mergeRuns", false);

  EEDetIdCollection0_ =  ps.getParameter<edm::InputTag>("EEDetIdCollection0");
  EEDetIdCollection1_ =  ps.getParameter<edm::InputTag>("EEDetIdCollection1");
  EEDetIdCollection2_ =  ps.getParameter<edm::InputTag>("EEDetIdCollection2");
  EEDetIdCollection3_ =  ps.getParameter<edm::InputTag>("EEDetIdCollection3");
  EcalElectronicsIdCollection1_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection1");
  EcalElectronicsIdCollection2_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection2");
  EcalElectronicsIdCollection3_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection3");
  EcalElectronicsIdCollection4_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection4");
  EcalElectronicsIdCollection5_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection5");
  EcalElectronicsIdCollection6_ = ps.getParameter<edm::InputTag>("EcalElectronicsIdCollection6");

  meIntegrityDCCSize = 0;
  for (int i = 0; i < 18; i++) {
    meIntegrityGain[i] = 0;
    meIntegrityChId[i] = 0;
    meIntegrityGainSwitch[i] = 0;
    meIntegrityTTId[i] = 0;
    meIntegrityTTBlockSize[i] = 0;
    meIntegrityMemChId[i] = 0;
    meIntegrityMemGain[i] = 0;
    meIntegrityMemTTId[i] = 0;
    meIntegrityMemTTBlockSize[i] = 0;
  }

}


EEIntegrityTask::~EEIntegrityTask(){

}

void EEIntegrityTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask");
    dqmStore_->rmdir(prefixME_ + "/EEIntegrityTask");
  }

  Numbers::initGeometry(c, false);

}

void EEIntegrityTask::beginRun(const Run& r, const EventSetup& c) {

  if ( ! mergeRuns_ ) this->reset();

}

void EEIntegrityTask::endRun(const Run& r, const EventSetup& c) {

}

void EEIntegrityTask::reset(void) {

  if ( meIntegrityDCCSize ) meIntegrityDCCSize->Reset();
  for (int i = 0; i < 18; i++) {
    if ( meIntegrityGain[i] ) meIntegrityGain[i]->Reset();
    if ( meIntegrityChId[i] ) meIntegrityChId[i]->Reset();
    if ( meIntegrityGainSwitch[i] ) meIntegrityGainSwitch[i]->Reset();
    if ( meIntegrityTTId[i] ) meIntegrityTTId[i]->Reset();
    if ( meIntegrityTTBlockSize[i] ) meIntegrityTTBlockSize[i]->Reset();
    if ( meIntegrityMemChId[i] ) meIntegrityMemChId[i]->Reset();
    if ( meIntegrityMemGain[i] ) meIntegrityMemGain[i]->Reset();
    if ( meIntegrityMemTTId[i] ) meIntegrityMemTTId[i]->Reset();
    if ( meIntegrityMemTTBlockSize[i] ) meIntegrityMemTTBlockSize[i]->Reset();
  }

}

void EEIntegrityTask::setup(void){

  init_ = true;

  char histo[200];

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask");

    // checking when number of towers in data different than expected from header
    sprintf(histo, "EEIT DCC size error");
    meIntegrityDCCSize = dqmStore_->book1D(histo, histo, 18, 1, 19.);
    for (int i = 0; i < 18; i++) {
      meIntegrityDCCSize->setBinLabel(i+1, Numbers::sEE(i+1).c_str(), 1);
    }

    // checking when the gain is 0
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/Gain");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT gain %s", Numbers::sEE(i+1).c_str());
      meIntegrityGain[i] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
      meIntegrityGain[i]->setAxisTitle("jx", 1);
      meIntegrityGain[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meIntegrityGain[i], i+1);
    }

    // checking when channel has unexpected or invalid ID
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/ChId");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT ChId %s", Numbers::sEE(i+1).c_str());
      meIntegrityChId[i] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
      meIntegrityChId[i]->setAxisTitle("jx", 1);
      meIntegrityChId[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meIntegrityChId[i], i+1);
    }

    // checking when channel has unexpected or invalid ID
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/GainSwitch");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT gain switch %s", Numbers::sEE(i+1).c_str());
      meIntegrityGainSwitch[i] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
      meIntegrityGainSwitch[i]->setAxisTitle("jx", 1);
      meIntegrityGainSwitch[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meIntegrityGainSwitch[i], i+1);
    }

    // checking when trigger tower has unexpected or invalid ID
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/TTId");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT TTId %s", Numbers::sEE(i+1).c_str());
      meIntegrityTTId[i] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
      meIntegrityTTId[i]->setAxisTitle("jx", 1);
      meIntegrityTTId[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meIntegrityTTId[i], i+1);
    }

    // checking when trigger tower has unexpected or invalid size
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/TTBlockSize");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT TTBlockSize %s", Numbers::sEE(i+1).c_str());
      meIntegrityTTBlockSize[i] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
      meIntegrityTTBlockSize[i]->setAxisTitle("jx", 1);
      meIntegrityTTBlockSize[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meIntegrityTTBlockSize[i], i+1);
    }

    // checking when mem channels have unexpected ID
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemChId");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT MemChId %s", Numbers::sEE(i+1).c_str());
      meIntegrityMemChId[i] = dqmStore_->book2D(histo, histo, 10, 0., 10., 5, 0., 5.);
      meIntegrityMemChId[i]->setAxisTitle("pseudo-strip", 1);
      meIntegrityMemChId[i]->setAxisTitle("channel", 2);
      dqmStore_->tag(meIntegrityMemChId[i], i+1);
    }

    // checking when mem samples have second bit encoding the gain different from 0
    // note: strictly speaking, this does not corrupt the mem sample gain value (since only first bit is considered)
    // but indicates that data are not completely correct
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemGain");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT MemGain %s", Numbers::sEE(i+1).c_str());
      meIntegrityMemGain[i] = dqmStore_->book2D(histo, histo, 10, 0., 10., 5, 0., 5.);
      meIntegrityMemGain[i]->setAxisTitle("pseudo-strip", 1);
      meIntegrityMemGain[i]->setAxisTitle("channel", 2);
      dqmStore_->tag(meIntegrityMemGain[i], i+1);
    }

    // checking when mem tower block has unexpected ID
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemTTId");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT MemTTId %s", Numbers::sEE(i+1).c_str());
      meIntegrityMemTTId[i] = dqmStore_->book2D(histo, histo, 2, 0., 2., 1, 0., 1.);
      meIntegrityMemTTId[i]->setAxisTitle("pseudo-strip", 1);
      meIntegrityMemTTId[i]->setAxisTitle("channel", 2);
      dqmStore_->tag(meIntegrityMemTTId[i], i+1);
    }

    // checking when mem tower block has invalid size
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemSize");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEIT MemSize %s", Numbers::sEE(i+1).c_str());
      meIntegrityMemTTBlockSize[i] = dqmStore_->book2D(histo, histo, 2, 0., 2., 1, 0., 1.);
      meIntegrityMemTTBlockSize[i]->setAxisTitle("pseudo-strip", 1);
      meIntegrityMemTTBlockSize[i]->setAxisTitle("channel", 2);
      dqmStore_->tag(meIntegrityMemTTBlockSize[i], i+1);
    }

  }

}

void EEIntegrityTask::cleanup(void){

  if ( ! init_ ) return;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask");

    if ( meIntegrityDCCSize ) dqmStore_->removeElement( meIntegrityDCCSize->getName() );
    meIntegrityDCCSize = 0;

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/Gain");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityGain[i] ) dqmStore_->removeElement( meIntegrityGain[i]->getName() );
      meIntegrityGain[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/ChId");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityChId[i] ) dqmStore_->removeElement( meIntegrityChId[i]->getName() );
      meIntegrityChId[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/GainSwitch");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityGainSwitch[i] ) dqmStore_->removeElement( meIntegrityGainSwitch[i]->getName() );
      meIntegrityGainSwitch[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/TTId");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityTTId[i] ) dqmStore_->removeElement( meIntegrityTTId[i]->getName() );
      meIntegrityTTId[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/TTBlockSize");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityTTBlockSize[i] ) dqmStore_->removeElement( meIntegrityTTBlockSize[i]->getName() );
      meIntegrityTTBlockSize[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemChId");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityMemChId[i] ) dqmStore_->removeElement( meIntegrityMemChId[i]->getName() );
      meIntegrityMemChId[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemGain");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityMemGain[i] ) dqmStore_->removeElement( meIntegrityMemGain[i]->getName() );
      meIntegrityMemGain[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemTTId");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityMemTTId[i] ) dqmStore_->removeElement( meIntegrityMemTTId[i]->getName() );
      meIntegrityMemTTId[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EEIntegrityTask/MemSize");
    for (int i = 0; i < 18; i++) {
      if ( meIntegrityMemTTBlockSize[i] ) dqmStore_->removeElement( meIntegrityMemTTBlockSize[i]->getName() );
      meIntegrityMemTTBlockSize[i] = 0;
    }

  }

  init_ = false;

}

void EEIntegrityTask::endJob(void){

  LogInfo("EEIntegrityTask") << "analyzed " << ievt_ << " events";

  if ( enableCleanup_ ) this->cleanup();

}

void EEIntegrityTask::analyze(const Event& e, const EventSetup& c){

  if ( ! enable ) return;

  if ( ! init_ ) this->setup();

  ievt_++;

  Handle<EEDetIdCollection> ids0;

  if ( e.getByLabel(EEDetIdCollection0_, ids0) ) {

    for ( EEDetIdCollection::const_iterator idItr = ids0->begin(); idItr != ids0->end(); ++idItr ) {

      int ism = Numbers::iSM( *idItr );

      float xism = ism - 0.5;

      if ( meIntegrityDCCSize ) meIntegrityDCCSize->Fill(xism);

    }

  } else {

//    LogWarning("EEIntegrityTask") << EEDetIdCollection0_ << " not available";

  }

  Handle<EEDetIdCollection> ids1;

  if ( e.getByLabel(EEDetIdCollection1_, ids1) ) {

    for ( EEDetIdCollection::const_iterator idItr = ids1->begin(); idItr != ids1->end(); ++idItr ) {

      EEDetId id = (*idItr);

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( meIntegrityGain[ism-1] ) meIntegrityGain[ism-1]->Fill(xix, xiy);

    }

  } else {

    LogWarning("EEIntegrityTask") << EEDetIdCollection1_ << " not available";

  }

  Handle<EEDetIdCollection> ids2;

  if ( e.getByLabel(EEDetIdCollection2_, ids2) ) {

    for ( EEDetIdCollection::const_iterator idItr = ids2->begin(); idItr != ids2->end(); ++idItr ) {

      EEDetId id = (*idItr);

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( meIntegrityChId[ism-1] ) meIntegrityChId[ism-1]->Fill(xix, xiy);

    }

  } else {

    LogWarning("EEIntegrityTask") << EEDetIdCollection2_ << " not available";

  }

  Handle<EEDetIdCollection> ids3;

  if ( e.getByLabel(EEDetIdCollection3_, ids3) ) {

    for ( EEDetIdCollection::const_iterator idItr = ids3->begin(); idItr != ids3->end(); ++idItr ) {

      EEDetId id = (*idItr);

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( meIntegrityGainSwitch[ism-1] ) meIntegrityGainSwitch[ism-1]->Fill(xix, xiy);

    }

  } else {

    LogWarning("EEIntegrityTask") << EEDetIdCollection3_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids4;

  if ( e.getByLabel(EcalElectronicsIdCollection1_, ids4) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids4->begin(); idItr != ids4->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      vector<DetId> crystals = Numbers::crystals( *idItr );

      for ( unsigned int i=0; i<crystals.size(); i++ ) {

      EEDetId id = crystals[i];

      int ix = id.ix();
      int iy = id.iy();

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( meIntegrityTTId[ism-1] ) meIntegrityTTId[ism-1]->Fill(xix, xiy);

      }

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection1_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids5;

  if ( e.getByLabel(EcalElectronicsIdCollection2_, ids5) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids5->begin(); idItr != ids5->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      vector<DetId> crystals = Numbers::crystals( *idItr );

      for ( unsigned int i=0; i<crystals.size(); i++ ) {

      EEDetId id = crystals[i];

      int ix = id.ix();
      int iy = id.iy();

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      if ( meIntegrityTTBlockSize[ism-1] ) meIntegrityTTBlockSize[ism-1]->Fill(xix, xiy);

      }

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection2_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids6;

  if ( e.getByLabel(EcalElectronicsIdCollection3_, ids6) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids6->begin(); idItr != ids6->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      int itt   = idItr->towerId();
      float iTt = itt + 0.5 - 69;

      if ( meIntegrityMemTTId[ism-1] ) meIntegrityMemTTId[ism-1]->Fill(iTt,0);

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection3_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids7;

  if ( e.getByLabel(EcalElectronicsIdCollection4_, ids7) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids7->begin(); idItr != ids7->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      int itt   = idItr->towerId();
      float iTt = itt + 0.5 - 69;

      if ( meIntegrityMemTTBlockSize[ism-1] ) meIntegrityMemTTBlockSize[ism-1]->Fill(iTt,0);

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection4_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids8;

  if (  e.getByLabel(EcalElectronicsIdCollection5_, ids8) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids8->begin(); idItr != ids8->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      int chid = idItr->channelId();
      int ie = EEIntegrityTask::chMemAbscissa[chid-1];
      int ip = EEIntegrityTask::chMemOrdinate[chid-1];

      int itt = idItr->towerId();
      ie += (itt-69)*5;

      float xix = ie - 0.5;
      float xiy = ip - 0.5;

      if ( meIntegrityMemChId[ism-1] ) meIntegrityMemChId[ism-1]->Fill(xix,xiy);

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection5_ << " not available";

  }

  Handle<EcalElectronicsIdCollection> ids9;

  if ( e.getByLabel(EcalElectronicsIdCollection6_, ids9) ) {

    for ( EcalElectronicsIdCollection::const_iterator idItr = ids9->begin(); idItr != ids9->end(); ++idItr ) {

      if ( Numbers::subDet( *idItr ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( *idItr );

      int chid = idItr->channelId();
      int ie = EEIntegrityTask::chMemAbscissa[chid-1];
      int ip = EEIntegrityTask::chMemOrdinate[chid-1];

      int itt = idItr->towerId();
      ie += (itt-69)*5;

      float xix = ie - 0.5;
      float xiy = ip - 0.5;

      if ( meIntegrityMemGain[ism-1] ) meIntegrityMemGain[ism-1]->Fill(xix,xiy);

    }

  } else {

    LogWarning("EEIntegrityTask") << EcalElectronicsIdCollection6_ << " not available";

  }

}//  end analyze

const int  EEIntegrityTask::chMemAbscissa [25] = {
    1, 1, 1, 1, 1,
    2, 2, 2, 2, 2,
    3, 3, 3, 3, 3,
    4, 4, 4, 4, 4,
    5, 5, 5, 5, 5
};

const int  EEIntegrityTask::chMemOrdinate [25] = {
    1, 2, 3, 4, 5,
    5, 4, 3, 2, 1,
    1, 2, 3, 4, 5,
    5, 4, 3, 2, 1,
    1, 2, 3, 4, 5
};

