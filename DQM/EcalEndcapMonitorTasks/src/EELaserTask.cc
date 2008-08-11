/*
 * \file EELaserTask.cc
 *
 * $Date: 2008/08/11 17:47:13 $
 * $Revision: 1.48 $
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

#include <DQM/EcalEndcapMonitorTasks/interface/EELaserTask.h>

using namespace cms;
using namespace edm;
using namespace std;

EELaserTask::EELaserTask(const ParameterSet& ps){

  init_ = false;

  dqmStore_ = Service<DQMStore>().operator->();

  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  mergeRuns_ = ps.getUntrackedParameter<bool>("mergeRuns", false);

  EcalRawDataCollection_ = ps.getParameter<edm::InputTag>("EcalRawDataCollection");
  EEDigiCollection_ = ps.getParameter<edm::InputTag>("EEDigiCollection");
  EcalPnDiodeDigiCollection_ = ps.getParameter<edm::InputTag>("EcalPnDiodeDigiCollection");
  EcalUncalibratedRecHitCollection_ = ps.getParameter<edm::InputTag>("EcalUncalibratedRecHitCollection");

  for (int i = 0; i < 18; i++) {
    meShapeMapL1A_[i] = 0;
    meAmplMapL1A_[i] = 0;
    meTimeMapL1A_[i] = 0;
    meAmplPNMapL1A_[i] = 0;
    meShapeMapL1B_[i] = 0;
    meAmplMapL1B_[i] = 0;
    meTimeMapL1B_[i] = 0;
    meAmplPNMapL1B_[i] = 0;
    mePnAmplMapG01L1_[i] = 0;
    mePnPedMapG01L1_[i] = 0;
    mePnAmplMapG16L1_[i] = 0;
    mePnPedMapG16L1_[i] = 0;

    meShapeMapL2A_[i] = 0;
    meAmplMapL2A_[i] = 0;
    meTimeMapL2A_[i] = 0;
    meAmplPNMapL2A_[i] = 0;
    meShapeMapL2B_[i] = 0;
    meAmplMapL2B_[i] = 0;
    meTimeMapL2B_[i] = 0;
    meAmplPNMapL2B_[i] = 0;
    mePnAmplMapG01L2_[i] = 0;
    mePnPedMapG01L2_[i] = 0;
    mePnAmplMapG16L2_[i] = 0;
    mePnPedMapG16L2_[i] = 0;

    meShapeMapL3A_[i] = 0;
    meAmplMapL3A_[i] = 0;
    meTimeMapL3A_[i] = 0;
    meAmplPNMapL3A_[i] = 0;
    meShapeMapL3B_[i] = 0;
    meAmplMapL3B_[i] = 0;
    meTimeMapL3B_[i] = 0;
    meAmplPNMapL3B_[i] = 0;
    mePnAmplMapG01L3_[i] = 0;
    mePnPedMapG01L3_[i] = 0;
    mePnAmplMapG16L3_[i] = 0;
    mePnPedMapG16L3_[i] = 0;

    meShapeMapL4A_[i] = 0;
    meAmplMapL4A_[i] = 0;
    meTimeMapL4A_[i] = 0;
    meAmplPNMapL4A_[i] = 0;
    meShapeMapL4B_[i] = 0;
    meAmplMapL4B_[i] = 0;
    meTimeMapL4B_[i] = 0;
    meAmplPNMapL4B_[i] = 0;
    mePnAmplMapG01L4_[i] = 0;
    mePnPedMapG01L4_[i] = 0;
    mePnAmplMapG16L4_[i] = 0;
    mePnPedMapG16L4_[i] = 0;
  }

}

EELaserTask::~EELaserTask(){

}

void EELaserTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask");
    dqmStore_->rmdir(prefixME_ + "/EELaserTask");
  }

  Numbers::initGeometry(c, false);

}

void EELaserTask::beginRun(const Run& r, const EventSetup& c) {

  if ( ! mergeRuns_ ) this->reset();

}

void EELaserTask::endRun(const Run& r, const EventSetup& c) {

}

void EELaserTask::reset(void) {

 for (int i = 0; i < 18; i++) {
    if ( meShapeMapL1A_[i] )  meShapeMapL1A_[i]->Reset();
    if ( meAmplMapL1A_[i] ) meAmplMapL1A_[i]->Reset();
    if ( meTimeMapL1A_[i] ) meTimeMapL1A_[i]->Reset();
    if ( meAmplPNMapL1A_[i] ) meAmplPNMapL1A_[i]->Reset();

    if ( meShapeMapL1B_[i] )  meShapeMapL1B_[i]->Reset();
    if ( meAmplMapL1B_[i] ) meAmplMapL1B_[i]->Reset();
    if ( meTimeMapL1B_[i] ) meTimeMapL1B_[i]->Reset();
    if ( meAmplPNMapL1B_[i] ) meAmplPNMapL1B_[i]->Reset();

    if ( meShapeMapL2A_[i] )  meShapeMapL2A_[i]->Reset();
    if ( meAmplMapL2A_[i] ) meAmplMapL2A_[i]->Reset();
    if ( meTimeMapL2A_[i] ) meTimeMapL2A_[i]->Reset();
    if ( meAmplPNMapL2A_[i] ) meAmplPNMapL2A_[i]->Reset();

    if ( meShapeMapL2B_[i] )  meShapeMapL2B_[i]->Reset();
    if ( meAmplMapL2B_[i] ) meAmplMapL2B_[i]->Reset();
    if ( meTimeMapL2B_[i] ) meTimeMapL2B_[i]->Reset();
    if ( meAmplPNMapL2B_[i] ) meAmplPNMapL2B_[i]->Reset();

    if ( meShapeMapL3A_[i] )  meShapeMapL3A_[i]->Reset();
    if ( meAmplMapL3A_[i] ) meAmplMapL3A_[i]->Reset();
    if ( meTimeMapL3A_[i] ) meTimeMapL3A_[i]->Reset();
    if ( meAmplPNMapL3A_[i] ) meAmplPNMapL3A_[i]->Reset();

    if ( meShapeMapL3B_[i] )  meShapeMapL3B_[i]->Reset();
    if ( meAmplMapL3B_[i] ) meAmplMapL3B_[i]->Reset();
    if ( meTimeMapL3B_[i] ) meTimeMapL3B_[i]->Reset();
    if ( meAmplPNMapL3B_[i] ) meAmplPNMapL3B_[i]->Reset();

    if ( meShapeMapL4A_[i] )  meShapeMapL4A_[i]->Reset();
    if ( meAmplMapL4A_[i] ) meAmplMapL4A_[i]->Reset();
    if ( meTimeMapL4A_[i] ) meTimeMapL4A_[i]->Reset();
    if ( meAmplPNMapL4A_[i] ) meAmplPNMapL4A_[i]->Reset();

    if ( meShapeMapL4B_[i] )  meShapeMapL4B_[i]->Reset();
    if ( meAmplMapL4B_[i] ) meAmplMapL4B_[i]->Reset();
    if ( meTimeMapL4B_[i] ) meTimeMapL4B_[i]->Reset();
    if ( meAmplPNMapL4B_[i] ) meAmplPNMapL4B_[i]->Reset();

    if ( mePnAmplMapG01L1_[i] ) mePnAmplMapG01L1_[i]->Reset();
    if ( mePnPedMapG01L1_[i] ) mePnPedMapG01L1_[i]->Reset();

    if ( mePnAmplMapG16L1_[i] ) mePnAmplMapG16L1_[i]->Reset();
    if ( mePnPedMapG16L1_[i] ) mePnPedMapG16L1_[i]->Reset();

    if ( mePnAmplMapG01L2_[i] ) mePnAmplMapG01L2_[i]->Reset();
    if ( mePnPedMapG01L2_[i] ) mePnPedMapG01L2_[i]->Reset();

    if ( mePnAmplMapG16L2_[i] ) mePnAmplMapG16L2_[i]->Reset();
    if ( mePnPedMapG16L2_[i] ) mePnPedMapG16L2_[i]->Reset();

    if ( mePnAmplMapG01L3_[i] ) mePnAmplMapG01L3_[i]->Reset();
    if ( mePnPedMapG01L3_[i] ) mePnPedMapG01L3_[i]->Reset();

    if ( mePnAmplMapG16L3_[i] ) mePnAmplMapG16L3_[i]->Reset();
    if ( mePnPedMapG16L3_[i] ) mePnPedMapG16L3_[i]->Reset();

    if ( mePnAmplMapG01L4_[i] ) mePnAmplMapG01L4_[i]->Reset();
    if ( mePnPedMapG01L4_[i] ) mePnPedMapG01L4_[i]->Reset();

    if ( mePnAmplMapG16L4_[i] ) mePnAmplMapG16L4_[i]->Reset();
    if ( mePnPedMapG16L4_[i] ) mePnPedMapG16L4_[i]->Reset();
  }

}

void EELaserTask::setup(void){

  init_ = true;

  char histo[200];

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EELT shape %s L1A", Numbers::sEE(i+1).c_str());
      meShapeMapL1A_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL1A_[i]->setAxisTitle("channel", 1);
      meShapeMapL1A_[i]->setAxisTitle("sample", 2);
      meShapeMapL1A_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL1A_[i], i+1);
      sprintf(histo, "EELT amplitude %s L1A", Numbers::sEE(i+1).c_str());
      meAmplMapL1A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL1A_[i]->setAxisTitle("jx", 1);
      meAmplMapL1A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL1A_[i], i+1);
      sprintf(histo, "EELT timing %s L1A", Numbers::sEE(i+1).c_str());
      meTimeMapL1A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL1A_[i]->setAxisTitle("jx", 1);
      meTimeMapL1A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL1A_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L1A", Numbers::sEE(i+1).c_str());
      meAmplPNMapL1A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL1A_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL1A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL1A_[i], i+1);

      sprintf(histo, "EELT shape %s L1B", Numbers::sEE(i+1).c_str());
      meShapeMapL1B_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL1B_[i]->setAxisTitle("channel", 1);
      meShapeMapL1B_[i]->setAxisTitle("sample", 2);
      meShapeMapL1B_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL1B_[i], i+1);
      sprintf(histo, "EELT amplitude %s L1B", Numbers::sEE(i+1).c_str());
      meAmplMapL1B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL1B_[i]->setAxisTitle("jx", 1);
      meAmplMapL1B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL1B_[i], i+1);
      sprintf(histo, "EELT timing %s L1B", Numbers::sEE(i+1).c_str());
      meTimeMapL1B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL1B_[i]->setAxisTitle("jx", 1);
      meTimeMapL1B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL1B_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L1B", Numbers::sEE(i+1).c_str());
      meAmplPNMapL1B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL1B_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL1B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL1B_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EELT shape %s L2A", Numbers::sEE(i+1).c_str());
      meShapeMapL2A_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL2A_[i]->setAxisTitle("channel", 1);
      meShapeMapL2A_[i]->setAxisTitle("sample", 2);
      meShapeMapL2A_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL2A_[i], i+1);
      sprintf(histo, "EELT amplitude %s L2A", Numbers::sEE(i+1).c_str());
      meAmplMapL2A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL2A_[i]->setAxisTitle("jx", 1);
      meAmplMapL2A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL2A_[i], i+1);
      sprintf(histo, "EELT timing %s L2A", Numbers::sEE(i+1).c_str());
      meTimeMapL2A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL2A_[i]->setAxisTitle("jx", 1);
      meTimeMapL2A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL2A_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L2A", Numbers::sEE(i+1).c_str());
      meAmplPNMapL2A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL2A_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL2A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL2A_[i], i+1);

      sprintf(histo, "EELT shape %s L2B", Numbers::sEE(i+1).c_str());
      meShapeMapL2B_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL2B_[i]->setAxisTitle("channel", 1);
      meShapeMapL2B_[i]->setAxisTitle("sample", 2);
      meShapeMapL2B_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL2B_[i], i+1);
      sprintf(histo, "EELT amplitude %s L2B", Numbers::sEE(i+1).c_str());
      meAmplMapL2B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL2B_[i]->setAxisTitle("jx", 1);
      meAmplMapL2B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL2B_[i], i+1);
      sprintf(histo, "EELT timing %s L2B", Numbers::sEE(i+1).c_str());
      meTimeMapL2B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL2B_[i]->setAxisTitle("jx", 1);
      meTimeMapL2B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL2B_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L2B", Numbers::sEE(i+1).c_str());
      meAmplPNMapL2B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL2B_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL2B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL2B_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EELT shape %s L3A", Numbers::sEE(i+1).c_str());
      meShapeMapL3A_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL3A_[i]->setAxisTitle("channel", 1);
      meShapeMapL3A_[i]->setAxisTitle("sample", 2);
      meShapeMapL3A_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL3A_[i], i+1);
      sprintf(histo, "EELT amplitude %s L3A", Numbers::sEE(i+1).c_str());
      meAmplMapL3A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL3A_[i]->setAxisTitle("jx", 1);
      meAmplMapL3A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL3A_[i], i+1);
      sprintf(histo, "EELT timing %s L3A", Numbers::sEE(i+1).c_str());
      meTimeMapL3A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL3A_[i]->setAxisTitle("jx", 1);
      meTimeMapL3A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL3A_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L3A", Numbers::sEE(i+1).c_str());
      meAmplPNMapL3A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL3A_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL3A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL3A_[i], i+1);

      sprintf(histo, "EELT shape %s L3B", Numbers::sEE(i+1).c_str());
      meShapeMapL3B_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL3B_[i]->setAxisTitle("channel", 1);
      meShapeMapL3B_[i]->setAxisTitle("sample", 2);
      meShapeMapL3B_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL3B_[i], i+1);
      sprintf(histo, "EELT amplitude %s L3B", Numbers::sEE(i+1).c_str());
      meAmplMapL3B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL3B_[i]->setAxisTitle("jx", 1);
      meAmplMapL3B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL3B_[i], i+1);
      sprintf(histo, "EELT timing %s L3B", Numbers::sEE(i+1).c_str());
      meTimeMapL3B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL3B_[i]->setAxisTitle("jx", 1);
      meTimeMapL3B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL3B_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L3B", Numbers::sEE(i+1).c_str());
      meAmplPNMapL3B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL3B_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL3B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL3B_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EELT shape %s L4A", Numbers::sEE(i+1).c_str());
      meShapeMapL4A_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL4A_[i]->setAxisTitle("channel", 1);
      meShapeMapL4A_[i]->setAxisTitle("sample", 2);
      meShapeMapL4A_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL4A_[i], i+1);
      sprintf(histo, "EELT amplitude %s L4A", Numbers::sEE(i+1).c_str());
      meAmplMapL4A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL4A_[i]->setAxisTitle("jx", 1);
      meAmplMapL4A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL4A_[i], i+1);
      sprintf(histo, "EELT timing %s L4A", Numbers::sEE(i+1).c_str());
      meTimeMapL4A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL4A_[i]->setAxisTitle("jx", 1);
      meTimeMapL4A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL4A_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L4A", Numbers::sEE(i+1).c_str());
      meAmplPNMapL4A_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL4A_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL4A_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL4A_[i], i+1);

      sprintf(histo, "EELT shape %s L4B", Numbers::sEE(i+1).c_str());
      meShapeMapL4B_[i] = dqmStore_->bookProfile2D(histo, histo, 850, 0., 850., 10, 0., 10., 4096, 0., 4096., "s");
      meShapeMapL4B_[i]->setAxisTitle("channel", 1);
      meShapeMapL4B_[i]->setAxisTitle("sample", 2);
      meShapeMapL4B_[i]->setAxisTitle("amplitude", 3);
      dqmStore_->tag(meShapeMapL4B_[i], i+1);
      sprintf(histo, "EELT amplitude %s L4B", Numbers::sEE(i+1).c_str());
      meAmplMapL4B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplMapL4B_[i]->setAxisTitle("jx", 1);
      meAmplMapL4B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplMapL4B_[i], i+1);
      sprintf(histo, "EELT timing %s L4B", Numbers::sEE(i+1).c_str());
      meTimeMapL4B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 250, 0., 10., "s");
      meTimeMapL4B_[i]->setAxisTitle("jx", 1);
      meTimeMapL4B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meTimeMapL4B_[i], i+1);
      sprintf(histo, "EELT amplitude over PN %s L4B", Numbers::sEE(i+1).c_str());
      meAmplPNMapL4B_[i] = dqmStore_->bookProfile2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50., 4096, 0., 4096.*12., "s");
      meAmplPNMapL4B_[i]->setAxisTitle("jx", 1);
      meAmplPNMapL4B_[i]->setAxisTitle("jy", 2);
      dqmStore_->tag(meAmplPNMapL4B_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G01 L1", Numbers::sEE(i+1).c_str());
      mePnAmplMapG01L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG01L1_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG01L1_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG01L1_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G01 L1", Numbers::sEE(i+1).c_str());
      mePnPedMapG01L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG01L1_[i]->setAxisTitle("channel", 1);
      mePnPedMapG01L1_[i]->setAxisTitle("pedestal", 2);
      dqmStore_->tag(mePnPedMapG01L1_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G16 L1", Numbers::sEE(i+1).c_str());
      mePnAmplMapG16L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG16L1_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG16L1_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG16L1_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G16 L1", Numbers::sEE(i+1).c_str());
      mePnPedMapG16L1_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG16L1_[i]->setAxisTitle("channel", 1);
      mePnPedMapG16L1_[i]->setAxisTitle("pedestal", 2); 
      dqmStore_->tag(mePnPedMapG16L1_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G01 L2", Numbers::sEE(i+1).c_str());
      mePnAmplMapG01L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG01L2_[i]->setAxisTitle("amplitude", 2);
      mePnAmplMapG01L2_[i]->setAxisTitle("channel", 1);
      dqmStore_->tag(mePnAmplMapG01L2_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G01 L2", Numbers::sEE(i+1).c_str());
      mePnPedMapG01L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG01L2_[i]->setAxisTitle("channel", 1);
      mePnPedMapG01L2_[i]->setAxisTitle("pedestal", 2);
      dqmStore_->tag(mePnPedMapG01L2_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G16 L2", Numbers::sEE(i+1).c_str());
      mePnAmplMapG16L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG16L2_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG16L2_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG16L2_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G16 L2", Numbers::sEE(i+1).c_str());
      mePnPedMapG16L2_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG16L2_[i]->setAxisTitle("channel", 1);
      mePnPedMapG16L2_[i]->setAxisTitle("pedestal", 2); 
      dqmStore_->tag(mePnPedMapG16L2_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G01 L3", Numbers::sEE(i+1).c_str());
      mePnAmplMapG01L3_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG01L3_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG01L3_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG01L3_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G01 L3", Numbers::sEE(i+1).c_str());
      mePnPedMapG01L3_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG01L3_[i]->setAxisTitle("channel", 1);
      mePnPedMapG01L3_[i]->setAxisTitle("pedestal", 2);
      dqmStore_->tag(mePnPedMapG01L3_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G16 L3", Numbers::sEE(i+1).c_str());
      mePnAmplMapG16L3_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG16L3_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG16L3_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG16L3_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G16 L3", Numbers::sEE(i+1).c_str());
      mePnPedMapG16L3_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG16L3_[i]->setAxisTitle("channel", 1);
      mePnPedMapG16L3_[i]->setAxisTitle("pedestal", 2); 
      dqmStore_->tag(mePnPedMapG16L3_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G01 L4", Numbers::sEE(i+1).c_str());
      mePnAmplMapG01L4_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG01L4_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG01L4_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG01L4_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G01 L4", Numbers::sEE(i+1).c_str());
      mePnPedMapG01L4_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG01L4_[i]->setAxisTitle("channel", 1);
      mePnPedMapG01L4_[i]->setAxisTitle("pedestal", 2);
      dqmStore_->tag(mePnPedMapG01L4_[i], i+1);
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EEPDT PNs amplitude %s G16 L4", Numbers::sEE(i+1).c_str());
      mePnAmplMapG16L4_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnAmplMapG16L4_[i]->setAxisTitle("channel", 1);
      mePnAmplMapG16L4_[i]->setAxisTitle("amplitude", 2);
      dqmStore_->tag(mePnAmplMapG16L4_[i], i+1);
      sprintf(histo, "EEPDT PNs pedestal %s G16 L4", Numbers::sEE(i+1).c_str());
      mePnPedMapG16L4_[i] = dqmStore_->bookProfile(histo, histo, 10, 0., 10., 4096, 0., 4096., "s");
      mePnPedMapG16L4_[i]->setAxisTitle("channel", 1);
      mePnPedMapG16L4_[i]->setAxisTitle("pedestal", 2); 
      dqmStore_->tag(mePnPedMapG16L4_[i], i+1);
    }

  }

}

void EELaserTask::cleanup(void){

  if ( ! init_ ) return;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1");
    for (int i = 0; i < 18; i++) {
      if ( meShapeMapL1A_[i] )  dqmStore_->removeElement( meShapeMapL1A_[i]->getName() );
      meShapeMapL1A_[i] = 0;
      if ( meAmplMapL1A_[i] ) dqmStore_->removeElement( meAmplMapL1A_[i]->getName() );
      meAmplMapL1A_[i] = 0;
      if ( meTimeMapL1A_[i] ) dqmStore_->removeElement( meTimeMapL1A_[i]->getName() );
      meTimeMapL1A_[i] = 0;
      if ( meAmplPNMapL1A_[i] ) dqmStore_->removeElement( meAmplPNMapL1A_[i]->getName() );
      meAmplPNMapL1A_[i] = 0;

      if ( meShapeMapL1B_[i] )  dqmStore_->removeElement( meShapeMapL1B_[i]->getName() );
      meShapeMapL1B_[i] = 0;
      if ( meAmplMapL1B_[i] ) dqmStore_->removeElement( meAmplMapL1B_[i]->getName() );
      meAmplMapL1B_[i] = 0;
      if ( meTimeMapL1B_[i] ) dqmStore_->removeElement( meTimeMapL1B_[i]->getName() );
      meTimeMapL1B_[i] = 0;
      if ( meAmplPNMapL1B_[i] ) dqmStore_->removeElement( meAmplPNMapL1B_[i]->getName() );
      meAmplPNMapL1B_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2");
    for (int i = 0; i < 18; i++) {
      if ( meShapeMapL2A_[i] )  dqmStore_->removeElement( meShapeMapL2A_[i]->getName() );
      meShapeMapL2A_[i] = 0;
      if ( meAmplMapL2A_[i] ) dqmStore_->removeElement( meAmplMapL2A_[i]->getName() );
      meAmplMapL2A_[i] = 0;
      if ( meTimeMapL2A_[i] ) dqmStore_->removeElement( meTimeMapL2A_[i]->getName() );
      meTimeMapL2A_[i] = 0;
      if ( meAmplPNMapL2A_[i] ) dqmStore_->removeElement( meAmplPNMapL2A_[i]->getName() );
      meAmplPNMapL2A_[i] = 0;

      if ( meShapeMapL2B_[i] )  dqmStore_->removeElement( meShapeMapL2B_[i]->getName() );
      meShapeMapL2B_[i] = 0;
      if ( meAmplMapL2B_[i] ) dqmStore_->removeElement( meAmplMapL2B_[i]->getName() );
      meAmplMapL2B_[i] = 0;
      if ( meTimeMapL2B_[i] ) dqmStore_->removeElement( meTimeMapL2B_[i]->getName() );
      meTimeMapL2B_[i] = 0;
      if ( meAmplPNMapL2B_[i] ) dqmStore_->removeElement( meAmplPNMapL2B_[i]->getName() );
      meAmplPNMapL2B_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3");
    for (int i = 0; i < 18; i++) {
      if ( meShapeMapL3A_[i] )  dqmStore_->removeElement( meShapeMapL3A_[i]->getName() );
      meShapeMapL3A_[i] = 0;
      if ( meAmplMapL3A_[i] ) dqmStore_->removeElement( meAmplMapL3A_[i]->getName() );
      meAmplMapL3A_[i] = 0;
      if ( meTimeMapL3A_[i] ) dqmStore_->removeElement( meTimeMapL3A_[i]->getName() );
      meTimeMapL3A_[i] = 0;
      if ( meAmplPNMapL3A_[i] ) dqmStore_->removeElement( meAmplPNMapL3A_[i]->getName() );
      meAmplPNMapL3A_[i] = 0;

      if ( meShapeMapL3B_[i] )  dqmStore_->removeElement( meShapeMapL3B_[i]->getName() );
      meShapeMapL3B_[i] = 0;
      if ( meAmplMapL3B_[i] ) dqmStore_->removeElement( meAmplMapL3B_[i]->getName() );
      meAmplMapL3B_[i] = 0;
      if ( meTimeMapL3B_[i] ) dqmStore_->removeElement( meTimeMapL3B_[i]->getName() );
      meTimeMapL3B_[i] = 0;
      if ( meAmplPNMapL3B_[i] ) dqmStore_->removeElement( meAmplPNMapL3B_[i]->getName() );
      meAmplPNMapL3B_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4");
    for (int i = 0; i < 18; i++) {
      if ( meShapeMapL4A_[i] )  dqmStore_->removeElement( meShapeMapL4A_[i]->getName() );
      meShapeMapL4A_[i] = 0;
      if ( meAmplMapL4A_[i] ) dqmStore_->removeElement( meAmplMapL4A_[i]->getName() );
      meAmplMapL4A_[i] = 0;
      if ( meTimeMapL4A_[i] ) dqmStore_->removeElement( meTimeMapL4A_[i]->getName() );
      meTimeMapL4A_[i] = 0;
      if ( meAmplPNMapL4A_[i] ) dqmStore_->removeElement( meAmplPNMapL4A_[i]->getName() );
      meAmplPNMapL4A_[i] = 0;

      if ( meShapeMapL4B_[i] )  dqmStore_->removeElement( meShapeMapL4B_[i]->getName() );
      meShapeMapL4B_[i] = 0;
      if ( meAmplMapL4B_[i] ) dqmStore_->removeElement( meAmplMapL4B_[i]->getName() );
      meAmplMapL4B_[i] = 0;
      if ( meTimeMapL4B_[i] ) dqmStore_->removeElement( meTimeMapL4B_[i]->getName() );
      meTimeMapL4B_[i] = 0;
      if ( meAmplPNMapL4B_[i] ) dqmStore_->removeElement( meAmplPNMapL4B_[i]->getName() );
      meAmplPNMapL4B_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG01L1_[i] ) dqmStore_->removeElement( mePnAmplMapG01L1_[i]->getName() );
      mePnAmplMapG01L1_[i] = 0;
      if ( mePnPedMapG01L1_[i] ) dqmStore_->removeElement( mePnPedMapG01L1_[i]->getName() );
      mePnPedMapG01L1_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser1/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG16L1_[i] ) dqmStore_->removeElement( mePnAmplMapG16L1_[i]->getName() );
      mePnAmplMapG16L1_[i] = 0;
      if ( mePnPedMapG16L1_[i] ) dqmStore_->removeElement( mePnPedMapG16L1_[i]->getName() );
      mePnPedMapG16L1_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG01L2_[i] ) dqmStore_->removeElement( mePnAmplMapG01L2_[i]->getName() );
      mePnAmplMapG01L2_[i] = 0;
      if ( mePnPedMapG01L2_[i] ) dqmStore_->removeElement( mePnPedMapG01L2_[i]->getName() );
      mePnPedMapG01L2_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser2/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG16L2_[i] ) dqmStore_->removeElement( mePnAmplMapG16L2_[i]->getName() );
      mePnAmplMapG16L2_[i] = 0;
      if ( mePnPedMapG16L2_[i] ) dqmStore_->removeElement( mePnPedMapG16L2_[i]->getName() );
      mePnPedMapG16L2_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG01L3_[i] ) dqmStore_->removeElement( mePnAmplMapG01L3_[i]->getName() );
      mePnAmplMapG01L3_[i] = 0;
      if ( mePnPedMapG01L3_[i] ) dqmStore_->removeElement( mePnPedMapG01L3_[i]->getName() );
      mePnPedMapG01L3_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser3/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG16L3_[i] ) dqmStore_->removeElement( mePnAmplMapG16L3_[i]->getName() );
      mePnAmplMapG16L3_[i] = 0;
      if ( mePnPedMapG16L3_[i] ) dqmStore_->removeElement( mePnPedMapG16L3_[i]->getName() );
      mePnPedMapG16L3_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN");

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN/Gain01");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG01L4_[i] ) dqmStore_->removeElement( mePnAmplMapG01L4_[i]->getName() );
      mePnAmplMapG01L4_[i] = 0;
      if ( mePnPedMapG01L4_[i] ) dqmStore_->removeElement( mePnPedMapG01L4_[i]->getName() );
      mePnPedMapG01L4_[i] = 0;
    }

    dqmStore_->setCurrentFolder(prefixME_ + "/EELaserTask/Laser4/PN/Gain16");
    for (int i = 0; i < 18; i++) {
      if ( mePnAmplMapG16L4_[i] ) dqmStore_->removeElement( mePnAmplMapG16L4_[i]->getName() );
      mePnAmplMapG16L4_[i] = 0;
      if ( mePnPedMapG16L4_[i] ) dqmStore_->removeElement( mePnPedMapG16L4_[i]->getName() );
      mePnPedMapG16L4_[i] = 0;
    }

  }

  init_ = false;

}

void EELaserTask::endJob(void){

  LogInfo("EELaserTask") << "analyzed " << ievt_ << " events";

  if ( enableCleanup_ ) this->cleanup();

}

void EELaserTask::analyze(const Event& e, const EventSetup& c){

  bool enable = false;
  map<int, EcalDCCHeaderBlock> dccMap;

  Handle<EcalRawDataCollection> dcchs;

  if ( e.getByLabel(EcalRawDataCollection_, dcchs) ) {

    for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin(); dcchItr != dcchs->end(); ++dcchItr ) {

      EcalDCCHeaderBlock dcch = (*dcchItr);

      if ( Numbers::subDet( dcch ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( dcch, EcalEndcap );

      map<int, EcalDCCHeaderBlock>::iterator i = dccMap.find( ism );
      if ( i != dccMap.end() ) continue;

      dccMap[ ism ] = dcch;

      if ( dcch.getRunType() == EcalDCCHeaderBlock::LASER_STD ||
           dcch.getRunType() == EcalDCCHeaderBlock::LASER_GAP ) enable = true;

    }

  } else {

    LogWarning("EELaserTask") << EcalRawDataCollection_ << " not available";

  }

  if ( ! enable ) return;

  if ( ! init_ ) this->setup();

  ievt_++;

  Handle<EEDigiCollection> digis;

  if ( e.getByLabel(EEDigiCollection_, digis) ) {

    int need = digis->size();
    LogDebug("EELaserTask") << "event " << ievt_ << " digi collection size " << need;

    for ( EEDigiCollection::const_iterator digiItr = digis->begin(); digiItr != digis->end(); ++digiItr ) {

      EEDataFrame dataframe = (*digiItr);
      EEDetId id = dataframe.id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      map<int, EcalDCCHeaderBlock>::iterator i = dccMap.find(ism);
      if ( i == dccMap.end() ) continue;

      if ( ! ( dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_STD ||
               dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_GAP ) ) continue;


      if ( dccMap[ism].getRtHalf() != Numbers::RtHalf(id) ) continue;

      LogDebug("EELaserTask") << " det id = " << id;
      LogDebug("EELaserTask") << " sm, ix, iy " << ism << " " << ix << " " << iy;

      int ic = Numbers::icEE(ism, ix, iy);

      for (int i = 0; i < 10; i++) {

        EcalMGPASample sample = dataframe.sample(i);
        int adc = sample.adc();
        float gain = 1.;

        MonitorElement* meShapeMap = 0;

        if ( sample.gainId() == 1 ) gain = 1./12.;
        if ( sample.gainId() == 2 ) gain = 1./ 6.;
        if ( sample.gainId() == 3 ) gain = 1./ 1.;

        if ( dccMap[ism].getRtHalf() == 0 ) {

          if ( dccMap[ism].getEventSettings().wavelength == 0 ) meShapeMap = meShapeMapL1A_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 1 ) meShapeMap = meShapeMapL2A_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 2 ) meShapeMap = meShapeMapL3A_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 3 ) meShapeMap = meShapeMapL4A_[ism-1];

        } else if ( dccMap[ism].getRtHalf() == 1 ) {

          if ( dccMap[ism].getEventSettings().wavelength == 0 ) meShapeMap = meShapeMapL1B_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 1 ) meShapeMap = meShapeMapL2B_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 2 ) meShapeMap = meShapeMapL3B_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 3 ) meShapeMap = meShapeMapL4B_[ism-1];

        } else {

          LogWarning("EELaserTask") << " RtHalf = " << dccMap[ism].getRtHalf();

        }

//        float xval = float(adc) * gain;
        float xval = float(adc);

        if ( meShapeMap ) meShapeMap->Fill(ic - 0.5, i + 0.5, xval);

      }

    }

  } else {

    LogWarning("EELaserTask") << EEDigiCollection_ << " not available";

  }

  float adcA[18];
  float adcB[18];

  for ( int i = 0; i < 18; i++ ) {
    adcA[i] = 0.;
    adcB[i] = 0.;
  }

  Handle<EcalPnDiodeDigiCollection> pns;

  if ( e.getByLabel(EcalPnDiodeDigiCollection_, pns) ) {

    int nep = pns->size();
    LogDebug("EELaserTask") << "event " << ievt_ << " pns collection size " << nep;

    for ( EcalPnDiodeDigiCollection::const_iterator pnItr = pns->begin(); pnItr != pns->end(); ++pnItr ) {

      EcalPnDiodeDigi pn = (*pnItr);
      EcalPnDiodeDetId id = pn.id();

      if ( Numbers::subDet( id ) != EcalEndcap ) continue;

      int ism = Numbers::iSM( id );

      int num = id.iPnId();

      map<int, EcalDCCHeaderBlock>::iterator i = dccMap.find(ism);
      if ( i == dccMap.end() ) continue;

      if ( ! ( dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_STD ||
               dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_GAP ) ) continue;

      LogDebug("EELaserTask") << " det id = " << id;
      LogDebug("EELaserTask") << " sm, num " << ism << " " << num;

      float xvalped = 0.;

      for (int i = 0; i < 4; i++) {

        EcalFEMSample sample = pn.sample(i);
        int adc = sample.adc();

        MonitorElement* mePNPed = 0;

        if ( sample.gainId() == 0 ) {
          if ( dccMap[ism].getEventSettings().wavelength == 0 ) mePNPed = mePnPedMapG01L1_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 1 ) mePNPed = mePnPedMapG01L2_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 2 ) mePNPed = mePnPedMapG01L3_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 3 ) mePNPed = mePnPedMapG01L4_[ism-1];
        }
        if ( sample.gainId() == 1 ) {
          if ( dccMap[ism].getEventSettings().wavelength == 0 ) mePNPed = mePnPedMapG16L1_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 1 ) mePNPed = mePnPedMapG16L2_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 2 ) mePNPed = mePnPedMapG16L3_[ism-1];
          if ( dccMap[ism].getEventSettings().wavelength == 3 ) mePNPed = mePnPedMapG16L4_[ism-1];
        }

        float xval = float(adc);

        if ( mePNPed ) mePNPed->Fill(num - 0.5, xval);

        xvalped = xvalped + xval;

      }

      xvalped = xvalped / 4;

      float xvalmax = 0.;

      MonitorElement* mePN = 0;

      for (int i = 0; i < 50; i++) {

        EcalFEMSample sample = pn.sample(i);
        int adc = sample.adc();

        float xval = float(adc);

        if ( xval >= xvalmax ) xvalmax = xval;

      }

      xvalmax = xvalmax - xvalped;

      if ( pn.sample(0).gainId() == 0 ) {
        if ( dccMap[ism].getEventSettings().wavelength == 0 ) mePN = mePnAmplMapG01L1_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 1 ) mePN = mePnAmplMapG01L2_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 2 ) mePN = mePnAmplMapG01L3_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 3 ) mePN = mePnAmplMapG01L4_[ism-1];
      }
      if ( pn.sample(0).gainId() == 1 ) {
        if ( dccMap[ism].getEventSettings().wavelength == 0 ) mePN = mePnAmplMapG16L1_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 1 ) mePN = mePnAmplMapG16L2_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 2 ) mePN = mePnAmplMapG16L3_[ism-1];
        if ( dccMap[ism].getEventSettings().wavelength == 3 ) mePN = mePnAmplMapG16L4_[ism-1];
      }

      if ( mePN ) mePN->Fill(num - 0.5, xvalmax);

      if ( num == 1 ) adcA[ism-1] = xvalmax;
      if ( num == 6 ) adcB[ism-1] = xvalmax;

    }

  } else {

    LogWarning("EELaserTask") << EcalPnDiodeDigiCollection_ << " not available";

  }

  Handle<EcalUncalibratedRecHitCollection> hits;

  if ( e.getByLabel(EcalUncalibratedRecHitCollection_, hits) ) {

    int neh = hits->size();
    LogDebug("EELaserTask") << "event " << ievt_ << " hits collection size " << neh;

    for ( EcalUncalibratedRecHitCollection::const_iterator hitItr = hits->begin(); hitItr != hits->end(); ++hitItr ) {

      EcalUncalibratedRecHit hit = (*hitItr);
      EEDetId id = hit.id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      map<int, EcalDCCHeaderBlock>::iterator i = dccMap.find(ism);
      if ( i == dccMap.end() ) continue;

      if ( ! ( dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_STD ||
               dccMap[ism].getRunType() == EcalDCCHeaderBlock::LASER_GAP ) ) continue;

      if ( dccMap[ism].getRtHalf() != Numbers::RtHalf(id) ) continue;

      LogDebug("EELaserTask") << " det id = " << id;
      LogDebug("EELaserTask") << " sm, ix, iy " << ism << " " << ix << " " << iy;

      MonitorElement* meAmplMap = 0;
      MonitorElement* meTimeMap = 0;
      MonitorElement* meAmplPNMap = 0;

      if ( dccMap[ism].getRtHalf() == 0 ) {

        if ( dccMap[ism].getEventSettings().wavelength == 0 ) {
          meAmplMap = meAmplMapL1A_[ism-1];
          meTimeMap = meTimeMapL1A_[ism-1];
          meAmplPNMap = meAmplPNMapL1A_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 1 ) {
          meAmplMap = meAmplMapL2A_[ism-1];
          meTimeMap = meTimeMapL2A_[ism-1];
          meAmplPNMap = meAmplPNMapL2A_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 2 ) {
          meAmplMap = meAmplMapL3A_[ism-1];
          meTimeMap = meTimeMapL3A_[ism-1];
          meAmplPNMap = meAmplPNMapL3A_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 3 ) {
          meAmplMap = meAmplMapL4A_[ism-1];
          meTimeMap = meTimeMapL4A_[ism-1];
          meAmplPNMap = meAmplPNMapL4A_[ism-1];
        }

      } else if ( dccMap[ism].getRtHalf() == 1 ) { 

        if ( dccMap[ism].getEventSettings().wavelength == 0 ) {
          meAmplMap = meAmplMapL1B_[ism-1];
          meTimeMap = meTimeMapL1B_[ism-1];
          meAmplPNMap = meAmplPNMapL1B_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 1 ) {
          meAmplMap = meAmplMapL2B_[ism-1];
          meTimeMap = meTimeMapL2B_[ism-1];
          meAmplPNMap = meAmplPNMapL2B_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 2 ) {
          meAmplMap = meAmplMapL3B_[ism-1];
          meTimeMap = meTimeMapL3B_[ism-1];
          meAmplPNMap = meAmplPNMapL3B_[ism-1];
        }
        if ( dccMap[ism].getEventSettings().wavelength == 3 ) {
          meAmplMap = meAmplMapL4B_[ism-1];
          meTimeMap = meTimeMapL4B_[ism-1];
          meAmplPNMap = meAmplPNMapL4B_[ism-1];
        }

      } else {

        LogWarning("EELaserTask") << " RtHalf = " << dccMap[ism].getRtHalf();

      }

      float xval = hit.amplitude();
      if ( xval <= 0. ) xval = 0.0;
      float yval = hit.jitter() + 5.0;
      if ( yval <= 0. ) yval = 0.0;
      float zval = hit.pedestal();
      if ( zval <= 0. ) zval = 0.0;

      LogDebug("EELaserTask") << " hit amplitude " << xval;
      LogDebug("EELaserTask") << " hit jitter " << yval;
      LogDebug("EELaserTask") << " hit pedestal " << zval;

      if ( meAmplMap ) meAmplMap->Fill(xix, xiy, xval);

      if ( xval > 16. ) {
        if ( meTimeMap ) meTimeMap->Fill(xix, xiy, yval);
      }

      float wval = 0.;

      if ( dccMap[ism].getRtHalf() == 0 ) {

        if ( adcA[ism-1] != 0. ) wval = xval / adcA[ism-1];

      } else if ( dccMap[ism].getRtHalf() == 1 ) {

        if ( adcB[ism-1] != 0. ) wval = xval / adcB[ism-1];

      } else {

        LogWarning("EELaserTask") << " RtHalf = " << dccMap[ism].getRtHalf();

      }

      LogDebug("EELaserTask") << " hit amplitude over PN " << wval;

      if ( meAmplPNMap ) meAmplPNMap->Fill(xix, xiy, wval);

    }

  } else {

    LogWarning("EELaserTask") << EcalUncalibratedRecHitCollection_ << " not available";

  }

}

