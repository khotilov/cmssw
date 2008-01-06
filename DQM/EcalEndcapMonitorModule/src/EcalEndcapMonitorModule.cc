/*
 * \file EcalEndcapMonitorModule.cc
 *
 * $Date: 2008/01/05 11:07:25 $
 * $Revision: 1.36 $
 * \author G. Della Ricca
 * \author G. Franzoni
 *
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <DQM/EcalCommon/interface/Numbers.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>

#include <DQM/EcalEndcapMonitorModule/interface/EcalEndcapMonitorModule.h>

using namespace cms;
using namespace edm;
using namespace std;

EcalEndcapMonitorModule::EcalEndcapMonitorModule(const ParameterSet& ps){

  init_ = false;

  EcalTBEventHeader_ = ps.getParameter<edm::InputTag>("EcalTBEventHeader");
  EcalRawDataCollection_ = ps.getParameter<edm::InputTag>("EcalRawDataCollection");
  EEDigiCollection_ = ps.getParameter<edm::InputTag>("EEDigiCollection");
  EcalUncalibratedRecHitCollection_ = ps.getParameter<edm::InputTag>("EcalUncalibratedRecHitCollection");
  EcalTrigPrimDigiCollection_ = ps.getParameter<edm::InputTag>("EcalTrigPrimDigiCollection");

  cout << endl;
  cout << " *** Ecal Endcap Generic Monitor ***" << endl;
  cout << endl;

  // this should come from the event header
  runNumber_ = ps.getUntrackedParameter<int>("runNumber", 0);

  fixedRunNumber_ = false;
  if ( runNumber_ != 0 ) fixedRunNumber_ = true;

  if ( fixedRunNumber_ ) {
    LogInfo("EcalBarrelMonitor") << " using fixed Run Number = " << runNumber_ << endl;
  }

  // this should come from the event header
  evtNumber_ = 0;

  // this should come from the EcalEndcap event header
  runType_ = ps.getUntrackedParameter<int>("runType", -1);
  evtType_ = runType_;

  fixedRunType_ = false;
  if ( runType_ != -1 ) fixedRunType_ = true;

  if ( fixedRunType_) {
    LogInfo("EcalBarrelMonitor") << " using fixed Run Type = " << runType_ << endl;
  }

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  if ( verbose_ ) {
    LogInfo("EcalEndcapMonitor") << " verbose switch is ON";
  } else {
    LogInfo("EcalEndcapMonitor") << " verbose switch is OFF";
  }

  // get hold of back-end interface
  dbe_ = Service<DaqMonitorBEInterface>().operator->();

  if ( dbe_ ) {
    if ( verbose_ ) {
      dbe_->setVerbose(1);
    } else {
      dbe_->setVerbose(0);
    }
  }

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  if ( enableCleanup_ ) {
    LogInfo("EcalBarrelMonitor") << " enableCleanup switch is ON";
  } else {
    LogInfo("EcalBarrelMonitor") << " enableCleanup switch is OFF";
  }

  // MonitorDaemon switch
  enableMonitorDaemon_ = ps.getUntrackedParameter<bool>("enableMonitorDaemon", true);

  if ( enableMonitorDaemon_ ) {
    LogInfo("EcalEndcapMonitor") << " enableMonitorDaemon switch is ON";
    Service<MonitorDaemon> daemon;
    daemon.operator->();
  } else {
    LogInfo("EcalEndcapMonitor") << " enableMonitorDaemon switch is OFF";
  }

  // EventDisplay switch
  enableEventDisplay_ = ps.getUntrackedParameter<bool>("enableEventDisplay", false);

  meStatus_ = 0;
  meRun_ = 0;
  meEvt_ = 0;
  meRunType_ = 0;
  meEvtType_ = 0;

  meEEDCC_ = 0;

  for (int i = 0; i < 2; i++) {
    meEEdigis_[i] = 0;
    meEEhits_[i] = 0;
    meEEtpdigis_[i] = 0;
  }

  for (int i = 0; i < 18; i++) {
    meEvent_[i] = 0;
  }

}

EcalEndcapMonitorModule::~EcalEndcapMonitorModule(){

}

void EcalEndcapMonitorModule::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EcalInfo");
    dbe_->rmdir("EcalEndcap/EcalInfo");
    if ( enableEventDisplay_ ) {
      dbe_->setCurrentFolder("EcalEndcap/EcalEvent");
      dbe_->rmdir("EcalEndcap/EcalEvent");
    }
  }

}

void EcalEndcapMonitorModule::setup(void){

  init_ = true;

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EcalInfo");

    meStatus_ = dbe_->bookInt("STATUS");

    meRun_ = dbe_->bookInt("RUN");
    meEvt_ = dbe_->bookInt("EVT");

    meRunType_ = dbe_->bookInt("RUNTYPE");
    meEvtType_ = dbe_->book1D("EVTTYPE", "EVTTYPE", 31, -1., 30.);
    meEvtType_->setBinLabel(1, "UNKNOWN", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::COSMIC, "COSMIC", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::BEAMH4, "BEAMH4", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::BEAMH2, "BEAMH2", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::MTCC, "MTCC", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LASER_STD, "LASER_STD", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LASER_POWER_SCAN, "LASER_POWER_SCAN", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LASER_DELAY_SCAN, "LASER_DELAY_SCAN", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::TESTPULSE_SCAN_MEM, "TESTPULSE_SCAN_MEM", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::TESTPULSE_MGPA, "TESTPULSE_MGPA", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PEDESTAL_STD, "PEDESTAL_STD", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PEDESTAL_OFFSET_SCAN, "PEDESTAL_OFFSET_SCAN", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PEDESTAL_25NS_SCAN, "PEDESTAL_25NS_SCAN", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LED_STD, "LED_STD", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PHYSICS_GLOBAL, "PHYSICS_GLOBAL", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::COSMICS_GLOBAL, "COSMICS_GLOBAL", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::HALO_GLOBAL, "HALO_GLOBAL", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LASER_GAP, "LASER_GAP", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::TESTPULSE_GAP, "TESTPULSE_GAP");
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PEDESTAL_GAP, "PEDESTAL_GAP");
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::LED_GAP, "LED_GAP", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::PHYSICS_LOCAL, "PHYSICS_LOCAL", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::COSMICS_LOCAL, "COSMICS_LOCAL", 1);
    meEvtType_->setBinLabel(2+EcalDCCHeaderBlock::HALO_LOCAL, "HALO_LOCAL", 1);
  }

  // unknown
  if ( meStatus_ ) meStatus_->Fill(-1);

  if ( meRun_ ) meRun_->Fill(-1);
  if ( meEvt_ ) meEvt_->Fill(-1);

  if ( meRunType_ ) meRunType_->Fill(-1);

  // this should give enough time to our control MEs to reach the Collector,
  // and then hopefully the Client

  if ( enableMonitorDaemon_ ) sleep(5);

  Char_t histo[20];

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EcalInfo");

    meEEDCC_ = dbe_->book1D("EEMM DCC", "EEMM DCC", 18, 1, 19.);
    for (int i = 0; i < 18; i++) {
      meEEDCC_->setBinLabel(i+1, Numbers::sEE(i+1).c_str(), 1);
    }

    meEEdigis_[0] = dbe_->book1D("EEMM digis number", "EEMM digis number", 100, 0., 13299.);

    meEEdigis_[1] = dbe_->bookProfile("EEMM digis number profile", "EEMM digis number profile", 18, 1, 19., 850, 0., 851., "s");
    for (int i = 0; i < 18; i++) {
      meEEdigis_[1]->setBinLabel(i+1, Numbers::sEE(i+1).c_str(), 1);
    }

    meEEhits_[0] = dbe_->book1D("EEMM hits number", "EEMM hits number", 100, 0., 13299.);

    meEEhits_[1] = dbe_->bookProfile("EEMM hits number profile", "EEMM hits number profile", 18, 1, 19., 850, 0., 851., "s");
    for (int i = 0; i < 18; i++) {
      meEEhits_[1]->setBinLabel(i+1, Numbers::sEE(i+1).c_str(), 1);
    }

    meEEtpdigis_[0] = dbe_->book1D("EEMM TP digis number", "EEMM TP digis number", 100, 0., 631.);

    meEEtpdigis_[1] = dbe_->bookProfile("EEMM TP digis number profile", "EEMM TP digis number profile", 18, 1, 19., 34, 0., 35., "s");
    for (int i = 0; i < 18; i++) {
      meEEtpdigis_[1]->setBinLabel(i+1, Numbers::sEE(i+1).c_str(), 1);
    }

    if ( enableEventDisplay_ ) {
      dbe_->setCurrentFolder("EcalEndcap/EcalEvent");
      for (int i = 0; i < 18; i++) {
        sprintf(histo, "EEMM event %s", Numbers::sEE(i+1).c_str());
        meEvent_[i] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(i+1)+0., Numbers::ix0EE(i+1)+50., 50, Numbers::iy0EE(i+1)+0., Numbers::iy0EE(i+1)+50.);
        meEvent_[i]->setAxisTitle("ix", 1);
        meEvent_[i]->setAxisTitle("iy", 2);
        dbe_->tag(meEvent_[i], i+1);
        if ( meEvent_[i] ) meEvent_[i]->setResetMe(true);
      }
    }

  }

}

void EcalEndcapMonitorModule::cleanup(void){

  if ( ! enableCleanup_ ) return;

  if ( dbe_ ) {

    dbe_->setCurrentFolder("EcalEndcap/EcalInfo");

    if ( meStatus_ ) dbe_->removeElement( meStatus_->getName() );
    meStatus_ = 0;

    if ( meRun_ ) dbe_->removeElement( meRun_->getName() );
    meRun_ = 0;

    if ( meEvt_ ) dbe_->removeElement( meEvt_->getName() );
    meEvt_ = 0;

    if ( meRunType_ ) dbe_->removeElement( meRunType_->getName() );
    meRunType_ = 0;

    if ( meEvtType_ ) dbe_->removeElement( meEvtType_->getName() );
    meEvtType_ = 0;

    if ( meEEDCC_ ) dbe_->removeElement( meEEDCC_->getName() );
    meEEDCC_ = 0;

    for (int i = 0; i < 2; i++) {

      if ( meEEdigis_[i] ) dbe_->removeElement( meEEdigis_[i]->getName() );
      meEEdigis_[i] = 0;

      if ( meEEhits_[i] ) dbe_->removeElement( meEEhits_[i]->getName() );
      meEEhits_[i] = 0;

      if ( meEEtpdigis_[i] ) dbe_->removeElement( meEEtpdigis_[i]->getName() );
      meEEtpdigis_[i] = 0;

    }

    if ( enableEventDisplay_ ) {

      dbe_->setCurrentFolder("EcalEndcap/EcalEvent");
      for (int i = 0; i < 18; i++) {

        if ( meEvent_[i] ) dbe_->removeElement( meEvent_[i]->getName() );
        meEvent_[i] = 0;

      }

    }

  }

  init_ = false;

}

void EcalEndcapMonitorModule::endJob(void) {

  LogInfo("EcalEndcapMonitor") << "analyzed " << ievt_ << " events";

  // end-of-run
  if ( meStatus_ ) meStatus_->Fill(2);

  if ( meRun_ ) meRun_->Fill(runNumber_);
  if ( meEvt_ ) meEvt_->Fill(evtNumber_);

  // this should give enough time to meStatus_ to reach the Collector,
  // and then hopefully the Client, and to allow the Client to complete

  // we should always sleep at least a little ...

  if ( enableMonitorDaemon_ ) sleep(5);

  if ( init_ ) this->cleanup();

}

void EcalEndcapMonitorModule::analyze(const Event& e, const EventSetup& c){

  Numbers::initGeometry(c);

  if ( ! init_ ) this->setup();

  ievt_++;

  LogInfo("EcalEndcapMonitor") << "processing event " << ievt_;

  if ( ! fixedRunNumber_ ) runNumber_ = e.id().run();

  evtNumber_ = e.id().event();

  map<int, EcalDCCHeaderBlock> dccMap;

  Handle<EcalRawDataCollection> dcchs;

  if ( e.getByLabel(EcalRawDataCollection_, dcchs) ) {

    int neec = 0;

    for ( EcalRawDataCollection::const_iterator dcchItr = dcchs->begin(); dcchItr != dcchs->end(); ++dcchItr ) {

      EcalDCCHeaderBlock dcch = (*dcchItr);

      if ( Numbers::subDet( dcch ) != EcalEndcap ) continue;

      neec++;

      map<int, EcalDCCHeaderBlock>::iterator i = dccMap.find(dcch.id());
      if ( i != dccMap.end() ) continue;

      dccMap[dcch.id()] = dcch;

      if ( meEEDCC_ ) meEEDCC_->Fill(Numbers::iSM( dcch, EcalEndcap )+0.5);

      if ( ! fixedRunNumber_ ) {
        runNumber_ = dcch.getRunNumber();
      }

      evtNumber_ = dcch.getLV1();

      if ( ! fixedRunType_ ) {
        runType_ = dcch.getRunType();
        evtType_ = runType_;
      }

    }

    LogDebug("EcalEndcapMonitor") << "event: " << ievt_ << " DCC headers collection size: " << neec;

  } else {

    LogWarning("EcalEndcapMonitorModule") << EcalRawDataCollection_ << " not available";

  }

  if ( meRunType_ ) meRunType_->Fill(runType_);

  if ( evtType_ >= 0 && evtType_ <= 22 ) {
    if ( meEvtType_ ) meEvtType_->Fill(evtType_+0.5);
  } else {
    LogWarning("EcalEndcapMonitor") << "Unknown event type = " << evtType_ << " at event: " << ievt_;
    if ( meEvtType_ ) meEvtType_->Fill( -1 );
  }

  if ( ievt_ == 1 ) {
    LogInfo("EcalEndcapMonitor") << "processing run " << runNumber_;
    // begin-of-run
    if ( meStatus_ ) meStatus_->Fill(0);
  } else {
    // running
    if ( meStatus_ ) meStatus_->Fill(1);
  }

  if ( meRun_ ) meRun_->Fill(runNumber_);
  if ( meEvt_ ) meEvt_->Fill(evtNumber_);

  // this should give enough time to all the MEs to reach the Collector,
  // and then hopefully the Client, even for short runs

  if ( ievt_ == 1 ) {
    if ( enableMonitorDaemon_ ) sleep(5);
  }

  // pause the shipping of monitoring elements
  dbe_->lock();

  Handle<EEDigiCollection> digis;

  if ( e.getByLabel(EEDigiCollection_, digis) ) {

    int need = digis->size();
    LogDebug("EcalEndcapMonitor") << "event " << ievt_ << " digi collection size " << need;

    int counter[18] = { 0 };

    if ( meEEdigis_[0] ) meEEdigis_[0]->Fill(float(need));

    for ( EEDigiCollection::const_iterator digiItr = digis->begin(); digiItr != digis->end(); ++digiItr ) {

      EEDataFrame dataframe = (*digiItr);
      EEDetId id = dataframe.id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      counter[ism-1]++;

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      LogDebug("EcalEndcapMonitor") << " det id = " << id;
      LogDebug("EcalEndcapMonitor") << " sm, ix, iy " << ism << " " << ix << " " << iy;

    }

    for (int i = 0; i < 18; i++) {

      if ( meEEdigis_[1] ) meEEdigis_[1]->Fill(i+1+0.5, counter[i]);

    }

  } else {

    LogWarning("EcalEndcapMonitorModule") << EEDigiCollection_ << " not available";

  }

  Handle<EcalUncalibratedRecHitCollection> hits;

  if ( e.getByLabel(EcalUncalibratedRecHitCollection_, hits) ) {

    int neeh = hits->size();
    LogDebug("EcalEndcapMonitor") << "event " << ievt_ << " hits collection size " << neeh;

    if ( meEEhits_[0] ) meEEhits_[0]->Fill(float(neeh));

    int counter[18] = { 0 };

    for ( EcalUncalibratedRecHitCollection::const_iterator hitItr = hits->begin(); hitItr != hits->end(); ++hitItr ) {

      EcalUncalibratedRecHit hit = (*hitItr);
      EEDetId id = hit.id();

      int ix = id.ix();
      int iy = id.iy();

      int ism = Numbers::iSM( id );

      counter[ism-1]++;

      if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;

      float xix = ix - 0.5;
      float xiy = iy - 0.5;

      LogDebug("EcalEndcapMonitor") << " det id = " << id;
      LogDebug("EcalEndcapMonitor") << " sm, ix, iy " << ism << " " << ix << " " << iy;

      float xval = hit.amplitude();

      LogDebug("EcalEndcapMonitor") << " hit amplitude " << xval;

      if ( enableEventDisplay_ ) {

        if ( xval >= 10 ) {
          if ( meEvent_[ism-1] ) meEvent_[ism-1]->Fill(xix, xiy, xval);
        }

      }

    }

    for (int i = 0; i < 18; i++) {

      if ( meEEhits_[1] ) meEEhits_[1]->Fill(i+1+0.5, counter[i]);

    }

  } else {

    LogWarning("EcalEndcapMonitorModule") << EcalUncalibratedRecHitCollection_ << " not available";

  }

  Handle<EcalTrigPrimDigiCollection> tpdigis;

  if ( e.getByLabel(EcalTrigPrimDigiCollection_, tpdigis) ) {

    int neetpd = 0;
    int counter[18] = { 0 };

    for ( EcalTrigPrimDigiCollection::const_iterator tpdigiItr = tpdigis->begin();
          tpdigiItr != tpdigis->end(); ++tpdigiItr ) {

      EcalTriggerPrimitiveDigi data = (*tpdigiItr);
      EcalTrigTowerDetId idt = data.id();

      if ( Numbers::subDet( idt ) != EcalEndcap ) continue;

      int ismt = Numbers::iSM( idt );

      neetpd++;
      counter[ismt-1]++;

    }

    LogDebug("EcalEndcapMonitor") << "event " << ievt_ << " TP digi collection size " << neetpd;
    if ( meEEtpdigis_[0] ) meEEtpdigis_[0]->Fill(float(neetpd));

    for (int i = 0; i < 18; i++) {

      if ( meEEtpdigis_[1] ) meEEtpdigis_[1]->Fill(i+1+0.5, counter[i]);

    }

  } else {

    LogWarning("EcalEndcapMonitorModule") << EcalTrigPrimDigiCollection_ << " not available";

  }

  // resume the shipping of monitoring elements
  dbe_->unlock();

}

