/*
 * \file EEOccupancyTask.cc
 *
 * $Date: 2007/04/02 16:23:13 $
 * $Revision: 1.1 $
 * \author G. Della Ricca
 * \author G. Franzoni
 *
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/Daemon/interface/MonitorDaemon.h"

#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include <DQM/EcalEndcapMonitorTasks/interface/EEOccupancyTask.h>

using namespace cms;
using namespace edm;
using namespace std;

EEOccupancyTask::EEOccupancyTask(const ParameterSet& ps){

  init_ = false;

  EBDigiCollection_ = ps.getParameter<edm::InputTag>("EBDigiCollection");
  EcalPnDiodeDigiCollection_ = ps.getParameter<edm::InputTag>("EcalPnDiodeDigiCollection");

  for (int i = 0; i < 36; i++) {
    meOccupancy_[i]    = 0;
    meOccupancyMem_[i] = 0;
  }

  // get hold of back-end interface
  dbe_ = Service<DaqMonitorBEInterface>().operator->();

}

EEOccupancyTask::~EEOccupancyTask(){

}

void EEOccupancyTask::beginJob(const EventSetup& c){

  ievt_ = 0;

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEOccupancyTask");
    dbe_->rmdir("EcalEndcap/EEOccupancyTask");
  }

}

void EEOccupancyTask::setup(void){

  init_ = true;

  Char_t histo[200];

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEOccupancyTask");

    for (int i = 0; i < 36; i++) {
      sprintf(histo, "EEOT occupancy SM%02d", i+1);
      meOccupancy_[i] = dbe_->book2D(histo, histo, 85, 0., 85., 20, 0., 20.);
      dbe_->tag(meOccupancy_[i], i+1);
    }
    for (int i = 0; i < 36; i++) {
      sprintf(histo, "EEOT MEM occupancy SM%02d", i+1);
      meOccupancyMem_[i] = dbe_->book2D(histo, histo, 10, 0., 10., 5, 0., 5.);
      dbe_->tag(meOccupancyMem_[i], i+1);
    }

  }

}

void EEOccupancyTask::cleanup(void){

  if ( dbe_ ) {
    dbe_->setCurrentFolder("EcalEndcap/EEOccupancyTask");

    for (int i = 0; i < 36; i++) {
      if ( meOccupancy_[i] ) dbe_->removeElement( meOccupancy_[i]->getName() );
      meOccupancy_[i] = 0;
      if ( meOccupancyMem_[i] ) dbe_->removeElement( meOccupancyMem_[i]->getName() );
      meOccupancyMem_[i] = 0;
    }

  }

  init_ = false;

}

void EEOccupancyTask::endJob(void) {

  LogInfo("EEOccupancyTask") << "analyzed " << ievt_ << " events";

  if ( init_ ) this->cleanup();

}

void EEOccupancyTask::analyze(const Event& e, const EventSetup& c){

  if ( ! init_ ) this->setup();

  ievt_++;

  try {

    Handle<EBDigiCollection> digis;
    e.getByLabel(EBDigiCollection_, digis);

    int nebd = digis->size();
    LogDebug("EEOccupancyTask") << "event " << ievt_ << " digi collection size " << nebd;

    for ( EBDigiCollection::const_iterator digiItr = digis->begin(); digiItr != digis->end(); ++digiItr ) {

      EBDataFrame dataframe = (*digiItr);
      EBDetId id = dataframe.id();

      int ic = id.ic();
      int ie = (ic-1)/20 + 1;
      int ip = (ic-1)%20 + 1;

      int ism = id.ism();

      float xie = ie - 0.5;
      float xip = ip - 0.5;

      LogDebug("EEOccupancyTask") << " det id = " << id;
      LogDebug("EEOccupancyTask") << " sm, eta, phi " << ism << " " << ie << " " << ip;

      if ( xie <= 0. || xie >= 85. || xip <= 0. || xip >= 20. ) {
        LogWarning("EEOccupancyTask") << " det id = " << id;
        LogWarning("EEOccupancyTask") << " sm, eta, phi " << ism << " " << ie << " " << ip;
        LogWarning("EEOccupancyTask") << " xie, xip " << xie << " " << xip;
      }

      if ( meOccupancy_[ism-1] ) meOccupancy_[ism-1]->Fill(xie, xip);

    }

  } catch ( exception& ex) {

    LogWarning("EEOccupancyTask") << EBDigiCollection_ << " not available";

  }

  try {

    Handle<EcalPnDiodeDigiCollection> PNs;
    e.getByLabel(EcalPnDiodeDigiCollection_, PNs);

    // filling mem occupancy only for the 5 channels belonging
    // to a fully reconstructed PN's

    for ( EcalPnDiodeDigiCollection::const_iterator pnItr = PNs->begin(); pnItr != PNs->end(); ++pnItr ) {

      int   ism   = (*pnItr).id().iDCCId();
      float PnId  = (*pnItr).id().iPnId();
      PnId        = PnId - 0.5;
      float st    = 0.0;

      for (int chInStrip = 1; chInStrip <= 5; chInStrip++){
        if ( meOccupancyMem_[ism-1] ) {
           st = chInStrip - 0.5;
           meOccupancyMem_[ism-1]->Fill(PnId, st);
        }
      }

    }

  } catch ( exception& ex) {

    LogWarning("EEOccupancyTask") << EcalPnDiodeDigiCollection_ << " not available";

  }

}

