/*
 * \file EETimingClient.cc
 *
 * $Date: 2010/04/14 16:13:40 $
 * $Revision: 1.99 $
 * \author G. Della Ricca
 *
*/

#include <memory>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMStore.h"

#ifdef WITH_ECAL_COND_DB
#include "OnlineDB/EcalCondDB/interface/MonTimingCrystalDat.h"
#include "OnlineDB/EcalCondDB/interface/RunCrystalErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/RunTTErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/EcalCondDBInterface.h"
#endif

#include "CondTools/Ecal/interface/EcalErrorDictionary.h"

#include "DQM/EcalCommon/interface/EcalErrorMask.h"
#include "DQM/EcalCommon/interface/UtilsClient.h"
#include "DQM/EcalCommon/interface/LogicID.h"
#include "DQM/EcalCommon/interface/Numbers.h"

#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include <DQM/EcalEndcapMonitorClient/interface/EETimingClient.h>

EETimingClient::EETimingClient(const edm::ParameterSet& ps) {

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  // verbose switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", true);

  // debug switch
  debug_ = ps.getUntrackedParameter<bool>("debug", false);

  // prefixME path
  prefixME_ = ps.getUntrackedParameter<std::string>("prefixME", "");

  // enableCleanup_ switch
  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  // vector of selected Super Modules (Defaults to all 18).
  superModules_.reserve(18);
  for ( unsigned int i = 1; i <= 18; i++ ) superModules_.push_back(i);
  superModules_ = ps.getUntrackedParameter<std::vector<int> >("superModules", superModules_);

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    h01_[ism-1] = 0;
    h02_[ism-1] = 0;

    meh01_[ism-1] = 0;
    meh02_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    meg01_[ism-1] = 0;

    mea01_[ism-1] = 0;

    mep01_[ism-1] = 0;

    mer01_[ism-1] = 0;

  }

  expectedMean_ = 0.0;
  discrepancyMean_ = 25.0;
  RMSThreshold_ = 125.;

}

EETimingClient::~EETimingClient() {

}

void EETimingClient::beginJob(void) {

  dqmStore_ = edm::Service<DQMStore>().operator->();

  if ( debug_ ) std::cout << "EETimingClient: beginJob" <<  std::endl;

  ievt_ = 0;
  jevt_ = 0;

}

void EETimingClient::beginRun(void) {

  if ( debug_ ) std::cout << "EETimingClient: beginRun" <<  std::endl;

  jevt_ = 0;

  this->setup();

}

void EETimingClient::endJob(void) {

  if ( debug_ ) std::cout << "EETimingClient: endJob, ievt = " << ievt_ <<  std::endl;

  this->cleanup();

}

void EETimingClient::endRun(void) {

  if ( debug_ ) std::cout << "EETimingClient: endRun, jevt = " << jevt_ <<  std::endl;

  this->cleanup();

}

void EETimingClient::setup(void) {

  char histo[200];

  dqmStore_->setCurrentFolder( prefixME_ + "/EETimingClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg01_[ism-1] ) dqmStore_->removeElement( meg01_[ism-1]->getName() );
    sprintf(histo, "EETMT timing quality %s", Numbers::sEE(ism).c_str());
    meg01_[ism-1] = dqmStore_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    meg01_[ism-1]->setAxisTitle("jx", 1);
    meg01_[ism-1]->setAxisTitle("jy", 2);

    if ( mea01_[ism-1] ) dqmStore_->removeElement( mea01_[ism-1]->getName() );
    sprintf(histo, "EETMT timing %s", Numbers::sEE(ism).c_str());
    mea01_[ism-1] = dqmStore_->book1D(histo, histo, 850, 0., 850.);
    mea01_[ism-1]->setAxisTitle("channel", 1);
    mea01_[ism-1]->setAxisTitle("time (ns)", 2);

    if ( mep01_[ism-1] ) dqmStore_->removeElement( mep01_[ism-1]->getName() );
    sprintf(histo, "EETMT timing mean %s", Numbers::sEE(ism).c_str());
    mep01_[ism-1] = dqmStore_->book1D(histo, histo, 100, -50., 50.);
    mep01_[ism-1]->setAxisTitle("mean (ns)", 1);

    if ( mer01_[ism-1] ) dqmStore_->removeElement( mer01_[ism-1]->getName() );
    sprintf(histo, "EETMT timing rms %s", Numbers::sEE(ism).c_str());
    mer01_[ism-1] = dqmStore_->book1D(histo, histo, 100, 0.0, 150.0);
    mer01_[ism-1]->setAxisTitle("rms (ns)", 1);

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg01_[ism-1] ) meg01_[ism-1]->Reset();

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, 6. );

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, 2. );
        }

      }
    }

    if ( mea01_[ism-1] ) mea01_[ism-1]->Reset();
    if ( mep01_[ism-1] ) mep01_[ism-1]->Reset();
    if ( mer01_[ism-1] ) mer01_[ism-1]->Reset();

  }

}

void EETimingClient::cleanup(void) {

  if ( ! enableCleanup_ ) return;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( cloneME_ ) {
      if ( h01_[ism-1] ) delete h01_[ism-1];
      if ( h02_[ism-1] ) delete h02_[ism-1];
    }

    h01_[ism-1] = 0;
    h02_[ism-1] = 0;

    meh01_[ism-1] = 0;
    meh02_[ism-1] = 0;

  }

  dqmStore_->setCurrentFolder( prefixME_ + "/EETimingClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg01_[ism-1] ) dqmStore_->removeElement( meg01_[ism-1]->getName() );
    meg01_[ism-1] = 0;

    if ( mea01_[ism-1] ) dqmStore_->removeElement( mea01_[ism-1]->getName() );
    mea01_[ism-1] = 0;

    if ( mep01_[ism-1] ) dqmStore_->removeElement( mep01_[ism-1]->getName() );
    mep01_[ism-1] = 0;

    if ( mer01_[ism-1] ) dqmStore_->removeElement( mer01_[ism-1]->getName() );
    mer01_[ism-1] = 0;

  }

}

#ifdef WITH_ECAL_COND_DB
bool EETimingClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov, bool& status) {

  status = true;

  EcalLogicID ecid;

  MonTimingCrystalDat t;
  map<EcalLogicID, MonTimingCrystalDat> dataset;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( verbose_ ) {
      std::cout << " " << Numbers::sEE(ism) << " (ism=" << ism << ")" <<  std::endl;
      std::cout <<  std::endl;
      UtilsClient::printBadChannels(meg01_[ism-1], h01_[ism-1]);
    }

    float num01;
    float mean01;
    float rms01;

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( ! Numbers::validEE(ism, jx, jy) ) continue;

        bool update01;

        update01 = UtilsClient::getBinStatistics(h01_[ism-1], ix, iy, num01, mean01, rms01);
        // Task timing map is shifted of +50 ns for graphical reasons. Shift back it.
        mean01 -= 50.;

        if ( update01 ) {

          if ( Numbers::icEE(ism, jx, jy) == 1 ) {

            if ( verbose_ ) {
              std::cout << "Preparing dataset for " << Numbers::sEE(ism) << " (ism=" << ism << ")" <<  std::endl;
              std::cout << "crystal (" << Numbers::ix0EE(i+1)+ix << "," << Numbers::iy0EE(i+1)+iy << ") " << num01  << " " << mean01 << " " << rms01  <<  std::endl;
              std::cout <<  std::endl;
            }

          }

          t.setTimingMean(mean01);
          t.setTimingRMS(rms01);

          if ( UtilsClient::getBinStatus(meg01_[ism-1], ix, iy) ) {
            t.setTaskStatus(true);
          } else {
            t.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQuality(meg01_[ism-1], ix, iy);

          int ic = Numbers::indexEE(ism, jx, jy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            ecid = LogicID::getEcalLogicID("EE_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
            dataset[ecid] = t;
          }

        }

      }
    }

  }

  if ( econn ) {
    try {
      if ( verbose_ ) std::cout << "Inserting MonTimingCrystalDat ..." <<  std::endl;
      if ( dataset.size() != 0 ) econn->insertDataArraySet(&dataset, moniov);
      if ( verbose_ ) std::cout << "done." <<  std::endl;
    } catch (runtime_error &e) {
      cerr << e.what() <<  std::endl;
    }
  }

  return true;

}
#endif

void EETimingClient::analyze(void) {

  ievt_++;
  jevt_++;
  if ( ievt_ % 10 == 0 ) {
    if ( debug_ ) std::cout << "EETimingClient: ievt/jevt = " << ievt_ << "/" << jevt_ <<  std::endl;
  }

  uint64_t bits01 = 0;
  bits01 |= EcalErrorDictionary::getMask("PHYSICS_MEAN_TIMING_WARNING");
  bits01 |= EcalErrorDictionary::getMask("PHYSICS_RMS_TIMING_WARNING");

  char histo[200];

  MonitorElement* me;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    sprintf(histo, (prefixME_ + "/EETimingTask/EETMT timing %s").c_str(), Numbers::sEE(ism).c_str());
    me = dqmStore_->get(histo);
    h01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h01_[ism-1] );
    meh01_[ism-1] = me;

    sprintf(histo, (prefixME_ + "/EETimingTask/EETMT timing vs amplitude %s").c_str(), Numbers::sEE(ism).c_str());
    me = dqmStore_->get(histo);
    h02_[ism-1] = UtilsClient::getHisto<TH2F*>( me, cloneME_, h02_[ism-1] );
    meh02_[ism-1] = me;

    if ( meg01_[ism-1] ) meg01_[ism-1]->Reset();
    if ( mea01_[ism-1] ) mea01_[ism-1]->Reset();
    if ( mep01_[ism-1] ) mep01_[ism-1]->Reset();
    if ( mer01_[ism-1] ) mer01_[ism-1]->Reset();

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent(ix, iy, 6.);

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, 2. );
        }

        bool update01;

        float num01;
        float mean01;
        float rms01;

        update01 = UtilsClient::getBinStatistics(h01_[ism-1], ix, iy, num01, mean01, rms01);
        // Task timing map is shifted of +50 ns for graphical reasons. Shift back it.
        mean01 -= 50.;

        if ( update01 ) {

          float val;

          val = 1.;
          if ( fabs(mean01 - expectedMean_) > discrepancyMean_ )
            val = 0.;
          if ( rms01 > RMSThreshold_ )
            val = 0.;
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent(ix, iy, val);

          int ic = Numbers::icEE(ism, jx, jy);

          if ( ic != -1 ) {
            if ( mea01_[ism-1] ) {
              mea01_[ism-1]->setBinContent(ic, mean01);
              mea01_[ism-1]->setBinError(ic, rms01);
            }

            if ( mep01_[ism-1] ) mep01_[ism-1]->Fill(mean01);
            if ( mer01_[ism-1] ) mer01_[ism-1]->Fill(rms01);
          }

        }

      }
    }

  }

#ifdef WITH_ECAL_COND_DB
  if ( EcalErrorMask::mapCrystalErrors_.size() != 0 ) {
    map<EcalLogicID, RunCrystalErrorsDat>::const_iterator m;
    for (m = EcalErrorMask::mapCrystalErrors_.begin(); m != EcalErrorMask::mapCrystalErrors_.end(); m++) {

      if ( (m->second).getErrorBits() & bits01 ) {
        EcalLogicID ecid = m->first;

        if ( strcmp(ecid.getMapsTo().c_str(), "EE_crystal_number") != 0 ) continue;

        int iz = ecid.getID1();
        int ix = ecid.getID2();
        int iy = ecid.getID3();

        for ( unsigned int i=0; i<superModules_.size(); i++ ) {
          int ism = superModules_[i];
          std::vector<int>::iterator iter = find(superModules_.begin(), superModules_.end(), ism);
          if (iter == superModules_.end()) continue;
          if ( iz == -1 && ( ism >=  1 && ism <=  9 ) ) {
            int jx = 101 - ix - Numbers::ix0EE(ism);
            int jy = iy - Numbers::iy0EE(ism);
            if ( Numbers::validEE(ism, ix, iy) ) UtilsClient::maskBinContent( meg01_[ism-1], jx, jy );
          }
          if ( iz == +1 && ( ism >= 10 && ism <= 18 ) ) {
            int jx = ix - Numbers::ix0EE(ism);
            int jy = iy - Numbers::iy0EE(ism);
            if ( Numbers::validEE(ism, ix, iy) ) UtilsClient::maskBinContent( meg01_[ism-1], jx, jy );
          }
        }

      }

    }
  }

  if ( EcalErrorMask::mapTTErrors_.size() != 0 ) {
    map<EcalLogicID, RunTTErrorsDat>::const_iterator m;
    for (m = EcalErrorMask::mapTTErrors_.begin(); m != EcalErrorMask::mapTTErrors_.end(); m++) {

      if ( (m->second).getErrorBits() & bits01 ) {
        EcalLogicID ecid = m->first;

        if ( strcmp(ecid.getMapsTo().c_str(), "EE_readout_tower") != 0 ) continue;

        int idcc = ecid.getID1() - 600;

        int ism = -1;
        if ( idcc >=   1 && idcc <=   9 ) ism = idcc;
        if ( idcc >=  46 && idcc <=  54 ) ism = idcc - 45 + 9;
        std::vector<int>::iterator iter = find(superModules_.begin(), superModules_.end(), ism);
        if (iter == superModules_.end()) continue;

        int itt = ecid.getID2();

        if ( itt > 70 ) continue;

        if ( itt >= 42 && itt <= 68 ) continue;

        if ( ( ism == 8 || ism == 17 ) && ( itt >= 18 && itt <= 24 ) ) continue;

        if ( itt >= 1 && itt <= 68 ) {
          std::vector<DetId>* crystals = Numbers::crystals( idcc, itt );
          for ( unsigned int i=0; i<crystals->size(); i++ ) {
            EEDetId id = (*crystals)[i];
            int ix = id.ix();
            int iy = id.iy();
            if ( ism >= 1 && ism <= 9 ) ix = 101 - ix;
            int jx = ix - Numbers::ix0EE(ism);
            int jy = iy - Numbers::iy0EE(ism);
            UtilsClient::maskBinContent( meg01_[ism-1], jx, jy );
          }
        }

      }

    }
  }
#endif

}

