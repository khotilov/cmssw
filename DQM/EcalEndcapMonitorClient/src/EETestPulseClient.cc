/*
 * \file EETestPulseClient.cc
 *
 * $Date: 2007/11/10 14:56:39 $
 * $Revision: 1.45 $
 * \author G. Della Ricca
 * \author F. Cossutti
 *
*/

#include <memory>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/UI/interface/MonitorUIRoot.h"

#include "OnlineDB/EcalCondDB/interface/RunTag.h"
#include "OnlineDB/EcalCondDB/interface/RunIOV.h"
#include "OnlineDB/EcalCondDB/interface/MonTestPulseDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPulseShapeDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPNMGPADat.h"
#include "OnlineDB/EcalCondDB/interface/RunCrystalErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/RunPNErrorsDat.h"

#include "CondTools/Ecal/interface/EcalErrorDictionary.h"

#include "DQM/EcalCommon/interface/EcalErrorMask.h"
#include <DQM/EcalCommon/interface/UtilsClient.h>
#include <DQM/EcalCommon/interface/LogicID.h>
#include <DQM/EcalCommon/interface/Numbers.h>

#include <DQM/EcalEndcapMonitorClient/interface/EETestPulseClient.h>

using namespace cms;
using namespace edm;
using namespace std;

EETestPulseClient::EETestPulseClient(const ParameterSet& ps){

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  // MonitorDaemon switch
  enableMonitorDaemon_ = ps.getUntrackedParameter<bool>("enableMonitorDaemon", true);

  // prefix to ME paths
  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  // vector of selected Super Modules (Defaults to all 18).
  superModules_.reserve(18);
  for ( unsigned int i = 1; i <= 18; i++ ) superModules_.push_back(i);
  superModules_ = ps.getUntrackedParameter<vector<int> >("superModules", superModules_);

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    ha01_[ism-1] = 0;
    ha02_[ism-1] = 0;
    ha03_[ism-1] = 0;

    hs01_[ism-1] = 0;
    hs02_[ism-1] = 0;
    hs03_[ism-1] = 0;

    i01_[ism-1] = 0;
    i02_[ism-1] = 0;
    i03_[ism-1] = 0;
    i04_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    meg01_[ism-1] = 0;
    meg02_[ism-1] = 0;
    meg03_[ism-1] = 0;

    meg04_[ism-1] = 0;
    meg05_[ism-1] = 0;

    mea01_[ism-1] = 0;
    mea02_[ism-1] = 0;
    mea03_[ism-1] = 0;

    mer04_[ism-1] = 0;
    mer05_[ism-1] = 0;

    mes01_[ism-1] = 0;
    mes02_[ism-1] = 0;
    mes03_[ism-1] = 0;

  }

  percentVariation_ = 0.2;
  RMSThreshold_ = 300.0;

  amplitudeThresholdPnG01_ = 200./16.;
  amplitudeThresholdPnG16_ = 200.;

  pedPnExpectedMean_[0] = 750.0;
  pedPnExpectedMean_[1] = 750.0;

  pedPnDiscrepancyMean_[0] = 100.0;
  pedPnDiscrepancyMean_[1] = 100.0;

  pedPnRMSThreshold_[0] = 1.0;
  pedPnRMSThreshold_[1] = 3.0;

}

EETestPulseClient::~EETestPulseClient(){

}

void EETestPulseClient::beginJob(MonitorUserInterface* mui){

  mui_ = mui;
  dbe_ = mui->getBEInterface();

  if ( verbose_ ) cout << "EETestPulseClient: beginJob" << endl;

  ievt_ = 0;
  jevt_ = 0;

}

void EETestPulseClient::beginRun(void){

  if ( verbose_ ) cout << "EETestPulseClient: beginRun" << endl;

  jevt_ = 0;

  this->setup();

  this->subscribe();

}

void EETestPulseClient::endJob(void) {

  if ( verbose_ ) cout << "EETestPulseClient: endJob, ievt = " << ievt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EETestPulseClient::endRun(void) {

  if ( verbose_ ) cout << "EETestPulseClient: endRun, jevt = " << jevt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EETestPulseClient::setup(void) {

  Char_t histo[200];

  dbe_->setCurrentFolder( "EcalEndcap/EETestPulseClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg01_[ism-1] ) dbe_->removeElement( meg01_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse quality G01 %s", Numbers::sEE(ism).c_str());
    meg01_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    meg01_[ism-1]->setAxisTitle("ix", 1);
    meg01_[ism-1]->setAxisTitle("iy", 2);
    if ( meg02_[ism-1] ) dbe_->removeElement( meg02_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse quality G06 %s", Numbers::sEE(ism).c_str());
    meg02_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    meg02_[ism-1]->setAxisTitle("ix", 1);
    meg02_[ism-1]->setAxisTitle("iy", 2);
    if ( meg03_[ism-1] ) dbe_->removeElement( meg03_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse quality G12 %s", Numbers::sEE(ism).c_str());
    meg03_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    meg03_[ism-1]->setAxisTitle("ix", 1);
    meg03_[ism-1]->setAxisTitle("iy", 2);

    if ( meg04_[ism-1] ) dbe_->removeElement( meg04_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse quality PNs G01 %s", Numbers::sEE(ism).c_str());
    meg04_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    meg04_[ism-1]->setAxisTitle("pseudo-strip", 1);
    meg04_[ism-1]->setAxisTitle("channel", 2);
    if ( meg05_[ism-1] ) dbe_->removeElement( meg05_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse quality PNs G16 %s", Numbers::sEE(ism).c_str());
    meg05_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    meg05_[ism-1]->setAxisTitle("pseudo-strip", 1);
    meg05_[ism-1]->setAxisTitle("channel", 2);

    if ( mea01_[ism-1] ) dbe_->removeElement( mea01_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse amplitude G01 %s", Numbers::sEE(ism).c_str());
    mea01_[ism-1] = dbe_->book1D(histo, histo, 850, 0., 850.);
    mea01_[ism-1]->setAxisTitle("channel", 1);
    mea01_[ism-1]->setAxisTitle("amplitude", 2);
    if ( mea02_[ism-1] ) dbe_->removeElement( mea02_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse amplitude G06 %s", Numbers::sEE(ism).c_str());
    mea02_[ism-1] = dbe_->book1D(histo, histo, 850, 0., 850.);
    mea02_[ism-1]->setAxisTitle("channel", 1);
    mea02_[ism-1]->setAxisTitle("amplitude", 2);
    if ( mea03_[ism-1] ) dbe_->removeElement( mea03_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse amplitude G12 %s", Numbers::sEE(ism).c_str());
    mea03_[ism-1] = dbe_->book1D(histo, histo, 850, 0., 850.);
    mea03_[ism-1]->setAxisTitle("channel", 1);
    mea03_[ism-1]->setAxisTitle("amplitude", 2);

    if ( mer04_[ism-1] ) dbe_->removeElement( mer04_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G01", Numbers::sEE(ism).c_str());
    mer04_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    mer04_[ism-1]->setAxisTitle("rms", 1);
    if ( mer05_[ism-1] ) dbe_->removeElement( mer05_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G16", Numbers::sEE(ism).c_str());
    mer05_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    mer05_[ism-1]->setAxisTitle("rms", 1);

    if ( mes01_[ism-1] ) dbe_->removeElement( mes01_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse shape G01 %s", Numbers::sEE(ism).c_str());
    mes01_[ism-1] = dbe_->book1D(histo, histo, 10, 0., 10.);
    mes01_[ism-1]->setAxisTitle("sample", 1);
    mes01_[ism-1]->setAxisTitle("amplitude", 2);
    if ( mes02_[ism-1] ) dbe_->removeElement( mes02_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse shape G06 %s", Numbers::sEE(ism).c_str());
    mes02_[ism-1] = dbe_->book1D(histo, histo, 10, 0., 10.);
    mes02_[ism-1]->setAxisTitle("sample", 1);
    mes02_[ism-1]->setAxisTitle("amplitude", 2);
    if ( mes03_[ism-1] ) dbe_->removeElement( mes03_[ism-1]->getName() );
    sprintf(histo, "EETPT test pulse shape G12 %s", Numbers::sEE(ism).c_str());
    mes03_[ism-1] = dbe_->book1D(histo, histo, 10, 0., 10.);
    mes03_[ism-1]->setAxisTitle("sample", 1);
    mes03_[ism-1]->setAxisTitle("amplitude", 2);

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    meg01_[ism-1]->Reset();
    meg02_[ism-1]->Reset();
    meg03_[ism-1]->Reset();

    meg04_[ism-1]->Reset();
    meg05_[ism-1]->Reset();

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        meg01_[ism-1]->setBinContent( ix, iy, -1. );
        meg02_[ism-1]->setBinContent( ix, iy, -1. );
        meg03_[ism-1]->setBinContent( ix, iy, -1. );

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          meg01_[ism-1]->setBinContent( ix, iy, 2. );
          meg02_[ism-1]->setBinContent( ix, iy, 2. );
          meg03_[ism-1]->setBinContent( ix, iy, 2. );
        }

      }
    }

    for ( int i = 1; i <= 10; i++ ) {

        meg04_[ism-1]->setBinContent( i, 1, 2. );
        meg05_[ism-1]->setBinContent( i, 1, 2. );

    }

    mea01_[ism-1]->Reset();
    mea02_[ism-1]->Reset();
    mea03_[ism-1]->Reset();

    mer04_[ism-1]->Reset();
    mer05_[ism-1]->Reset();

    mes01_[ism-1]->Reset();
    mes02_[ism-1]->Reset();
    mes03_[ism-1]->Reset();

  }

}

void EETestPulseClient::cleanup(void) {

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( cloneME_ ) {
      if ( ha01_[ism-1] ) delete ha01_[ism-1];
      if ( ha02_[ism-1] ) delete ha02_[ism-1];
      if ( ha03_[ism-1] ) delete ha03_[ism-1];

      if ( hs01_[ism-1] ) delete hs01_[ism-1];
      if ( hs02_[ism-1] ) delete hs02_[ism-1];
      if ( hs03_[ism-1] ) delete hs03_[ism-1];

      if ( i01_[ism-1] ) delete i01_[ism-1];
      if ( i02_[ism-1] ) delete i02_[ism-1];
      if ( i03_[ism-1] ) delete i03_[ism-1];
      if ( i04_[ism-1] ) delete i04_[ism-1];
    }

    ha01_[ism-1] = 0;
    ha02_[ism-1] = 0;
    ha03_[ism-1] = 0;

    hs01_[ism-1] = 0;
    hs02_[ism-1] = 0;
    hs03_[ism-1] = 0;

    i01_[ism-1] = 0;
    i02_[ism-1] = 0;
    i03_[ism-1] = 0;
    i04_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    dbe_->setCurrentFolder( "EcalEndcap/EETestPulseClient" );

    if ( meg01_[ism-1] ) dbe_->removeElement( meg01_[ism-1]->getName() );
    meg01_[ism-1] = 0;
    if ( meg02_[ism-1] ) dbe_->removeElement( meg02_[ism-1]->getName() );
    meg02_[ism-1] = 0;
    if ( meg03_[ism-1] ) dbe_->removeElement( meg03_[ism-1]->getName() );
    meg03_[ism-1] = 0;

    if ( meg04_[ism-1] ) dbe_->removeElement( meg04_[ism-1]->getName() );
    meg04_[ism-1] = 0;
    if ( meg05_[ism-1] ) dbe_->removeElement( meg05_[ism-1]->getName() );
    meg05_[ism-1] = 0;

    if ( mea01_[ism-1] ) dbe_->removeElement( mea01_[ism-1]->getName() );
    mea01_[ism-1] = 0;
    if ( mea02_[ism-1] ) dbe_->removeElement( mea02_[ism-1]->getName() );
    mea02_[ism-1] = 0;
    if ( mea03_[ism-1] ) dbe_->removeElement( mea03_[ism-1]->getName() );
    mea03_[ism-1] = 0;

    if ( mer04_[ism-1] ) dbe_->removeElement( mer04_[ism-1]->getName() );
    mer04_[ism-1] = 0;
    if ( mer05_[ism-1] ) dbe_->removeElement( mer05_[ism-1]->getName() );
    mer05_[ism-1] = 0;

    if ( mes01_[ism-1] ) dbe_->removeElement( mes01_[ism-1]->getName() );
    mes01_[ism-1] = 0;
    if ( mes02_[ism-1] ) dbe_->removeElement( mes02_[ism-1]->getName() );
    mes02_[ism-1] = 0;
    if ( mes03_[ism-1] ) dbe_->removeElement( mes03_[ism-1]->getName() );
    mes03_[ism-1] = 0;

  }

}

bool EETestPulseClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov) {

  bool status = true;

  EcalLogicID ecid;

  MonTestPulseDat adc;
  map<EcalLogicID, MonTestPulseDat> dataset1;
  MonPulseShapeDat shape;
  map<EcalLogicID, MonPulseShapeDat> dataset2;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    cout << " SM=" << ism << endl;
    cout << endl;

    UtilsClient::printBadChannels(meg01_[ism-1], ha01_[ism-1]);
    UtilsClient::printBadChannels(meg02_[ism-1], ha02_[ism-1]);
    UtilsClient::printBadChannels(meg03_[ism-1], ha03_[ism-1]);

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( ! Numbers::validEE(ism, jx, jy) ) continue;

        bool update01;
        bool update02;
        bool update03;

        float num01, num02, num03;
        float mean01, mean02, mean03;
        float rms01, rms02, rms03;

        update01 = UtilsClient::getBinStats(ha01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(ha02_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(ha03_[ism-1], ix, iy, num03, mean03, rms03);

        if ( update01 || update02 || update03 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;
            cout << "G01 (" << ix << "," << iy << ") " << num01 << " " << mean01 << " " << rms01 << endl;
            cout << "G06 (" << ix << "," << iy << ") " << num02 << " " << mean02 << " " << rms02 << endl;
            cout << "G12 (" << ix << "," << iy << ") " << num03 << " " << mean03 << " " << rms03 << endl;

            cout << endl;

          }

          adc.setADCMeanG1(mean01);
          adc.setADCRMSG1(rms01);

          adc.setADCMeanG6(mean02);
          adc.setADCRMSG6(rms02);

          adc.setADCMeanG12(mean03);
          adc.setADCRMSG12(rms03);

          if ( meg01_[ism-1] && int(meg01_[ism-1]->getBinContent( ix, iy )) % 3 == 1. &&
               meg02_[ism-1] && int(meg02_[ism-1]->getBinContent( ix, iy )) % 3 == 1. &&
               meg03_[ism-1] && int(meg03_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            adc.setTaskStatus(true);
          } else {
            adc.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg01_[ism-1], ix, iy) &&
                             UtilsClient::getBinQual(meg02_[ism-1], ix, iy) &&
                             UtilsClient::getBinQual(meg03_[ism-1], ix, iy);

          if ( ix == 1 && iy == 1 ) {

            vector<float> sample01, sample02, sample03;

            sample01.clear();
            sample02.clear();
            sample03.clear();

            if ( mes01_[ism-1] ) {
              for ( int i = 1; i <= 10; i++ ) {
                sample01.push_back(int(mes01_[ism-1]->getBinContent(i)));
              }
            } else {
              for ( int i = 1; i <= 10; i++ ) { sample01.push_back(-1.); }
            }

            if ( mes02_[ism-1] ) {
              for ( int i = 1; i <= 10; i++ ) {
                sample02.push_back(int(mes02_[ism-1]->getBinContent(i)));
              }
            } else {
              for ( int i = 1; i <= 10; i++ ) { sample02.push_back(-1.); }
            }

            if ( mes03_[ism-1] ) {
              for ( int i = 1; i <= 10; i++ ) {
                sample03.push_back(int(mes03_[ism-1]->getBinContent(i)));
              }
            } else {
              for ( int i = 1; i <= 10; i++ ) { sample03.push_back(-1.); }
            }

            cout << "sample01 = " << flush;
            for ( unsigned int i = 0; i < sample01.size(); i++ ) {
              cout << sample01[i] << " " << flush;
            }
            cout << endl;

            cout << "sample02 = " << flush;
            for ( unsigned int i = 0; i < sample02.size(); i++ ) {
              cout << sample02[i] << " " << flush;
            }
            cout << endl;

            cout << "sample03 = " << flush;
            for ( unsigned int i = 0; i < sample03.size(); i++ ) {
              cout << sample03[i] << " " << flush;
            }
            cout << endl;

            cout << endl;

            shape.setSamples(sample01,  1);
            shape.setSamples(sample02,  6);
            shape.setSamples(sample03, 12);

          }

          int ic = Numbers::indexEE(ism, ix, iy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EE_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
              dataset1[ecid] = adc;
              if ( ix == 1 && iy == 1 ) dataset2[ecid] = shape;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

      }
    }

  }

  if ( econn ) {
    try {
      cout << "Inserting MonTestPulseDat ... " << flush;
      if ( dataset1.size() != 0 ) econn->insertDataArraySet(&dataset1, moniov);
      if ( dataset2.size() != 0 ) econn->insertDataSet(&dataset2, moniov);
      cout << "done." << endl;
    } catch (runtime_error &e) {
      cerr << e.what() << endl;
    }
  }

  cout << endl;

  MonPNMGPADat pn;
  map<EcalLogicID, MonPNMGPADat> dataset3;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    cout << " SM=" << ism << endl;
    cout << endl;

    UtilsClient::printBadChannels(meg04_[ism-1], i01_[ism-1]);
    UtilsClient::printBadChannels(meg04_[ism-1], i03_[ism-1]);
    UtilsClient::printBadChannels(meg05_[ism-1], i02_[ism-1]);
    UtilsClient::printBadChannels(meg05_[ism-1], i04_[ism-1]);

    for ( int i = 1; i <= 10; i++ ) {

      bool update01;
      bool update02;
      bool update03;
      bool update04;

      float num01, num02, num03, num04;
      float mean01, mean02, mean03, mean04;
      float rms01, rms02, rms03, rms04;

      update01 = UtilsClient::getBinStats(i01_[ism-1], i, 0, num01, mean01, rms01);
      update02 = UtilsClient::getBinStats(i02_[ism-1], i, 0, num02, mean02, rms02);
      update03 = UtilsClient::getBinStats(i03_[ism-1], i, 0, num03, mean03, rms03);
      update04 = UtilsClient::getBinStats(i04_[ism-1], i, 0, num04, mean04, rms04);

      if ( update01 || update02 || update03 || update04 ) {

        if ( i == 1 ) {

          cout << "Preparing dataset for SM=" << ism << endl;

          cout << "PNs (" << i << ") G01 " << num01  << " " << mean01 << " " << rms01 << " " << num03 << " " << mean03 << " " << rms03 << endl;
          cout << "PNs (" << i << ") G16 " << num02  << " " << mean02 << " " << rms02 << " " << num04 << " " << mean04 << " " << rms04 << endl;

          cout << endl;

        }

        pn.setADCMeanG1(mean01);
        pn.setADCRMSG1(rms01);

        pn.setPedMeanG1(mean03);
        pn.setPedRMSG1(rms03);

        pn.setADCMeanG16(mean02);
        pn.setADCRMSG16(rms02);

        pn.setPedMeanG16(mean04);
        pn.setPedRMSG16(rms04);

        if ( meg04_[ism-1] && int(meg04_[ism-1]->getBinContent( i, 1 )) % 3 == 1. &&
             meg05_[ism-1] && int(meg05_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ) {
          pn.setTaskStatus(true);
        } else {
          pn.setTaskStatus(false);
        }

        status = status && UtilsClient::getBinQual(meg04_[ism-1], i, 1) &&
                           UtilsClient::getBinQual(meg05_[ism-1], i, 1);

        if ( econn ) {
          try {
            ecid = LogicID::getEcalLogicID("EE_LM_PN", Numbers::iSM(ism, EcalEndcap), i-1);
            dataset3[ecid] = pn;
          } catch (runtime_error &e) {
            cerr << e.what() << endl;
          }
        }

      }

    }

  }

  if ( econn ) {
    try {
      cout << "Inserting MonPNMGPADat ... " << flush;
      if ( dataset3.size() != 0 ) econn->insertDataArraySet(&dataset3, moniov);
      cout << "done." << endl;
    } catch (runtime_error &e) {
      cerr << e.what() << endl;
    }
  }

  return status;

}

void EETestPulseClient::subscribe(void){

  if ( verbose_ ) cout << "EETestPulseClient: subscribe" << endl;

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT shape %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT amplitude %s G06", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT shape %s G06", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT amplitude %s G12", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT shape %s G12", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs amplitude %s G16", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs pedestal %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs pedestal %s G16", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

  }

}

void EETestPulseClient::subscribeNew(void){

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT shape %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT amplitude %s G06", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT shape %s G06", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT amplitude %s G12", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT shape %s G12", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs amplitude %s G16", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs pedestal %s G01", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs pedestal %s G16", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

  }

}

void EETestPulseClient::unsubscribe(void){

  if ( verbose_ ) cout << "EETestPulseClient: unsubscribe" << endl;

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain01/EETPT shape %s G01", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT amplitude %s G06", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain06/EETPT shape %s G06", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT amplitude %s G12", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/Gain12/EETPT shape %s G12", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs amplitude %s G01", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs amplitude %s G16", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs pedestal %s G01", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs pedestal %s G16", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

  }

}

void EETestPulseClient::softReset(void){

}

void EETestPulseClient::analyze(void){

  ievt_++;
  jevt_++;
  if ( ievt_ % 10 == 0 ) {
    if ( verbose_ ) cout << "EETestPulseClient: ievt/jevt = " << ievt_ << "/" << jevt_ << endl;
  }

  uint64_t bits01 = 0;
  bits01 |= EcalErrorDictionary::getMask("TESTPULSE_LOW_GAIN_MEAN_WARNING");
  bits01 |= EcalErrorDictionary::getMask("TESTPULSE_LOW_GAIN_RMS_WARNING");

  uint64_t bits02 = 0;
  bits02 |= EcalErrorDictionary::getMask("TESTPULSE_MIDDLE_GAIN_MEAN_WARNING");
  bits02 |= EcalErrorDictionary::getMask("TESTPULSE_MIDDLE_GAIN_RMS_WARNING");

  uint64_t bits03 = 0;
  bits03 |= EcalErrorDictionary::getMask("TESTPULSE_HIGH_GAIN_MEAN_WARNING");
  bits03 |= EcalErrorDictionary::getMask("TESTPULSE_HIGH_GAIN_RMS_WARNING");

  uint64_t bits04 = 0;
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_MEAN_WARNING");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_RMS_WARNING");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_MEAN_ERROR");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_RMS_ERROR");

  uint64_t bits05 = 0;
  bits05 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_MEAN_WARNING");
  bits05 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_RMS_WARNING");
  bits05 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_MEAN_ERROR");
  bits05 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_RMS_ERROR");

  uint64_t bits06 = 0;
  bits06 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_MEAN_WARNING");
  bits06 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_RMS_WARNING");
  bits06 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_MEAN_ERROR");
  bits06 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_RMS_ERROR");

  map<EcalLogicID, RunCrystalErrorsDat> mask1;
  map<EcalLogicID, RunPNErrorsDat> mask2;

  EcalErrorMask::fetchDataSet(&mask1);
  EcalErrorMask::fetchDataSet(&mask2);

  Char_t histo[200];

  MonitorElement* me;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain01/EETPT amplitude %s G01").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    ha01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, ha01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain06/EETPT amplitude %s G06").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    ha02_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, ha02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain12/EETPT amplitude %s G12").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    ha03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, ha03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain01/EETPT shape %s G01").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain06/EETPT shape %s G06").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs02_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/Gain12/EETPT shape %s G12").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs amplitude %s G01").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i01_[ism-1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, i01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs amplitude %s G16").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i02_[ism-1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, i02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/PN/Gain01/EEPDT PNs pedestal %s G01").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i03_[ism-1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, i03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EETestPulseTask/PN/Gain16/EEPDT PNs pedestal %s G16").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i04_[ism-1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, i04_[ism-1] );

    meg01_[ism-1]->Reset();
    meg02_[ism-1]->Reset();
    meg03_[ism-1]->Reset();

    meg04_[ism-1]->Reset();
    meg05_[ism-1]->Reset();

    mea01_[ism-1]->Reset();
    mea02_[ism-1]->Reset();
    mea03_[ism-1]->Reset();

    mer04_[ism-1]->Reset();
    mer05_[ism-1]->Reset();

    mes01_[ism-1]->Reset(); 
    mes02_[ism-1]->Reset();
    mes03_[ism-1]->Reset();

    float meanAmpl01, meanAmpl02, meanAmpl03;

    int nCry01, nCry02, nCry03;

    meanAmpl01 = meanAmpl02 = meanAmpl03 = 0.;

    nCry01 = nCry02 = nCry03 = 0;

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        bool update01;
        bool update02;
        bool update03;

        float num01, num02, num03;
        float mean01, mean02, mean03;
        float rms01, rms02, rms03;

        update01 = UtilsClient::getBinStats(ha01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(ha02_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(ha03_[ism-1], ix, iy, num03, mean03, rms03);

        if ( update01 ) {
          meanAmpl01 += mean01;
          nCry01++;
        }

        if ( update02 ) {
          meanAmpl02 += mean02;
          nCry02++;
        }

        if ( update03 ) {
          meanAmpl03 += mean03;
          nCry03++;
        }

      }
    }

    if ( nCry01 > 0 ) meanAmpl01 /= float (nCry01);
    if ( nCry02 > 0 ) meanAmpl02 /= float (nCry02);
    if ( nCry03 > 0 ) meanAmpl03 /= float (nCry03);

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, -1. );
        if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, -1. );
        if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, -1. );

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, 2. );
          if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, 2. );
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, 2. );
        }

        float numEventsinCry[3] = {0., 0., 0.};

        if ( ha01_[ism-1] ) numEventsinCry[0] = ha01_[ism-1]->GetBinEntries(ha01_[ism-1]->GetBin(ix, iy));
        if ( ha02_[ism-1] ) numEventsinCry[1] = ha02_[ism-1]->GetBinEntries(ha02_[ism-1]->GetBin(ix, iy));
        if ( ha03_[ism-1] ) numEventsinCry[2] = ha03_[ism-1]->GetBinEntries(ha03_[ism-1]->GetBin(ix, iy));

        bool update01;
        bool update02;
        bool update03;

        float num01, num02, num03;
        float mean01, mean02, mean03;
        float rms01, rms02, rms03;

        update01 = UtilsClient::getBinStats(ha01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(ha02_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(ha03_[ism-1], ix, iy, num03, mean03, rms03);

        if ( update01 ) {

          float val;

          val = 1.;
          if ( fabs(mean01 - meanAmpl01) > fabs(percentVariation_ * meanAmpl01) )
            val = 0.;
          if ( rms01 > RMSThreshold_ )
            val = 0.;
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea01_[ism-1] ) {
              if ( mean01 > 0. ) {
                mea01_[ism-1]->setBinContent( ic, mean01 );
                mea01_[ism-1]->setBinError( ic, rms01 );
              } else {
                mea01_[ism-1]->setEntries( 1.+mea01_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update02 ) {

          float val;

          val = 1.;
          if ( fabs(mean02 - meanAmpl02) > fabs(percentVariation_ * meanAmpl02) )
            val = 0.;
          if ( rms02 > RMSThreshold_ )
            val = 0.;
          if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea02_[ism-1] ) {
              if ( mean02 > 0. ) {
                mea02_[ism-1]->setBinContent( ic, mean02 );
                mea02_[ism-1]->setBinError( ic, rms02 );
              } else {
                mea02_[ism-1]->setEntries( 1.+mea02_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update03 ) {

          float val;

          val = 1.;
          if ( fabs(mean03 - meanAmpl03) > fabs(percentVariation_ * meanAmpl03) )
            val = 0.;
          if ( rms03 > RMSThreshold_ )
            val = 0.;
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea03_[ism-1] ) {
              if ( mean03 > 0. ) {
                mea03_[ism-1]->setBinContent( ic, mean03 );
                mea03_[ism-1]->setBinError( ic, rms03 );
              } else {
                mea03_[ism-1]->setEntries( 1.+mea03_[ism-1]->getEntries() );
              }
            }
          }

        }

        // masking

        if ( mask1.size() != 0 ) {
          map<EcalLogicID, RunCrystalErrorsDat>::const_iterator m;
          for (m = mask1.begin(); m != mask1.end(); m++) {

            int jx = ix + Numbers::ix0EE(ism);
            int jy = iy + Numbers::iy0EE(ism);

            if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

            if ( ! Numbers::validEE(ism, jx, jy) ) continue;

            int ic = Numbers::indexEE(ism, ix, iy);

            if ( ic == -1 ) continue;

            EcalLogicID ecid = m->first;

            if ( ecid.getID1() == Numbers::iSM(ism, EcalEndcap) && ecid.getID2() == ic ) {
              if ( (m->second).getErrorBits() & bits01 ) {
                if ( meg01_[ism-1] ) {
                  float val = int(meg01_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg01_[ism-1]->setBinContent( ix, iy, val+3 );
                }
              }
              if ( (m->second).getErrorBits() & bits02 ) {
                if ( meg02_[ism-1] ) {
                  float val = int(meg02_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg02_[ism-1]->setBinContent( ix, iy, val+3 );
                }
              }
              if ( (m->second).getErrorBits() & bits03 ) {
                if ( meg03_[ism-1] ) {
                  float val = int(meg03_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg03_[ism-1]->setBinContent( ix, iy, val+3 );
                }
              }
            }

          }
        }

      }
    }

    for ( int i = 1; i <= 10; i++ ) {

      if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg05_[ism-1] ) meg05_[ism-1]->setBinContent( i, 1, 2. );

      bool update01;
      bool update02;
      bool update03;
      bool update04;

      float num01, num02, num03, num04;
      float mean01, mean02, mean03, mean04;
      float rms01, rms02, rms03, rms04;

      update01 = UtilsClient::getBinStats(i01_[ism-1], i, 0, num01, mean01, rms01);
      update02 = UtilsClient::getBinStats(i02_[ism-1], i, 0, num02, mean02, rms02);
      update03 = UtilsClient::getBinStats(i03_[ism-1], i, 0, num03, mean03, rms03);
      update04 = UtilsClient::getBinStats(i04_[ism-1], i, 0, num04, mean04, rms04);

      if ( mer04_[ism-1] ) mer04_[ism-1]->Fill(rms03);
      if ( mer05_[ism-1] ) mer05_[ism-1]->Fill(rms04);

      if ( update01 && update03 ) {

        float val;

        val = 1.;
        if ( mean01 < amplitudeThresholdPnG01_ )
          val = 0.;
        if ( mean03 <  pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0] ||
          pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0] < mean03)
          val = 0.;
        if ( rms03 > pedPnRMSThreshold_[0] )
          val = 0.;
        if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent(i, 1, val);

      }

      if ( update02 && update04 ) {

        float val;

        val = 1.;
        if ( mean02 < amplitudeThresholdPnG16_ )
          val = 0.;
        if ( mean04 <  pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1] ||
          pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1] < mean04)
          val = 0.;
        if ( rms04 > pedPnRMSThreshold_[1] )
          val = 0.;
        if ( meg05_[ism-1] ) meg05_[ism-1]->setBinContent(i, 1, val);

      }

      // masking

      if ( mask2.size() != 0 ) {
        map<EcalLogicID, RunPNErrorsDat>::const_iterator m;
        for (m = mask2.begin(); m != mask2.end(); m++) {

          EcalLogicID ecid = m->first;

          if ( ecid.getID1() == Numbers::iSM(ism, EcalEndcap) && ecid.getID2() == i-1 ) {
            if ( (m->second).getErrorBits() & (bits01|bits04) ) {
              if ( meg04_[ism-1] ) {
                float val = int(meg04_[ism-1]->getBinContent(i, 1)) % 3;
                meg04_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits03|bits06) ) {
              if ( meg05_[ism-1] ) {
                float val = int(meg05_[ism-1]->getBinContent(i, 1)) % 3;
                meg05_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
          }

        }
      }

    }

    for ( int i = 1; i <= 10; i++ ) {

      if ( hs01_[ism-1] ) {
        mes01_[ism-1]->setBinContent( i, hs01_[ism-1]->GetBinContent(1, i) );
        mes01_[ism-1]->setBinError( i, hs01_[ism-1]->GetBinError(1, i) );
      }

      if ( hs02_[ism-1] ) {
        mes02_[ism-1]->setBinContent( i, hs02_[ism-1]->GetBinContent(1, i) );
        mes02_[ism-1]->setBinError( i, hs02_[ism-1]->GetBinError(1, i) );
      }

      if ( hs03_[ism-1] ) {
        mes03_[ism-1]->setBinContent( i, hs03_[ism-1]->GetBinContent(1, i) );
        mes03_[ism-1]->setBinError( i, hs03_[ism-1]->GetBinError(1, i) );
      }

    }

  }

}

void EETestPulseClient::htmlOutput(int run, string htmlDir, string htmlName){

  cout << "Preparing EETestPulseClient html output ..." << endl;

  ofstream htmlFile;

  htmlFile.open((htmlDir + htmlName).c_str());

  // html page header
  htmlFile << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">  " << endl;
  htmlFile << "<html>  " << endl;
  htmlFile << "<head>  " << endl;
  htmlFile << "  <meta content=\"text/html; charset=ISO-8859-1\"  " << endl;
  htmlFile << " http-equiv=\"content-type\">  " << endl;
  htmlFile << "  <title>Monitor:TestPulseTask output</title> " << endl;
  htmlFile << "</head>  " << endl;
  htmlFile << "<style type=\"text/css\"> td { font-weight: bold } </style>" << endl;
  htmlFile << "<body>  " << endl;
  //htmlFile << "<br>  " << endl;
  htmlFile << "<a name=""top""></a>" << endl;
  htmlFile << "<h2>Run:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" << endl;
  htmlFile << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">" << run << "</span></h2>" << endl;
  htmlFile << "<h2>Monitoring task:&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">TEST PULSE</span></h2> " << endl;
  htmlFile << "<hr>" << endl;
  htmlFile << "<table border=1><tr><td bgcolor=red>channel has problems in this task</td>" << endl;
  htmlFile << "<td bgcolor=lime>channel has NO problems</td>" << endl;
  htmlFile << "<td bgcolor=yellow>channel is missing</td></table>" << endl;
  htmlFile << "<br>" << endl;
  htmlFile << "<table border=1>" << std::endl;
  for ( unsigned int i=0; i<superModules_.size(); i ++ ) {
    htmlFile << "<td bgcolor=white><a href=""#"
             << Numbers::sEE(superModules_[i]).c_str() << ">"
             << setfill( '0' ) << setw(2) << superModules_[i] << "</a></td>";
  }
  htmlFile << std::endl << "</table>" << std::endl;

  // Produce the plots to be shown as .png files from existing histograms

  const int csize = 250;

//  const double histMax = 1.e15;

  int pCol3[6] = { 301, 302, 303, 304, 305, 306 };

  TH2S labelGrid("labelGrid","label grid", 100, -2., 98., 100, -2., 98.);
  for ( short j=0; j<400; j++ ) {
    int x = 5*(1 + j%20);
    int y = 5*(1 + j/20);
    labelGrid.SetBinContent(x, y, Numbers::inTowersEE[j]);
  }
  labelGrid.SetMarkerSize(1);
  labelGrid.SetMinimum(0.1);

  TH2C dummy1( "dummy1", "dummy1 for sm mem", 10, 0, 10, 5, 0, 5 );
  for ( short i=0; i<2; i++ ) {
    int a = 2 + i*5;
    int b = 2;
    dummy1.Fill( a, b, i+1+68 );
  }
  dummy1.SetMarkerSize(2);
  dummy1.SetMinimum(0.1);

  string imgNameQual[3], imgNameAmp[3], imgNameShape[3], imgNameMEPnQual[2], imgNameMEPn[2], imgNameMEPnPed[2], imgNameMEPnPedRms[2], imgName, meName;

  TCanvas* cQual = new TCanvas("cQual", "Temp", 2*csize, 2*csize);
  TCanvas* cQualPN = new TCanvas("cQualPN", "Temp", 2*csize, csize);
  TCanvas* cAmp = new TCanvas("cAmp", "Temp", csize, csize);
  TCanvas* cShape = new TCanvas("cShape", "Temp", csize, csize);
  TCanvas* cPed = new TCanvas("cPed", "Temp", csize, csize);

  TH2F* obj2f;
  TH1F* obj1f;
  TProfile* objp;

  // Loop on barrel supermodules

  for ( unsigned int i=0; i<superModules_.size(); i ++ ) {

    int ism = superModules_[i];

    // Loop on gains

    for ( int iCanvas = 1 ; iCanvas <= 3 ; iCanvas++ ) {

      // Quality plots

      imgNameQual[iCanvas-1] = "";

      obj2f = 0;
      switch ( iCanvas ) {
        case 1:
          obj2f = UtilsClient::getHisto<TH2F*>( meg01_[ism-1] );
          break;
        case 2:
          obj2f = UtilsClient::getHisto<TH2F*>( meg02_[ism-1] );
          break;
        case 3:
          obj2f = UtilsClient::getHisto<TH2F*>( meg03_[ism-1] );
          break;
        default:
          break;
      }

      if ( obj2f ) {

        meName = obj2f->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1, "_");
          }
        }
        imgNameQual[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameQual[iCanvas-1];

        cQual->cd();
        gStyle->SetOptStat(" ");
        gStyle->SetPalette(6, pCol3);
        cQual->SetGridx();
        cQual->SetGridy();
        obj2f->GetXaxis()->SetLabelSize(0.02);
        obj2f->GetXaxis()->SetTitleSize(0.02);
        obj2f->GetYaxis()->SetLabelSize(0.02);
        obj2f->GetYaxis()->SetTitleSize(0.02);
        obj2f->SetMinimum(-0.00000001);
        obj2f->SetMaximum(6.0);
        obj2f->Draw("col");
        int x1 = labelGrid.GetXaxis()->FindBin(Numbers::ix0EE(ism)+0.);
        int x2 = labelGrid.GetXaxis()->FindBin(Numbers::ix0EE(ism)+50.);
        int y1 = labelGrid.GetYaxis()->FindBin(Numbers::iy0EE(ism)+0.);
        int y2 = labelGrid.GetYaxis()->FindBin(Numbers::iy0EE(ism)+50.);
        labelGrid.GetXaxis()->SetRange(x1, x2);
        labelGrid.GetYaxis()->SetRange(y1, y2);
        labelGrid.Draw("text,same");
        cQual->SetBit(TGraph::kClipFrame);
        TLine l;
        l.SetLineWidth(1);
        for ( int i=0; i<201; i=i+1){
          if ( (Numbers::ixSectorsEE[i]!=0 || Numbers::iySectorsEE[i]!=0) && (Numbers::ixSectorsEE[i+1]!=0 || Numbers::iySectorsEE[i+1]!=0) ) {
            l.DrawLine(Numbers::ixSectorsEE[i], Numbers::iySectorsEE[i], Numbers::ixSectorsEE[i+1], Numbers::iySectorsEE[i+1]);
          }
        }
        cQual->Update();
        cQual->SaveAs(imgName.c_str());

      }

      // Amplitude distributions

      imgNameAmp[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( mea01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( mea02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( mea03_[ism-1] );
          break;
        default:
          break;
      }

      if ( obj1f ) {

        meName = obj1f->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameAmp[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameAmp[iCanvas-1];

        cAmp->cd();
        gStyle->SetOptStat("euo");
        obj1f->SetStats(kTRUE);
//        if ( obj1f->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1f->SetMinimum(0.0);
        obj1f->Draw();
        cAmp->Update();
        gPad->SetLogy(0);
        cAmp->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Shape distributions

      imgNameShape[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( mes01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( mes02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( mes03_[ism-1] );  
          break;
        default: 
          break;
      }

      if ( obj1f ) {

        meName = obj1f->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1, "_");
          }
        }
        imgNameShape[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameShape[iCanvas-1];

        cShape->cd();
        gStyle->SetOptStat("euo");
        obj1f->SetStats(kTRUE);
//        if ( obj1f->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1f->Draw();
        cShape->Update();
        cShape->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

    }

    // Loop on gain

    for ( int iCanvas = 1 ; iCanvas <= 2 ; iCanvas++ ) {

      // Monitoring elements plots

      imgNameMEPnQual[iCanvas-1] = "";

      obj2f = 0;
      switch ( iCanvas ) {
      case 1:
        obj2f = UtilsClient::getHisto<TH2F*>( meg04_[ism-1] );
        break;
      case 2:
        obj2f = UtilsClient::getHisto<TH2F*>( meg05_[ism-1] );
        break;
      default:
        break;
      }

      if ( obj2f ) {

        meName = obj2f->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1, "_");
          }
        }
        imgNameMEPnQual[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnQual[iCanvas-1];

        cQualPN->cd();
        gStyle->SetOptStat(" ");
        gStyle->SetPalette(6, pCol3);
        obj2f->GetXaxis()->SetNdivisions(10);
        obj2f->GetYaxis()->SetNdivisions(5);
        cQualPN->SetGridx();
        cQualPN->SetGridy(0);
        obj2f->SetMinimum(-0.00000001);
        obj2f->SetMaximum(6.0);
        obj2f->Draw("col");
        dummy1.Draw("text,same");
        cQualPN->Update();
        cQualPN->SaveAs(imgName.c_str());

      }

      imgNameMEPn[iCanvas-1] = "";

      objp = 0;
      switch ( iCanvas ) {
        case 1:
          objp = i01_[ism-1];
          break;
        case 2:
          objp = i02_[ism-1];
          break;
        default:
          break;
      }

      if ( objp ) {

        meName = objp->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPn[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPn[iCanvas-1];

        cAmp->cd();
        gStyle->SetOptStat("euo");
        objp->SetStats(kTRUE);
//        if ( objp->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        objp->SetMinimum(0.0);
        objp->Draw();
        cAmp->Update();
        cAmp->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Monitoring elements plots

      imgNameMEPnPed[iCanvas-1] = "";

      objp = 0;
      switch ( iCanvas ) {
        case 1:
          objp = i03_[ism-1];
          break;
        case 2:
          objp = i04_[ism-1];
          break;
        default:
          break;
      }

      if ( objp ) {

        meName = objp->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnPed[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnPed[iCanvas-1];

        cPed->cd();
        gStyle->SetOptStat("euo");
        objp->SetStats(kTRUE);
//        if ( objp->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        objp->SetMinimum(0.0);
        objp->Draw();
        cPed->Update();
        cPed->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      imgNameMEPnPedRms[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
      case 1:
        if ( mer04_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mer04_[ism-1]);
        break;
      case 2:
        if ( mer05_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mer05_[ism-1]);
        break;
      default:
        break;
      }

      if ( obj1f ) {

        meName = obj1f->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnPedRms[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnPedRms[iCanvas-1];

        cPed->cd();
        gStyle->SetOptStat("euo");
        obj1f->SetStats(kTRUE);
        //        if ( obj1f->GetMaximum(histMax) > 0. ) {
        //          gPad->SetLogy(1);
        //        } else {
        //          gPad->SetLogy(0);
        //        }
        obj1f->SetMinimum(0.0);
        obj1f->Draw();
        cPed->Update();
        cPed->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

    }

    if( i>0 ) htmlFile << "<a href=""#top"">Top</a>" << std::endl;
    htmlFile << "<hr>" << std::endl;
    htmlFile << "<h3><a name="""
             << Numbers::sEE(ism).c_str() << """></a><strong>"
             << Numbers::sEE(ism).c_str() << "</strong></h3>" << endl;
    htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
    htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 3 ; iCanvas++ ) {

      if ( imgNameQual[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameQual[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;
    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 3 ; iCanvas++ ) {

      if ( imgNameAmp[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameAmp[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameShape[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameShape[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr align=\"center\"><td colspan=\"2\">Gain 1</td><td colspan=\"2\">Gain 6</td><td colspan=\"2\">Gain 12</td></tr>" << endl;
    htmlFile << "</table>" << endl;
    htmlFile << "<br>" << endl;

    htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
    htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 2 ; iCanvas++ ) {

      if ( imgNameMEPnQual[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPnQual[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;
    htmlFile << "</table>" << endl;
    htmlFile << "<br>" << endl;

    htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
    htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 2 ; iCanvas++ ) {

      if ( imgNameMEPnPed[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPnPed[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPnPedRms[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPnPedRms[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPn[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPn[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr align=\"center\">  <td colspan=\"2\"> </td>  <td colspan=\"2\">Gain 1</td>  <td colspan=\"2\"> </td> <td colspan=\"2\"> </td> <td colspan=\"2\">Gain 16</td></tr>" << endl;
    htmlFile << "</table>" << endl;
    htmlFile << "<br>" << endl;

  }

  delete cQual;
  delete cQualPN;
  delete cAmp;
  delete cShape;
  delete cPed;

  // html page footer
  htmlFile << "</body> " << endl;
  htmlFile << "</html> " << endl;

  htmlFile.close();

}

