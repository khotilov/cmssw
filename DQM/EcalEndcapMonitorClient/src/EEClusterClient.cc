/*
 * \file EEClusterClient.cc
 *
 * $Date: 2007/05/12 09:48:34 $
 * $Revision: 1.7 $
 * \author G. Della Ricca
 * \author F. Cossutti
 * \author E. Di Marco
 *
*/

#include <memory>
#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TGaxis.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DaqMonitorBEInterface.h"
#include "DQMServices/UI/interface/MonitorUIRoot.h"
#include "DQMServices/Core/interface/QTestStatus.h"
#include "DQMServices/QualityTests/interface/QCriterionRoot.h"

#include "OnlineDB/EcalCondDB/interface/RunTag.h"
#include "OnlineDB/EcalCondDB/interface/RunIOV.h"
#include "OnlineDB/EcalCondDB/interface/MonPedestalsOnlineDat.h"

#include <DQM/EcalEndcapMonitorClient/interface/EEClusterClient.h>
#include <DQM/EcalCommon/interface/UtilsClient.h>

using namespace cms;
using namespace edm;
using namespace std;

EEClusterClient::EEClusterClient(const ParameterSet& ps){

  // collateSources switch
  collateSources_ = ps.getUntrackedParameter<bool>("collateSources", false);

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  // enableQT switch
  enableQT_ = ps.getUntrackedParameter<bool>("enableQT", true);

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  // MonitorDaemon switch
  enableMonitorDaemon_ = ps.getUntrackedParameter<bool>("enableMonitorDaemon", true);

  // prefix to ME paths
  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  // vector of selected Super Modules (Defaults to all 18).
  superModules_.reserve(18);
  for ( unsigned int i = 1; i < 19; i++ ) superModules_.push_back(i);
  superModules_ = ps.getUntrackedParameter<vector<int> >("superModules", superModules_);

  h01_[0] = 0;
  h01_[1] = 0;
  h01_[2] = 0;

  h02_[0] = 0;
  h02_[1] = 0;

  h03_ = 0;
  h04_ = 0;

  i01_[0] = 0;
  i01_[1] = 0;
  i01_[2] = 0;

  i02_[0] = 0;
  i02_[1] = 0;

  i03_ = 0;
  i04_ = 0;

  s01_[0] = 0;
  s01_[1] = 0;

}

EEClusterClient::~EEClusterClient(){

}

void EEClusterClient::beginJob(MonitorUserInterface* mui){

  mui_ = mui;

  if ( verbose_ ) cout << "EEClusterClient: beginJob" << endl;

  ievt_ = 0;
  jevt_ = 0;

  if ( enableQT_ ) {


  }

}

void EEClusterClient::beginRun(void){

  if ( verbose_ ) cout << "EEClusterClient: beginRun" << endl;

  jevt_ = 0;

  this->setup();

  this->subscribe();

}

void EEClusterClient::endJob(void) {

  if ( verbose_ ) cout << "EEClusterClient: endJob, ievt = " << ievt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EEClusterClient::endRun(void) {

  if ( verbose_ ) cout << "EEClusterClient: endRun, jevt = " << jevt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EEClusterClient::setup(void) {

  mui_->setCurrentFolder( "EcalEndcap/EEClusterClient" );

}

void EEClusterClient::cleanup(void) {

  if ( cloneME_ ) {
    if ( h01_[0] ) delete h01_[0];
    if ( h01_[1] ) delete h01_[1];
    if ( h01_[2] ) delete h01_[2];

    if ( h02_[0] ) delete h02_[0];
    if ( h02_[1] ) delete h02_[1];

    if ( h03_ ) delete h03_;
    if ( h04_ ) delete h04_;

    if ( i01_[0] ) delete i01_[0];
    if ( i01_[1] ) delete i01_[1];
    if ( i01_[2] ) delete i01_[2];

    if ( i02_[0] ) delete i02_[0];
    if ( i02_[1] ) delete i02_[1];

    if ( i03_ ) delete i03_;
    if ( i04_ ) delete i04_;

    if ( s01_[0] ) delete s01_[0];
    if ( s01_[1] ) delete s01_[1];

  }

  h01_[0] = 0;
  h01_[1] = 0;
  h01_[2] = 0;

  h02_[0] = 0;
  h02_[1] = 0;

  h03_ = 0;
  h04_ = 0;

  i01_[0] = 0;
  i01_[1] = 0;
  i01_[2] = 0;

  i02_[0] = 0;
  i02_[1] = 0;

  i03_ = 0;
  i04_ = 0;

  s01_[0] = 0;
  s01_[1] = 0;

}

bool EEClusterClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov, int ism) {

  bool status = true;

  return status;

}

void EEClusterClient::subscribe(void){

  if ( verbose_ ) cout << "EEClusterClient: subscribe" << endl;

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster crystals");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster ET map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster size map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster ET map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size map");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT hybrid S1toE");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
  mui_->subscribe(histo);

  if ( collateSources_ ) {

    if ( verbose_ ) cout << "EEClusterClient: collate" << endl;

    sprintf(histo, "EECLT island basic cluster energy");
    me_h01_[0] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy");
    mui_->add(me_h01_[0], histo);

    sprintf(histo, "EECLT island basic cluster number");
    me_h01_[1] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number");
    mui_->add(me_h01_[1], histo);

    sprintf(histo, "EECLT island basic cluster crystals");
    me_h01_[2] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster crystals");
    mui_->add(me_h01_[2], histo);

    sprintf(histo, "EECLT island basic cluster energy map");
    me_h02_[0] = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy map");
    mui_->add(me_h02_[0], histo);

    sprintf(histo, "EECLT island basic cluster ET map");
    me_h02_[1] = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster ET map");
    mui_->add(me_h02_[1], histo);

    sprintf(histo, "EECLT island basic cluster number map");
    me_h03_ = mui_->collate2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number map");
    mui_->add(me_h03_, histo);

    sprintf(histo, "EECLT island basic cluster size map");
    me_h04_ = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster size map");
    mui_->add(me_h04_, histo);

    sprintf(histo, "EECLT island super cluster energy");
    me_i01_[0] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy");
    mui_->add(me_i01_[0], histo);

    sprintf(histo, "EECLT island super cluster number");
    me_i01_[1] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number");
    mui_->add(me_i01_[1], histo);

    sprintf(histo, "EECLT island super cluster size");
    me_i01_[2] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size");
    mui_->add(me_i01_[2], histo);

    sprintf(histo, "EECLT island super cluster energy map");
    me_i02_[0] = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy map");
    mui_->add(me_i02_[0], histo);

    sprintf(histo, "EECLT island super cluster ET map");
    me_i02_[1] = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster ET map");
    mui_->add(me_i02_[1], histo);

    sprintf(histo, "EECLT island super cluster number map");
    me_i03_ = mui_->collate2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number map");
    mui_->add(me_i03_, histo);

    sprintf(histo, "EECLT island super cluster size map");
    me_i04_ = mui_->collateProf2D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size map");
    mui_->add(me_i04_, histo);

    sprintf(histo, "EECLT hybrid S1toE");
    me_s01_[0] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT hybrid S1toE");
    mui_->add(me_s01_[0], histo);

    sprintf(histo, "EECLT dicluster invariant mass");
    me_s01_[1] = mui_->collate1D(histo, histo, "EcalEndcap/Sums/EEClusterTask");
    sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
    mui_->add(me_s01_[1], histo);
  }

}

void EEClusterClient::subscribeNew(void){

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster crystals");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster ET map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster size map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster ET map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size map");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT hybrid S1toE");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
  mui_->subscribeNew(histo);

}

void EEClusterClient::unsubscribe(void){

  if ( verbose_ ) cout << "EEClusterClient: unsubscribe" << endl;

  if ( collateSources_ ) {

    if ( verbose_ ) cout << "EEClusterClient: uncollate" << endl;

    if ( mui_ ) {

      mui_->removeCollate(me_h01_[0]);
      mui_->removeCollate(me_h01_[1]);
      mui_->removeCollate(me_h01_[2]);
      mui_->removeCollate(me_h02_[0]);
      mui_->removeCollate(me_h02_[1]);
      mui_->removeCollate(me_h03_);
      mui_->removeCollate(me_h04_);

      mui_->removeCollate(me_i01_[0]);
      mui_->removeCollate(me_i01_[1]);
      mui_->removeCollate(me_i01_[2]);
      mui_->removeCollate(me_i02_[0]);
      mui_->removeCollate(me_i02_[1]);
      mui_->removeCollate(me_i03_);
      mui_->removeCollate(me_i04_);

      mui_->removeCollate(me_s01_[0]);
      mui_->removeCollate(me_s01_[1]);

    }

  }

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster crystals");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster energy map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster ET map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster number map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island basic cluster size map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster energy map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster ET map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster number map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island super cluster size map");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT hybrid S1toE");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
  mui_->unsubscribe(histo);

}

void EEClusterClient::softReset(void){

}

void EEClusterClient::analyze(void){

  ievt_++;
  jevt_++;
  if ( ievt_ % 10 == 0 ) {
    if ( verbose_ ) cout << "EEClusterClient: ievt/jevt = " << ievt_ << "/" << jevt_ << endl;
  }

  Char_t histo[200];

  MonitorElement* me;

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster energy");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster energy").c_str());
  }
  me = mui_->get(histo);
  h01_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, h01_[0] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster number");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster number").c_str());
  }
  me = mui_->get(histo);
  h01_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, h01_[1] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster crystals");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster crystals").c_str());
  }
  me = mui_->get(histo);
  h01_[2] = UtilsClient::getHisto<TH1F*>( me, cloneME_, h01_[2] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster energy map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster energy map").c_str());
  }
  me = mui_->get(histo);
  h02_[0] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h02_[0] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster ET map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster ET map").c_str());
  }
  me = mui_->get(histo);
  h02_[1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h02_[1] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster number map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster number map").c_str());
  }
  me = mui_->get(histo);
  h03_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, h03_ );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island basic cluster size map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island basic cluster size map").c_str());
  }
  me = mui_->get(histo);
  h04_ = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h04_ );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster energy");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster energy").c_str());
  }
  me = mui_->get(histo);
  i01_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, i01_[0] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster number");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster number").c_str());
  }
  me = mui_->get(histo);
  i01_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, i01_[1] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster size");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster size").c_str());
  }
  me = mui_->get(histo);
  i01_[2] = UtilsClient::getHisto<TH1F*>( me, cloneME_, i01_[2] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster energy map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster energy map").c_str());
  }
  me = mui_->get(histo);
  i02_[0] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i02_[0] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster ET map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster ET map").c_str());
  }
  me = mui_->get(histo);
  i02_[1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i02_[1] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster number map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster number map").c_str());
  }
  me = mui_->get(histo);
  i03_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, i03_ );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT island super cluster size map");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island super cluster size map").c_str());
  }
  me = mui_->get(histo);
  i04_ = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i04_ );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT hybrid S1toE");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT hybrid S1toE").c_str());
  }
  me = mui_->get(histo);
  s01_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, s01_[0] );

  if ( collateSources_ ) {
    sprintf(histo, "EcalEndcap/Sums/EEClusterTask/EECLT dicluster invariant mass");
  } else {
    sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT dicluster invariant mass").c_str());
  }
  me = mui_->get(histo);
  s01_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, s01_[1] );

}

void EEClusterClient::htmlOutput(int run, string htmlDir, string htmlName){

  cout << "Preparing EEClusterClient html output ..." << endl;

  ofstream htmlFile;

  htmlFile.open((htmlDir + htmlName).c_str());

  // html page header
  htmlFile << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">  " << endl;
  htmlFile << "<html>  " << endl;
  htmlFile << "<head>  " << endl;
  htmlFile << "  <meta content=\"text/html; charset=ISO-8859-1\"  " << endl;
  htmlFile << " http-equiv=\"content-type\">  " << endl;
  htmlFile << "  <title>Monitor:ClusterTask output</title> " << endl;
  htmlFile << "</head>  " << endl;
  htmlFile << "<style type=\"text/css\"> td { font-weight: bold } </style>" << endl;
  htmlFile << "<body>  " << endl;
  htmlFile << "<br>  " << endl;
  htmlFile << "<h2>Run:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" << endl;
  htmlFile << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">" << run << "</span></h2>" << endl;
  htmlFile << "<h2>Monitoring task:&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">CLUSTER</span></h2> " << endl;
  htmlFile << "<hr>" << endl;
  //  htmlFile << "<table border=1><tr><td bgcolor=red>channel has problems in this task</td>" << endl;
  //  htmlFile << "<td bgcolor=lime>channel has NO problems</td>" << endl;
  //  htmlFile << "<td bgcolor=yellow>channel is missing</td></table>" << endl;
  //  htmlFile << "<hr>" << endl;

  htmlFile <<  "<a href=\"#bc_plots\"> Basic Clusters plots </a>" << endl;
  htmlFile << "<p>" << endl;
  htmlFile <<  "<a href=\"#sc_plots\"> Super Clusters plots </a>" << endl;
  htmlFile << "<p>" << endl;
  htmlFile <<  "<a href=\"#hl_plots\"> Higher Level Quantities plots </a>" << endl;
  htmlFile << "<p>" << endl;

  htmlFile << "<hr>" << endl;
  htmlFile << "<p>" << endl;

  // Produce the plots to be shown as .png files from existing histograms

  const int csize1D = 250;
  const int csize2D = 300;

  //  const double histMax = 1.e15;

  int pCol4[10];
  for ( int i = 0; i < 10; i++ ) pCol4[i] = 401+i;

  // dummy histogram labelling the SM's
  TH2C labelGrid("labelGrid","label grid for SM", 9, -M_PI*(9+0.5)/9, M_PI*(9-0.5)/9, 2, -1.479, 1.479);
  for ( short sm=0; sm<18; sm++ ) {
    int x = 1 + sm%9;
    int y = 2 - sm/9;
    int z = x + 4;
    if ( z > 9 ) z = z - 9;
    if ( y == 1 ) {
      labelGrid.SetBinContent(x, y, -z);
    } else {
      labelGrid.SetBinContent(x, y, +z);
    }
  }
  labelGrid.SetMarkerSize(2);
  labelGrid.SetMinimum(-9.01);

  TGaxis axis(-M_PI*(9+0.5)/9, -1.479, M_PI*(9-0.5)/9, -1.479, -M_PI*(9+0.5)/9, M_PI*(9-0.5)/9, 40306, "N");

  string imgNameB[3], imgNameBMap[4], imgNameS[3], imgNameSMap[4];
  string imgNameBXproj[4], imgNameBYproj[4], imgNameSXproj[4], imgNameSYproj[4];
  string imgNameHL[2], imgName, meName;

  TCanvas* cEne = new TCanvas("cEne", "Temp", csize1D, csize1D);
  TCanvas* cMap = new TCanvas("cMap", "Temp", int(360./170.*csize2D), csize2D);

  TH1F* obj1f = 0;
  TProfile2D* objp;
  TH2F* obj2f = 0;
  TH1D* obj1dX;
  TH1D* obj1dY;

  // ==========================================================================
  // basic clusters
  // ==========================================================================

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameB[iCanvas-1] = "";

    obj1f = h01_[iCanvas-1];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }
      imgNameB[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameB[iCanvas-1];

      cEne->cd();
      gStyle->SetOptStat("euomr");
      obj1f->SetStats(kTRUE);
      obj1f->Draw();
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  ///*****************************************************************************///
  htmlFile << "<br>" << endl;
  htmlFile <<  "<a name=\"bc_plots\"> <B> Basic Clusters plots </B> </a> " << endl;
  htmlFile << "</br>" << endl;
  ///*****************************************************************************///

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    if ( imgNameB[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameB[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameBMap[iCanvas-1] = "";

    objp = (iCanvas!=3) ? h02_[iCanvas-1] : h04_;

    if ( objp ) {

      meName = objp->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1 ,"_" );
        }
      }
      imgNameBMap[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameBMap[iCanvas-1];

      cMap->cd();
      gStyle->SetOptStat(" ");
      gStyle->SetPalette(10, pCol4);
      objp->GetXaxis()->SetNdivisions( 40109, kFALSE);
      objp->GetYaxis()->SetNdivisions(170102, kFALSE);
      cMap->SetGridx();
      cMap->SetGridy();
      objp->Draw("colz");
      labelGrid.Draw("text,same");
      axis.Draw();
      cMap->Update();
      objp->GetXaxis()->SetLabelColor(0);
      cMap->SaveAs(imgName.c_str());
      objp->GetXaxis()->SetLabelColor(1);

      char projXName[100];
      char projYName[100];
      sprintf(projXName,"%s_px",meName.c_str());
      imgNameBXproj[iCanvas-1] = string(projXName) + ".png";
      sprintf(projYName,"%s_py",meName.c_str());
      imgNameBYproj[iCanvas-1] = string(projYName) + ".png";

      obj1dX = objp->ProjectionX(projXName,1,objp->GetNbinsY(),"e");
      obj1dY = objp->ProjectionY(projYName,1,objp->GetNbinsX(),"e");

      cEne->cd();
      gStyle->SetOptStat("emr");
      obj1dX->GetXaxis()->SetNdivisions(40306, kFALSE);
      obj1dY->GetXaxis()->SetNdivisions(6, kFALSE);

      imgName = htmlDir + imgNameBXproj[iCanvas-1];
      obj1dX->SetStats(kTRUE);
      obj1dX->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());

      imgName = htmlDir + imgNameBYproj[iCanvas-1];
      obj1dY->SetStats(kTRUE);
      obj1dY->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  imgNameBMap[3] = "";

  obj2f = h03_;

  if ( obj2f ) {

    meName = obj2f->GetName();

    for ( unsigned int i = 0; i < meName.size(); i++ ) {
      if ( meName.substr(i, 1) == " " )  {
        meName.replace(i, 1 ,"_" );
      }
    }
    imgNameBMap[3] = meName + ".png";
    imgName = htmlDir + imgNameBMap[3];

    cMap->cd();
    gStyle->SetOptStat(" ");
    gStyle->SetPalette(10, pCol4);
    obj2f->GetXaxis()->SetNdivisions( 40109, kFALSE);
    obj2f->GetYaxis()->SetNdivisions(170102, kFALSE);
    cMap->SetGridx();
    cMap->SetGridy();
    obj2f->Draw("colz");
    labelGrid.Draw("text,same");
    axis.Draw();
    cMap->Update();
    obj2f->GetXaxis()->SetLabelColor(0);
    cMap->SaveAs(imgName.c_str());
    obj2f->GetXaxis()->SetLabelColor(1);

    char projXName[100];
    char projYName[100];
    sprintf(projXName,"%s_px",meName.c_str());
    imgNameBXproj[3] = string(projXName) + ".png";
    sprintf(projYName,"%s_py",meName.c_str());
    imgNameBYproj[3] = string(projYName) + ".png";

    obj1dX = obj2f->ProjectionX(projXName,1,obj2f->GetNbinsY(),"e");
    obj1dY = obj2f->ProjectionY(projYName,1,obj2f->GetNbinsX(),"e");

    cEne->cd();
    gStyle->SetOptStat("emr");
    obj1dX->GetXaxis()->SetNdivisions(40306, kFALSE);
    obj1dY->GetXaxis()->SetNdivisions(6, kFALSE);

    imgName = htmlDir + imgNameBXproj[3];
    obj1dX->SetStats(kTRUE);
    obj1dX->Draw("pe");
    cEne->Update();
    cEne->SaveAs(imgName.c_str());

    imgName = htmlDir + imgNameBYproj[3];
    obj1dY->SetStats(kTRUE);
    obj1dY->Draw("pe");
    cEne->Update();
    cEne->SaveAs(imgName.c_str());

  }

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 2; iCanvas++ ) {

    if ( imgNameBMap[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameBMap[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 3; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameBMap[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameBMap[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // projections X...
  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameBXproj[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameBXproj[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // projections Y...
  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameBYproj[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameBYproj[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // ====================================================================
  // super clusters
  // ====================================================================

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameS[iCanvas-1] = "";

    obj1f = i01_[iCanvas-1];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }
      imgNameS[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameS[iCanvas-1];

      cEne->cd();
      gStyle->SetOptStat("euomr");
      obj1f->SetStats(kTRUE);
      obj1f->Draw();
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  ///*****************************************************************************///
  htmlFile << "<br>" << endl;
  htmlFile <<  "<a name=\"sc_plots\"> <B> Super Clusters plots </B> </a> " << endl;
  htmlFile << "</br>" << endl;
  ///*****************************************************************************///

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    if ( imgNameS[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameS[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameSMap[iCanvas-1] = "";

    objp = (iCanvas!=3) ? i02_[iCanvas-1] : i04_;

    if ( objp ) {

      meName = objp->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1 ,"_" );
        }
      }
      imgNameSMap[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameSMap[iCanvas-1];

      cMap->cd();
      gStyle->SetOptStat(" ");
      gStyle->SetPalette(10, pCol4);
      objp->GetXaxis()->SetNdivisions( 40109, kFALSE);
      objp->GetYaxis()->SetNdivisions(170102, kFALSE);
      cMap->SetGridx();
      cMap->SetGridy();
      objp->Draw("colz");
      labelGrid.Draw("text,same");
      axis.Draw();
      cMap->Update();
      objp->GetXaxis()->SetLabelColor(0);
      cMap->SaveAs(imgName.c_str());
      objp->GetXaxis()->SetLabelColor(1);

      char projXName[100];
      char projYName[100];
      sprintf(projXName,"%s_px",meName.c_str());
      imgNameSXproj[iCanvas-1] = string(projXName) + ".png";
      sprintf(projYName,"%s_py",meName.c_str());
      imgNameSYproj[iCanvas-1] = string(projYName) + ".png";

      obj1dX = objp->ProjectionX(projXName,1,objp->GetNbinsY(),"e");
      obj1dY = objp->ProjectionY(projYName,1,objp->GetNbinsX(),"e");

      cEne->cd();
      gStyle->SetOptStat("emr");
      obj1dX->GetXaxis()->SetNdivisions(40306, kFALSE);
      obj1dY->GetXaxis()->SetNdivisions(6, kFALSE);

      imgName = htmlDir + imgNameSXproj[iCanvas-1];
      obj1dX->SetStats(kTRUE);
      obj1dX->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());

      imgName = htmlDir + imgNameSYproj[iCanvas-1];
      obj1dY->SetStats(kTRUE);
      obj1dY->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  imgNameSMap[3] = "";

  obj2f = i03_;

  if ( obj2f ) {

    meName = obj2f->GetName();

    for ( unsigned int i = 0; i < meName.size(); i++ ) {
      if ( meName.substr(i, 1) == " " )  {
        meName.replace(i, 1 ,"_" );
      }
    }
    imgNameSMap[3] = meName + ".png";
    imgName = htmlDir + imgNameSMap[3];

    cMap->cd();
    gStyle->SetOptStat(" ");
    gStyle->SetPalette(10, pCol4);
    obj2f->GetXaxis()->SetNdivisions( 40109, kFALSE);
    obj2f->GetYaxis()->SetNdivisions(170102, kFALSE);
    cMap->SetGridx();
    cMap->SetGridy();
    obj2f->Draw("colz");
    labelGrid.Draw("text,same");
    axis.Draw();
    cMap->Update();
    obj2f->GetXaxis()->SetLabelColor(0);
    cMap->SaveAs(imgName.c_str());
    obj2f->GetXaxis()->SetLabelColor(1);

    char projXName[100];
    char projYName[100];
    sprintf(projXName,"%s_px",meName.c_str());
    imgNameSXproj[3] = string(projXName) + ".png";
    sprintf(projYName,"%s_py",meName.c_str());
    imgNameSYproj[3] = string(projYName) + ".png";

    obj1dX = obj2f->ProjectionX(projXName,1,obj2f->GetNbinsY(),"e");
    obj1dY = obj2f->ProjectionY(projYName,1,obj2f->GetNbinsX(),"e");

    cEne->cd();
    gStyle->SetOptStat("emr");
    obj1dX->GetXaxis()->SetNdivisions(40306, kFALSE);
    obj1dY->GetXaxis()->SetNdivisions(6, kFALSE);

    imgName = htmlDir + imgNameSXproj[3];
    obj1dX->SetStats(kTRUE);
    obj1dX->Draw("pe");
    cEne->Update();
    cEne->SaveAs(imgName.c_str());

    imgName = htmlDir + imgNameSYproj[3];
    obj1dY->SetStats(kTRUE);
    obj1dY->Draw("pe");
    cEne->Update();
    cEne->SaveAs(imgName.c_str());

  }

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 2; iCanvas++ ) {

    if ( imgNameSMap[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameSMap[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 3; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameSMap[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameSMap[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // projections X...
  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameSXproj[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameSXproj[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // projections Y...
  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for ( int iCanvas = 1; iCanvas <= 4; iCanvas++ ) {

    if ( imgNameSYproj[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameSYproj[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // ===========================================================================
  // Higher Level variables
  // ===========================================================================

  for( int iCanvas = 1; iCanvas <= 2; iCanvas++ ) {

    imgNameHL[iCanvas-1] = "";

    obj1f = s01_[iCanvas-1];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }
      imgNameHL[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameHL[iCanvas-1];

      cEne->cd();
      gStyle->SetOptStat("euomr");
      obj1f->SetStats(kTRUE);
      obj1f->Draw();
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  ///*****************************************************************************///
  htmlFile << "<br>" << endl;
  htmlFile <<  "<a name=\"hl_plots\"> <B> Higher Level Quantities plots </B> </a> " << endl;
  htmlFile << "</br>" << endl;
  ///*****************************************************************************///

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  // for now only E1/Etot
  for ( int iCanvas = 1; iCanvas <= 2; iCanvas++ ) {

    if ( imgNameHL[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameHL[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  delete cEne;
  delete cMap;

  // html page footer
  htmlFile << "</body> " << endl;
  htmlFile << "</html> " << endl;

  htmlFile.close();

}

