/*
 * \file EEClusterClient.cc
 *
 * $Date: 2007/11/10 14:09:11 $
 * $Revision: 1.25 $
 * \author G. Della Ricca
 * \author E. Di Marco
 *
*/

#include <memory>
#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLine.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/UI/interface/MonitorUIRoot.h"

#include "OnlineDB/EcalCondDB/interface/RunTag.h"
#include "OnlineDB/EcalCondDB/interface/RunIOV.h"
#include "OnlineDB/EcalCondDB/interface/MonPedestalsOnlineDat.h"

#include "OnlineDB/EcalCondDB/interface/EcalCondDBInterface.h"

#include "CondTools/Ecal/interface/EcalErrorDictionary.h"

#include "DQM/EcalCommon/interface/EcalErrorMask.h"
#include "DQM/EcalCommon/interface/UtilsClient.h"
#include "DQM/EcalCommon/interface/LogicID.h"
#include "DQM/EcalCommon/interface/Numbers.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include <DQM/EcalEndcapMonitorClient/interface/EEClusterClient.h>

using namespace cms;
using namespace edm;
using namespace std;

EEClusterClient::EEClusterClient(const ParameterSet& ps){

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  // verbosity switch
  verbose_ = ps.getUntrackedParameter<bool>("verbose", false);

  // MonitorDaemon switch
  enableMonitorDaemon_ = ps.getUntrackedParameter<bool>("enableMonitorDaemon", true);

  // prefix to ME paths
  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  hBC1D_[0] = 0;
  hBC1D_[1] = 0;
  hBC1D_[2] = 0;

  for(int iEE=0;iEE<2;++iEE) {
    for(int i=0;i<3;++i) {
      hProfMap_[i][iEE] = 0;
      hProfMapProjR_[i][iEE] = 0;
      hProfMapProjPhi_[i][iEE] = 0;
    }
  }

  hOccMap_[0] = 0;
  hOccMapProjR_[0] = 0;
  hOccMapProjPhi_[0] = 0;

  hOccMap_[1] = 0;
  hOccMapProjR_[1] = 0;
  hOccMapProjPhi_[1] = 0;

  hSC1D_[0] = 0;
  hSC1D_[1] = 0;
  hSC1D_[2] = 0;

  s01_[0] = 0;
  s01_[1] = 0;
  s01_[2] = 0;

}

EEClusterClient::~EEClusterClient(){

}

void EEClusterClient::beginJob(MonitorUserInterface* mui){

  mui_ = mui;
  dbe_ = mui->getBEInterface();

  if ( verbose_ ) cout << "EEClusterClient: beginJob" << endl;

  ievt_ = 0;
  jevt_ = 0;

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

  dbe_->setCurrentFolder( "EcalEndcap/EEClusterClient" );

}

void EEClusterClient::cleanup(void) {

  if ( cloneME_ ) {
    if ( hBC1D_[0] ) delete hBC1D_[0];
    if ( hBC1D_[1] ) delete hBC1D_[1];
    if ( hBC1D_[2] ) delete hBC1D_[2];

    for(int iEE=0;iEE<2;++iEE) {
      for(int i=0;i<3;++i) {
        if(hProfMap_[i][iEE]) delete hProfMap_[i][iEE];
        if(hProfMapProjR_[i][iEE]) delete hProfMapProjR_[i][iEE];
        if(hProfMapProjPhi_[i][iEE]) delete hProfMapProjPhi_[i][iEE];
      }
    }

    if(hOccMap_[0]) delete hOccMap_[0];
    if(hOccMapProjR_[0]) delete hOccMapProjR_[0];
    if(hOccMapProjPhi_[0]) delete hOccMapProjPhi_[0];

    if(hOccMap_[1]) delete hOccMap_[1];
    if(hOccMapProjR_[1]) delete hOccMapProjR_[1];
    if(hOccMapProjPhi_[1]) delete hOccMapProjPhi_[1];

    if(hSC1D_[0]) delete hSC1D_[0];
    if(hSC1D_[1]) delete hSC1D_[1];
    if(hSC1D_[2]) delete hSC1D_[2];

    if(s01_[0]) delete s01_[0];
    if(s01_[1]) delete s01_[1];
    if(s01_[2]) delete s01_[2];

  }

  hBC1D_[0] = 0;
  hBC1D_[1] = 0;
  hBC1D_[2] = 0;

  for(int iEE=0;iEE<2;++iEE) {
    for(int i=0;i<3;++i) {
      hProfMap_[i][iEE] = 0;
      hProfMapProjR_[i][iEE] = 0;
      hProfMapProjPhi_[i][iEE] = 0;
    }
  }

  hOccMap_[0] = 0;
  hOccMapProjR_[0] = 0;
  hOccMapProjPhi_[0] = 0;

  hOccMap_[1] = 0;
  hOccMapProjR_[1] = 0;
  hOccMapProjPhi_[1] = 0;

  hSC1D_[0] = 0;
  hSC1D_[1] = 0;
  hSC1D_[2] = 0;

  s01_[0] = 0;
  s01_[1] = 0;
  s01_[2] = 0;

}

bool EEClusterClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov) {

  bool status = true;

  return status;

}

void EEClusterClient::subscribe(void){

  if ( verbose_ ) cout << "EEClusterClient: subscribe" << endl;

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE +");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE -");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC energy");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC size");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC number");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s1s9");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s9s25");
  mui_->subscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
  mui_->subscribe(histo);

}

void EEClusterClient::subscribeNew(void){

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE +");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE -");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC energy");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC size");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC number");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s1s9");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s9s25");
  mui_->subscribeNew(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT dicluster invariant mass");
  mui_->subscribeNew(histo);


}

void EEClusterClient::unsubscribe(void){

  if ( verbose_ ) cout << "EEClusterClient: unsubscribe" << endl;

  Char_t histo[200];

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE +");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy map EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number map EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET map EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size map EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection R EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection R EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC number projection phi EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection R EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection R EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT BC size projection phi EE -");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC energy");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC size");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT SC number");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s1s9");
  mui_->unsubscribe(histo);

  sprintf(histo, "*/EcalEndcap/EEClusterTask/EECLT island s9s25");
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

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy").c_str());
  me = dbe_->get(histo);
  hBC1D_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hBC1D_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size").c_str());
  me = dbe_->get(histo);
  hBC1D_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hBC1D_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number").c_str());
  me = dbe_->get(histo);
  hBC1D_[2] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hBC1D_[2] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy map EE +").c_str());
  me = dbe_->get(histo);
  hProfMap_[0][0] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[0][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number map EE +").c_str());
  me = dbe_->get(histo);
  hOccMap_[0] = UtilsClient::getHisto<TH2F*>( me, cloneME_, hOccMap_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET map EE +").c_str());
  me = dbe_->get(histo);
  hProfMap_[1][0] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[1][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size map EE +").c_str());
  me = dbe_->get(histo);
  hProfMap_[2][0] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[2][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy projection R EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[0][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[0][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[0][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[0][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number projection R EE +").c_str());
  me = dbe_->get(histo);
  hOccMapProjR_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hOccMapProjR_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number projection phi EE +").c_str());
  me = dbe_->get(histo);
  hOccMapProjPhi_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hOccMapProjPhi_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET projection R EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[1][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[1][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[1][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[1][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size projection R EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[2][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[2][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size projection phi EE +").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[2][0] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[2][0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy map EE -").c_str());
  me = dbe_->get(histo);
  hProfMap_[0][1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[0][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number map EE -").c_str());
  me = dbe_->get(histo);
  hOccMap_[1] = UtilsClient::getHisto<TH2F*>( me, cloneME_, hOccMap_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET map EE -").c_str());
  me = dbe_->get(histo);
  hProfMap_[1][1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[1][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size map EE -").c_str());
  me = dbe_->get(histo);
  hProfMap_[2][1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hProfMap_[2][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy projection R EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[0][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[0][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC energy projection phi EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[0][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[0][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number projection R EE -").c_str());
  me = dbe_->get(histo);
  hOccMapProjR_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hOccMapProjR_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC number projection phi EE -").c_str());
  me = dbe_->get(histo);
  hOccMapProjPhi_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hOccMapProjPhi_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET projection R EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[1][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[1][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC ET projection phi EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[1][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[1][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size projection R EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjR_[2][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjR_[2][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT BC size projection phi EE -").c_str());
  me = dbe_->get(histo);
  hProfMapProjPhi_[2][1] = UtilsClient::getHisto<TProfile*>( me, cloneME_, hProfMapProjPhi_[2][1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT SC energy").c_str());
  me = dbe_->get(histo);
  hSC1D_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hSC1D_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT SC size").c_str());
  me = dbe_->get(histo);
  hSC1D_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hSC1D_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT SC number").c_str());
  me = dbe_->get(histo);
  hSC1D_[2] = UtilsClient::getHisto<TH1F*>( me, cloneME_, hSC1D_[2] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island s1s9").c_str());
  me = dbe_->get(histo);
  s01_[0] = UtilsClient::getHisto<TH1F*>( me, cloneME_, s01_[0] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT island s9s25").c_str());
  me = dbe_->get(histo);
  s01_[1] = UtilsClient::getHisto<TH1F*>( me, cloneME_, s01_[1] );

  sprintf(histo, (prefixME_+"EcalEndcap/EEClusterTask/EECLT dicluster invariant mass").c_str());
  me = dbe_->get(histo);
  s01_[2] = UtilsClient::getHisto<TH1F*>( me, cloneME_, s01_[2] );

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
  const int csize2D = 500;

  const double histMax = 1.e15;

  int pCol4[10];
  for ( int i = 0; i < 10; i++ ) pCol4[i] = 401+i;

  TH2C labelGrid1("labelGrid1","label grid for EE -", 10, -150.0, 150.0, 10, -150.0, 150.0);

  for ( int i=1; i<=10; i++) {
    for ( int j=1; j<=10; j++) {
      labelGrid1.SetBinContent(i, j, -10);
    }
  }

  labelGrid1.SetBinContent(2, 5, -3);
  labelGrid1.SetBinContent(2, 7, -2);
  labelGrid1.SetBinContent(4, 9, -1);
  labelGrid1.SetBinContent(7, 9, -9);
  labelGrid1.SetBinContent(9, 7, -8);
  labelGrid1.SetBinContent(9, 5, -7);
  labelGrid1.SetBinContent(8, 3, -6);
  labelGrid1.SetBinContent(5, 2, -5);
  labelGrid1.SetBinContent(3, 3, -4);

  labelGrid1.SetMarkerSize(2);
  labelGrid1.SetMinimum(-9.01);
  labelGrid1.SetMaximum(-0.01);

  TH2C labelGrid2("labelGrid2","label grid for EE +", 10, -150.0, 150.0, 10, -150.0, 150.0);

  for ( int i=1; i<=10; i++) {
    for ( int j=1; j<=10; j++) {
      labelGrid2.SetBinContent(i, j, -10);
    }
  }

  labelGrid2.SetBinContent(2, 5, +7);
  labelGrid2.SetBinContent(2, 7, +8);
  labelGrid2.SetBinContent(4, 9, +9);
  labelGrid2.SetBinContent(7, 9, +1);
  labelGrid2.SetBinContent(9, 7, +2);
  labelGrid2.SetBinContent(9, 5, +3);
  labelGrid2.SetBinContent(8, 3, +4);
  labelGrid2.SetBinContent(6, 2, +5);
  labelGrid2.SetBinContent(3, 3, +6);

  labelGrid2.SetMarkerSize(2);
  labelGrid2.SetMinimum(+0.01);
  labelGrid2.SetMaximum(+9.01);

  TGaxis Xaxis(-150.0, -150.0,  150.0, -150.0, -150.0, 150.0, 50210, "N");
  TGaxis Yaxis(-150.0, -150.0, -150.0,  150.0, -150.0, 150.0, 50210, "N");
  Xaxis.SetLabelSize(0.02);
  Yaxis.SetLabelSize(0.02);

  string imgNameAll[3], imgNameEneMap[3][2], imgNameNumMap[2];
  string imgNameEneXproj[3][2], imgNameNumXproj[2], imgNameEneYproj[3][2], imgNameNumYproj[2];
  string imgNameHL[3], imgName, meName;

  TCanvas* cEne = new TCanvas("cEne", "Temp", csize1D, csize1D);
  TCanvas* cMap = new TCanvas("cMap", "Temp", csize2D, csize2D);

  TH1F* obj1f;
  TProfile2D* objp;
  TH2F* objf;
  TProfile* obj1pX;
  TProfile* obj1pY;
  TH1F* obj1dX;
  TH1F* obj1dY;

  gStyle->SetPaintTextFormat("+g");

  // ====> B A S I C     C L U S T E R S <===
  // ==========================================================================
  // all Ecal Endcap 1D plots
  // ==========================================================================

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameAll[iCanvas-1] = "";

    obj1f = hBC1D_[iCanvas-1];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }
      imgNameAll[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameAll[iCanvas-1];

      cEne->cd();
      gStyle->SetOptStat("euomr");
      obj1f->SetStats(kTRUE);
      if ( obj1f->GetMaximum(histMax) > 0. ) {
        gPad->SetLogy(1);
      } else {
        gPad->SetLogy(0);
      }
      obj1f->Draw();
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
      gPad->SetLogy(0);
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

    if ( imgNameAll[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameAll[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // Energy profiles
  for(int iVar=0; iVar<3; ++iVar) {
    for(int iEE=0; iEE<2; ++iEE) {

      imgNameEneMap[iVar][iEE] = "";

      objp = hProfMap_[iVar][iEE];

      if ( objp ) {

        meName = objp->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameEneMap[iVar][iEE] = meName + ".png";
        imgName = htmlDir + imgNameEneMap[iVar][iEE];

        cMap->cd();
        gStyle->SetOptStat(" ");
        gStyle->SetPalette(10, pCol4);
        objp->GetXaxis()->SetNdivisions(10, kFALSE);
        objp->GetXaxis()->SetLabelSize(0.02);
        objp->GetXaxis()->SetTitleSize(0.02);
        objp->GetYaxis()->SetNdivisions(10, kFALSE);
        objp->GetYaxis()->SetLabelSize(0.02);
        objp->GetYaxis()->SetTitleSize(0.02);
        objp->GetZaxis()->SetLabelSize(0.02);
        cMap->SetGridx();
        cMap->SetGridy();
        objp->Draw("colz");
        if ( iEE == 0 ) labelGrid1.Draw("text,same");
        if ( iEE == 1 ) labelGrid2.Draw("text,same");
        Xaxis.Draw();
        Yaxis.Draw();
        cMap->SetBit(TGraph::kClipFrame);
        TLine l;
        l.SetLineWidth(1);
        for ( int i=0; i<201; i=i+1){
          if ( (Numbers::ixSectorsEE[i]!=0 || Numbers::iySectorsEE[i]!=0) && (Numbers::ixSectorsEE[i+1]!=0 || Numbers::iySectorsEE[i+1]!=0) ) {
            l.DrawLine(3.0*(Numbers::ixSectorsEE[i]-50), 3.0*(Numbers::iySectorsEE[i]-50), 3.0*(Numbers::ixSectorsEE[i+1]-50), 3.0*(Numbers::iySectorsEE[i+1]-50));
          }
        }
        cMap->Update();
        objp->GetXaxis()->SetLabelColor(0);
        objp->GetYaxis()->SetLabelColor(0);
        cMap->SaveAs(imgName.c_str());
        objp->GetXaxis()->SetLabelColor(1);
        objp->GetYaxis()->SetLabelColor(1);
      }

      char projXName[100];
      char projYName[100];
      sprintf(projXName,"%s_px",meName.c_str());
      imgNameEneXproj[iVar][iEE] = string(projXName) + ".png";
      sprintf(projYName,"%s_py",meName.c_str());
      imgNameEneYproj[iVar][iEE] = string(projYName) + ".png";

      obj1pX = hProfMapProjR_[iVar][iEE];
      obj1pY = hProfMapProjPhi_[iVar][iEE];

      if(obj1pX && obj1pY) {
        cEne->cd();
        gStyle->SetOptStat("emr");
        obj1pX->GetXaxis()->SetNdivisions(50205, kFALSE);
        obj1pY->GetXaxis()->SetNdivisions(50206, kFALSE);

        imgName = htmlDir + imgNameEneXproj[iVar][iEE];
        obj1pX->SetStats(kTRUE);
        obj1pX->Draw("pe");
        cEne->Update();
        cEne->SaveAs(imgName.c_str());

        imgName = htmlDir + imgNameEneYproj[iVar][iEE];
        obj1pY->SetStats(kTRUE);
        obj1pY->Draw("pe");
        cEne->Update();
        cEne->SaveAs(imgName.c_str());
      }
    }
  }

  // Cluster occupancy profiles
  for ( int iEE=0; iEE<2; ++iEE ) {

    imgNameNumMap[iEE] = "";

    objf = hOccMap_[iEE];

    if ( objf ) {

      meName = objf->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1 ,"_" );
        }
      }
      imgNameNumMap[iEE] = meName + ".png";
      imgName = htmlDir + imgNameNumMap[iEE];

      cMap->cd();
      gStyle->SetOptStat(" ");
      gStyle->SetPalette(10, pCol4);
      objf->GetXaxis()->SetNdivisions(10, kFALSE);
      objf->GetXaxis()->SetLabelSize(0.02);
      objf->GetXaxis()->SetTitleSize(0.02);
      objf->GetYaxis()->SetNdivisions(10, kFALSE);
      objf->GetYaxis()->SetLabelSize(0.02);
      objf->GetYaxis()->SetTitleSize(0.02);
      objf->GetZaxis()->SetLabelSize(0.02);
      cMap->SetGridx();
      cMap->SetGridy();
      objf->Draw("colz");
      if ( iEE == 0 ) labelGrid1.Draw("text,same");
      if ( iEE == 1 ) labelGrid2.Draw("text,same");
      Xaxis.Draw();
      Yaxis.Draw();
      cMap->SetBit(TGraph::kClipFrame);
      TLine l;
      l.SetLineWidth(1);
      for ( int i=0; i<201; i=i+1){
        if ( (Numbers::ixSectorsEE[i]!=0 || Numbers::iySectorsEE[i]!=0) && (Numbers::ixSectorsEE[i+1]!=0 || Numbers::iySectorsEE[i+1]!=0) ) {
          l.DrawLine(3.0*(Numbers::ixSectorsEE[i]-50), 3.0*(Numbers::iySectorsEE[i]-50), 3.0*(Numbers::ixSectorsEE[i+1]-50), 3.0*(Numbers::iySectorsEE[i+1]-50));
        }
      }
      cMap->Update();
      objf->GetXaxis()->SetLabelColor(0);
      objf->GetYaxis()->SetLabelColor(0);
      cMap->SaveAs(imgName.c_str());
      objf->GetXaxis()->SetLabelColor(1);
      objf->GetYaxis()->SetLabelColor(1);
    }

    char projXName[100];
    char projYName[100];
    sprintf(projXName,"%s_px",meName.c_str());
    imgNameNumXproj[iEE] = string(projXName) + ".png";
    sprintf(projYName,"%s_py",meName.c_str());
    imgNameNumYproj[iEE] = string(projYName) + ".png";

    obj1dX = hOccMapProjR_[iEE];
    obj1dY = hOccMapProjPhi_[iEE];

    if(obj1dX && obj1dY) {
      cEne->cd();
      gStyle->SetOptStat("emr");
      obj1dX->GetXaxis()->SetNdivisions(50205, kFALSE);
      obj1dY->GetXaxis()->SetNdivisions(50206, kFALSE);

      imgName = htmlDir + imgNameNumXproj[iEE];
      obj1dX->SetStats(kTRUE);
      obj1dX->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());

      imgName = htmlDir + imgNameNumYproj[iEE];
      obj1dY->SetStats(kTRUE);
      obj1dY->Draw("pe");
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
    }
  }

  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for(int iVar=0; iVar<3; ++iVar) {
    for(int iEE=0; iEE<2; ++iEE) {
      if ( imgNameEneMap[iVar][iEE].size() != 0)
        htmlFile << "<td><img src=\"" << imgNameEneMap[iVar][iEE] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
    }
    htmlFile << "</tr>" << endl;
  }

  for ( int iEE = 0; iEE<2; ++iEE ) {
    if ( imgNameNumMap[iEE].size() != 0)
      htmlFile << "<td><img src=\"" << imgNameNumMap[iEE] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // projections...
  htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
  htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
  htmlFile << "<tr align=\"center\">" << endl;

  for(int iVar=0; iVar<3; ++iVar) {
    for(int iEE=0; iEE<2; ++iEE) {
      if ( imgNameEneXproj[iVar][iEE].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameEneXproj[iVar][iEE] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
      if ( imgNameEneYproj[iVar][iEE].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameEneYproj[iVar][iEE] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
    }
    htmlFile << "</tr>" << endl;
  }

  for(int iEE=0; iEE<2; ++iEE) {
    if ( imgNameNumXproj[iEE].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameNumXproj[iEE] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
    if ( imgNameNumYproj[iEE].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameNumYproj[iEE] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // ====>  S U P E R   C L U S T E R S   <====

  // ==========================================================================
  // all Ecal Endcap 1D plots
  // ==========================================================================

  for ( int iCanvas = 1; iCanvas <= 3; iCanvas++ ) {

    imgNameAll[iCanvas-1] = "";

    obj1f = hSC1D_[iCanvas-1];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }
      imgNameAll[iCanvas-1] = meName + ".png";
      imgName = htmlDir + imgNameAll[iCanvas-1];

      cEne->cd();
      gStyle->SetOptStat("euomr");
      obj1f->SetStats(kTRUE);
      if ( obj1f->GetMaximum(histMax) > 0. ) {
        gPad->SetLogy(1);
      } else {
        gPad->SetLogy(0);
      }
      obj1f->Draw();
      cEne->Update();
      cEne->SaveAs(imgName.c_str());
      gPad->SetLogy(0);
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

    if ( imgNameAll[iCanvas-1].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameAll[iCanvas-1] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  // ===========================================================================
  // Higher Level variables
  // ===========================================================================

  for(int iVar=0; iVar<3; ++iVar) {

    imgNameHL[iVar] = "";

    obj1f = s01_[iVar];

    if ( obj1f ) {

      meName = obj1f->GetName();

      for ( unsigned int i = 0; i < meName.size(); i++ ) {
        if ( meName.substr(i, 1) == " " )  {
          meName.replace(i, 1, "_");
        }
      }

      imgNameHL[iVar] = meName + ".png";
      imgName = htmlDir + imgNameHL[iVar];

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

  for(int iVar=0; iVar<3; ++iVar) {
    if ( imgNameHL[iVar].size() != 0 )
      htmlFile << "<td><img src=\"" << imgNameHL[iVar] << "\"></td>" << endl;
    else
      htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
  }

  htmlFile << "</tr>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "<br>" << endl;

  delete cEne;
  delete cMap;

  gStyle->SetPaintTextFormat();

  // html page footer
  htmlFile << "</body> " << endl;
  htmlFile << "</html> " << endl;

  htmlFile.close();

}

