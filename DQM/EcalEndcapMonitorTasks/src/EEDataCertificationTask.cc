#include <iostream>
#include <algorithm>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include <DataFormats/EcalDetId/interface/EEDetId.h>

#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
#include "CondFormats/RunInfo/interface/RunSummary.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"

#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include <DQM/EcalCommon/interface/Numbers.h>
#include "DQM/EcalCommon/interface/UtilsClient.h"

#include "DQM/EcalEndcapMonitorTasks/interface/EEDataCertificationTask.h"

using namespace cms;
using namespace edm;
using namespace std;

EEDataCertificationTask::EEDataCertificationTask(const ParameterSet& ps) {

  // cloneME switch
  cloneME_ = ps.getUntrackedParameter<bool>("cloneME", true);

  dqmStore_ = Service<DQMStore>().operator->();

  prefixME_ = ps.getUntrackedParameter<string>("prefixME", "");

  enableCleanup_ = ps.getUntrackedParameter<bool>("enableCleanup", false);

  mergeRuns_ = ps.getUntrackedParameter<bool>("mergeRuns", false);

  meEEDataCertificationSummary_ = 0;
  meEEDataCertificationSummaryMap_ = 0;
  for (int i = 0; i < 18; i++) {
    meEEDataCertification_[i] = 0;
  }

  hDQM_ = 0;
  hDAQ_ = 0;
  hIntegrityByLumi_ = 0;
  hFrontendByLumi_ = 0;

}

EEDataCertificationTask::~EEDataCertificationTask() {

}

void EEDataCertificationTask::beginJob(void){

  char histo[200];
  
  if ( dqmStore_ ) {

    dqmStore_->setCurrentFolder(prefixME_ + "/EventInfo");
    
    sprintf(histo, "CertificationSummary");
    meEEDataCertificationSummary_ = dqmStore_->bookFloat(histo);
    meEEDataCertificationSummary_->Fill(-1.0);

    sprintf(histo, "CertificationSummaryMap");
    meEEDataCertificationSummaryMap_ = dqmStore_->book2D(histo,histo, 200, 0., 200., 100, 0., 100.);
    meEEDataCertificationSummaryMap_->setAxisTitle("jx", 1);
    meEEDataCertificationSummaryMap_->setAxisTitle("jy", 2);

    dqmStore_->setCurrentFolder(prefixME_ + "/EventInfo/CertificationContents");

    for (int i = 0; i < 18; i++) {
      sprintf(histo, "EcalEndcap_%s", Numbers::sEE(i+1).c_str());
      meEEDataCertification_[i] = dqmStore_->bookFloat(histo);
      meEEDataCertification_[i]->Fill(-1.0);
    }

  }

}

void EEDataCertificationTask::endJob(void) {

  if ( enableCleanup_ ) this->cleanup();

}

void EEDataCertificationTask::beginLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const  edm::EventSetup& iSetup){

}

void EEDataCertificationTask::endLuminosityBlock(const edm::LuminosityBlock&  lumiBlock, const  edm::EventSetup& iSetup) {

  this->reset();

  char histo[200];

  MonitorElement* me;

  // evaluate the DQM quality of observables checked by lumi
  float DQMVal[18];
  for (int i = 0; i < 18; i++) {
    DQMVal[i] = 0.;
  }

  // evaluate the quality of observables checked by lumi
  sprintf(histo, (prefixME_ + "/EEIntegrityTask/EEIT weighted integrity errors by lumi").c_str());
  me = dqmStore_->get(histo);
  hIntegrityByLumi_ = UtilsClient::getHisto<TH1F*>( me, cloneME_, hIntegrityByLumi_ );

  sprintf(histo, (prefixME_ + "/EEStatusFlagsTask/FEStatus/EESFT weighted frontend errors by lumi").c_str());
  me = dqmStore_->get(histo);
  hFrontendByLumi_ = UtilsClient::getHisto<TH1F*>( me, cloneME_, hFrontendByLumi_ );
  
  float integrityErrSum, frontendErrSum;
  integrityErrSum = frontendErrSum = 0.;
  for ( int i=0; i<18; i++) {
    float ismIntegrityQual = 1.0;
    if( hIntegrityByLumi_ && hIntegrityByLumi_->GetBinContent(0) > 0 ) {
      float errors = hIntegrityByLumi_->GetBinContent(i+1);
      ismIntegrityQual = 1.0 - errors/hIntegrityByLumi_->GetBinContent(0);
      integrityErrSum += errors;
    }
    float ismFrontendQual = 1.0;
    if( hFrontendByLumi_ && hFrontendByLumi_->GetBinContent(0) > 0 ) {
      float errors = hFrontendByLumi_->GetBinContent(i+1);
      ismFrontendQual = 1.0 - errors/hFrontendByLumi_->GetBinContent(0);
      frontendErrSum += errors;
    }
    DQMVal[i] = min(ismIntegrityQual,ismFrontendQual);

    sprintf(histo, "EcalEndcap_%s", Numbers::sEE(i+1).c_str());
    me = dqmStore_->get(prefixME_ + "/EventInfo/reportSummaryContents/" + histo);
    if( me ) me->Fill(DQMVal[i]);

    sprintf(histo, "reportSummaryMap");
    me = dqmStore_->get(prefixME_ + "/EventInfo/" + histo );

    if( me ) {
      for ( int iz = -1; iz < 2; iz+=2 ) {
        for ( int ix = 1; ix <= 100; ix++ ) {
          for ( int iy = 1; iy <= 100; iy++ ) {
            int jx = (iz==1) ? 100 + ix : ix;
            int jy = iy;
            if( EEDetId::validDetId(ix, iy, iz) ) me->setBinContent( jx, jy, DQMVal[i] );
            else me->setBinContent( jx, jy, -1.0 );
          }
        }
      }
    }

  }

  float totDQMVal = 0.;
  
  float integrityQual = 1.0;
  if( hIntegrityByLumi_ && hIntegrityByLumi_->GetBinContent(0) > 0 ) integrityQual = 1.0 - integrityErrSum/hIntegrityByLumi_->GetBinContent(0);
  float frontendQual = 1.0;
  if( hFrontendByLumi_ && hFrontendByLumi_->GetBinContent(0) > 0 ) frontendQual = 1.0 - frontendErrSum/hFrontendByLumi_->GetBinContent(0);
  
  totDQMVal = min(integrityQual,frontendQual);

  sprintf(histo, (prefixME_ + "/EventInfo/reportSummary").c_str());
  me = dqmStore_->get(histo);
  if( me ) me->Fill(totDQMVal);

  // now combine reduced DQM with DCS and DAQ
  sprintf(histo, (prefixME_ + "/EventInfo/DAQSummaryMap").c_str());
  me = dqmStore_->get(histo);
  hDAQ_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, hDAQ_ );

  sprintf(histo, (prefixME_ + "/EventInfo/DCSSummaryMap").c_str());
  me = dqmStore_->get(histo);
  hDCS_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, hDCS_ );

  float sumCert = 0.;
  float sumCertEE[18];
  int nValidChannels = 0;
  int nValidChannelsEE[18];

  for (int i = 0; i < 18; i++) {
    sumCertEE[i] = 0;
    nValidChannelsEE[i] = 0;
  }

  for ( int iz = -1; iz < 2; iz+=2 ) {
    for ( int ix = 1; ix <= 100; ix++ ) {
      for ( int iy = 1; iy <= 100; iy++ ) {
        int jx = (iz==1) ? 100 + ix : ix;
        int jy = iy;
        if( EEDetId::validDetId(ix, iy, iz) ) {

          // map the 1-18 index to the correct SM
          int ism = 0;
          int firstSec = ( iz < 0 ) ? 1 : 10;
          int lastSec = ( iz < 0 ) ? 9 : 18;
          for (int i = firstSec; i <= lastSec; i++) {
            if ( Numbers::validEE(i, ix, iy) ) ism = i;
          }

          float xvalDQM = DQMVal[ism-1];

          float xvalDAQ, xvalDCS;
          xvalDAQ = xvalDCS = -1.;
          float xcert = -1;
          
          if ( hDAQ_ ) xvalDAQ = hDAQ_->GetBinContent( jx, jy );
          if ( hDCS_ ) xvalDCS = hDCS_->GetBinContent( jx, jy );

          // all white means problems: DAQ and DCS not available and DQM empty
          if ( xvalDQM == -1 && xvalDAQ == -1 && xvalDCS == -1) xcert = 0.0;
          else {
            // do not consider the white value of DAQ and DCS (problems with DB)
            xcert = fabs(xvalDQM) * fabs(xvalDAQ) * fabs(xvalDCS);
          }

          if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->setBinContent( jx, jy, xcert );

          sumCertEE[ism-1] += xcert;
          nValidChannelsEE[ism-1]++;

          sumCert += xcert;
          nValidChannels++;

        } else {
          if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->setBinContent( jx, jy, -1.0 );
        }
      }
    }
  }

  if( meEEDataCertificationSummary_ ) { 
    if( nValidChannels>0 ) meEEDataCertificationSummary_->Fill( sumCert/nValidChannels );
    else meEEDataCertificationSummary_->Fill( 0.0 );
  }

  for (int i = 0; i < 18; i++) {
    if( meEEDataCertification_[i] ) {
      if( nValidChannelsEE[i]>0 ) meEEDataCertification_[i]->Fill( sumCertEE[i]/nValidChannelsEE[i] );
      else meEEDataCertification_[i]->Fill( 0.0 );
    }
  }

}

void EEDataCertificationTask::beginRun(const Run& r, const EventSetup& c) {

  if ( ! mergeRuns_ ) this->reset();

}

void EEDataCertificationTask::endRun(const Run& r, const EventSetup& c) {

  this->reset();

  char histo[200];
  
  MonitorElement* me;

  sprintf(histo, (prefixME_ + "/EventInfo/reportSummaryMap").c_str());
  me = dqmStore_->get(histo);
  hDQM_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, hDQM_ );

  sprintf(histo, (prefixME_ + "/EventInfo/DAQSummaryMap").c_str());
  me = dqmStore_->get(histo);
  hDAQ_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, hDAQ_ );

  sprintf(histo, (prefixME_ + "/EventInfo/DCSSummaryMap").c_str());
  me = dqmStore_->get(histo);
  hDCS_ = UtilsClient::getHisto<TH2F*>( me, cloneME_, hDCS_ );

  float sumCert = 0.;
  float sumCertEE[18];
  int nValidChannels = 0;
  int nValidChannelsEE[18];

  for (int i = 0; i < 18; i++) {
    sumCertEE[i] = 0;
    nValidChannelsEE[i] = 0;
  }

  for ( int iz = -1; iz < 2; iz+=2 ) {
    for ( int ix = 1; ix <= 100; ix++ ) {
      for ( int iy = 1; iy <= 100; iy++ ) {
        int jx = (iz==1) ? 100 + ix : ix;
        int jy = iy;
        if( EEDetId::validDetId(ix, iy, iz) ) {

          float xvalDQM, xvalDAQ, xvalDCS;
          xvalDQM = xvalDAQ = xvalDCS = -1.;
          float xcert = -1;
          
          if ( hDQM_ ) xvalDQM = hDQM_->GetBinContent( jx, jy );
          if ( hDAQ_ ) xvalDAQ = hDAQ_->GetBinContent( jx, jy );
          if ( hDCS_ ) xvalDCS = hDCS_->GetBinContent( jx, jy );

          // all white means problems: DAQ and DCS not available and DQM empty
          if ( xvalDQM == -1 && xvalDAQ == -1 && xvalDCS == -1) xcert = 0.0;
          else {
            // do not consider the white value of DAQ and DCS (problems with DB)
            xcert = fabs(xvalDQM) * fabs(xvalDAQ) * fabs(xvalDCS);
          }

          if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->setBinContent( jx, jy, xcert );

          // map the 1-18 index to the correct SM
          int firstSec = ( iz < 0 ) ? 1 : 10;
          int lastSec = ( iz < 0 ) ? 9 : 18;
          for (int i = firstSec; i <= lastSec; i++) {
            if ( Numbers::validEE(i, ix, iy) ) {
              sumCertEE[i-1] += xcert;
              nValidChannelsEE[i-1]++;
            }
          }

          sumCert += xcert;
          nValidChannels++;

        } else {
          if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->setBinContent( jx, jy, -1.0 );
        }
      }
    }
  }

  if( meEEDataCertificationSummary_ ) { 
    if( nValidChannels>0 ) meEEDataCertificationSummary_->Fill( sumCert/nValidChannels );
    else meEEDataCertificationSummary_->Fill( 0.0 );
  }

  for (int i = 0; i < 18; i++) {
    if( meEEDataCertification_[i] ) {
      if( nValidChannelsEE[i]>0 ) meEEDataCertification_[i]->Fill( sumCertEE[i]/nValidChannelsEE[i] );
      else meEEDataCertification_[i]->Fill( 0.0 );
    }
  }

}

void EEDataCertificationTask::reset(void) {

  if ( meEEDataCertificationSummary_ ) meEEDataCertificationSummary_->Reset();

  for (int i = 0; i < 18; i++) {
    if ( meEEDataCertification_[i] ) meEEDataCertification_[i]->Reset();
  }

  if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->Reset();

  if ( meEEDataCertificationSummaryMap_ ) meEEDataCertificationSummaryMap_->Reset();
  
}


void EEDataCertificationTask::cleanup(void){
  
  if ( cloneME_ ) {
    if( hDQM_ ) delete hDQM_;
    if( hDAQ_ ) delete hDAQ_;
    if( hDCS_ ) delete hDCS_;
    if( hIntegrityByLumi_ ) delete hIntegrityByLumi_;
    if( hFrontendByLumi_ ) delete hFrontendByLumi_;
  }
  hDQM_ = 0;
  hDAQ_ = 0;
  hDCS_ = 0;
  hIntegrityByLumi_ = 0;
  hFrontendByLumi_ = 0;

  if ( dqmStore_ ) {
    dqmStore_->setCurrentFolder(prefixME_ + "/EventInfo");
    if ( meEEDataCertificationSummary_ ) dqmStore_->removeElement( meEEDataCertificationSummary_->getName() );
    if ( meEEDataCertificationSummaryMap_ ) dqmStore_->removeElement( meEEDataCertificationSummaryMap_->getName() );

    dqmStore_->setCurrentFolder(prefixME_ + "/EventInfo/CertificationContents");
    for (int i = 0; i < 18; i++) {
      if ( meEEDataCertification_[i] ) dqmStore_->removeElement( meEEDataCertification_[i]->getName() );
    }    
  }

}

void EEDataCertificationTask::analyze(const Event& e, const EventSetup& c){ 

}
