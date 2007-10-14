/*
 * \file EELaserClient.cc
 *
 * $Date: 2007/10/14 12:31:59 $
 * $Revision: 1.35 $
 * \author G. Della Ricca
 * \author G. Franzoni
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
#include "DQMServices/Core/interface/QTestStatus.h"
#include "DQMServices/QualityTests/interface/QCriterionRoot.h"

#include "OnlineDB/EcalCondDB/interface/RunTag.h"
#include "OnlineDB/EcalCondDB/interface/RunIOV.h"
#include "OnlineDB/EcalCondDB/interface/MonLaserBlueDat.h"
#include "OnlineDB/EcalCondDB/interface/MonLaserGreenDat.h"
#include "OnlineDB/EcalCondDB/interface/MonLaserIRedDat.h"
#include "OnlineDB/EcalCondDB/interface/MonLaserRedDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPNBlueDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPNGreenDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPNIRedDat.h"
#include "OnlineDB/EcalCondDB/interface/MonPNRedDat.h"
#include "OnlineDB/EcalCondDB/interface/RunCrystalErrorsDat.h"
#include "OnlineDB/EcalCondDB/interface/RunPNErrorsDat.h"

#include "CondTools/Ecal/interface/EcalErrorDictionary.h"

#include "DQM/EcalCommon/interface/EcalErrorMask.h"
#include <DQM/EcalCommon/interface/UtilsClient.h>
#include <DQM/EcalCommon/interface/LogicID.h>
#include <DQM/EcalCommon/interface/Numbers.h>

#include <DQM/EcalEndcapMonitorClient/interface/EELaserClient.h>

using namespace cms;
using namespace edm;
using namespace std;

EELaserClient::EELaserClient(const ParameterSet& ps){

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

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    h01_[ism-1] = 0;
    h02_[ism-1] = 0;
    h03_[ism-1] = 0;
    h04_[ism-1] = 0;
    h05_[ism-1] = 0;
    h06_[ism-1] = 0;
    h07_[ism-1] = 0;
    h08_[ism-1] = 0;

    h09_[ism-1] = 0;
    h10_[ism-1] = 0;
    h11_[ism-1] = 0;
    h12_[ism-1] = 0;

    h13_[ism-1] = 0;
    h14_[ism-1] = 0;
    h15_[ism-1] = 0;
    h16_[ism-1] = 0;
    h17_[ism-1] = 0;
    h18_[ism-1] = 0;
    h19_[ism-1] = 0;
    h20_[ism-1] = 0;

    h21_[ism-1] = 0;
    h22_[ism-1] = 0;
    h23_[ism-1] = 0;
    h24_[ism-1] = 0;

    hs01_[ism-1] = 0;
    hs02_[ism-1] = 0;
    hs03_[ism-1] = 0;
    hs04_[ism-1] = 0;

    hs05_[ism-1] = 0;
    hs06_[ism-1] = 0;
    hs07_[ism-1] = 0;
    hs08_[ism-1] = 0;

    i01_[ism-1] = 0;
    i02_[ism-1] = 0;
    i03_[ism-1] = 0;
    i04_[ism-1] = 0;
    i05_[ism-1] = 0;
    i06_[ism-1] = 0;
    i07_[ism-1] = 0;
    i08_[ism-1] = 0;

    i09_[ism-1] = 0;
    i10_[ism-1] = 0;
    i11_[ism-1] = 0;
    i12_[ism-1] = 0;
    i13_[ism-1] = 0;
    i14_[ism-1] = 0;
    i15_[ism-1] = 0;
    i16_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    meg01_[ism-1] = 0;
    meg02_[ism-1] = 0;
    meg03_[ism-1] = 0;
    meg04_[ism-1] = 0;

    meg05_[ism-1] = 0;
    meg06_[ism-1] = 0;
    meg07_[ism-1] = 0;
    meg08_[ism-1] = 0;
    meg09_[ism-1] = 0;
    meg10_[ism-1] = 0;
    meg11_[ism-1] = 0;
    meg12_[ism-1] = 0;

    mea01_[ism-1] = 0;
    mea02_[ism-1] = 0;
    mea03_[ism-1] = 0;
    mea04_[ism-1] = 0;
    mea05_[ism-1] = 0;
    mea06_[ism-1] = 0;
    mea07_[ism-1] = 0;
    mea08_[ism-1] = 0;

    met01_[ism-1] = 0;
    met02_[ism-1] = 0;
    met03_[ism-1] = 0;
    met04_[ism-1] = 0;
    met05_[ism-1] = 0;
    met06_[ism-1] = 0;
    met07_[ism-1] = 0;
    met08_[ism-1] = 0;

    metav01_[ism-1] = 0;
    metav02_[ism-1] = 0;
    metav03_[ism-1] = 0;
    metav04_[ism-1] = 0;
    metav05_[ism-1] = 0;
    metav06_[ism-1] = 0;
    metav07_[ism-1] = 0;
    metav08_[ism-1] = 0;

    metrms01_[ism-1] = 0;
    metrms02_[ism-1] = 0;
    metrms03_[ism-1] = 0;
    metrms04_[ism-1] = 0;
    metrms05_[ism-1] = 0;
    metrms06_[ism-1] = 0;
    metrms07_[ism-1] = 0;
    metrms08_[ism-1] = 0;

    meaopn01_[ism-1] = 0;
    meaopn02_[ism-1] = 0;
    meaopn03_[ism-1] = 0;
    meaopn04_[ism-1] = 0;
    meaopn05_[ism-1] = 0;
    meaopn06_[ism-1] = 0;
    meaopn07_[ism-1] = 0;
    meaopn08_[ism-1] = 0;
    
    mepnprms01_[ism-1] = 0;
    mepnprms02_[ism-1] = 0;
    mepnprms03_[ism-1] = 0;
    mepnprms04_[ism-1] = 0;
    mepnprms05_[ism-1] = 0;
    mepnprms06_[ism-1] = 0;
    mepnprms07_[ism-1] = 0;
    mepnprms08_[ism-1] = 0;
    
    qth01_[ism-1] = 0;
    qth02_[ism-1] = 0;
    qth03_[ism-1] = 0;
    qth04_[ism-1] = 0;
    qth05_[ism-1] = 0;
    qth06_[ism-1] = 0;
    qth07_[ism-1] = 0;
    qth08_[ism-1] = 0;

    qth09_[ism-1] = 0;
    qth10_[ism-1] = 0;
    qth11_[ism-1] = 0;
    qth12_[ism-1] = 0;
    qth13_[ism-1] = 0;
    qth14_[ism-1] = 0;
    qth15_[ism-1] = 0;
    qth16_[ism-1] = 0;
    qth17_[ism-1] = 0;
    qth18_[ism-1] = 0;
    qth19_[ism-1] = 0;
    qth20_[ism-1] = 0;
    qth21_[ism-1] = 0;
    qth22_[ism-1] = 0;
    qth23_[ism-1] = 0;
    qth24_[ism-1] = 0;

    qtg01_[ism-1] = 0;
    qtg02_[ism-1] = 0;
    qtg03_[ism-1] = 0;
    qtg04_[ism-1] = 0;

    qtg05_[ism-1] = 0;
    qtg06_[ism-1] = 0;
    qtg07_[ism-1] = 0;
    qtg08_[ism-1] = 0;
    qtg09_[ism-1] = 0;
    qtg10_[ism-1] = 0;
    qtg11_[ism-1] = 0;
    qtg12_[ism-1] = 0;

  }

  percentVariation_ = 0.4;

  amplitudeThresholdPnG01_ = 50.;
  amplitudeThresholdPnG16_ = 50.;
  
  pedPnExpectedMean_[0] = 750.0;
  pedPnExpectedMean_[1] = 750.0;
  
  pedPnDiscrepancyMean_[0] = 100.0;
  pedPnDiscrepancyMean_[1] = 100.0;
  
  pedPnRMSThreshold_[0] = 1.0; // value at h4; expected nominal: 0.5
  pedPnRMSThreshold_[1] = 3.0; // value at h4; expected nominal: 1.6
  
}

EELaserClient::~EELaserClient(){

}

void EELaserClient::beginJob(MonitorUserInterface* mui){

  mui_ = mui;
  dbe_ = mui->getBEInterface();

  if ( verbose_ ) cout << "EELaserClient: beginJob" << endl;

  ievt_ = 0;
  jevt_ = 0;

  if ( enableQT_ ) {

    Char_t qtname[200];

    for ( unsigned int i=0; i<superModules_.size(); i++ ) {

      int ism = superModules_[i];

      sprintf(qtname, "EELT laser quality %s L1A", Numbers::sEE(ism).c_str());
      qth01_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L2A", Numbers::sEE(ism).c_str());
      qth02_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L3A", Numbers::sEE(ism).c_str());
      qth03_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L4A", Numbers::sEE(ism).c_str());
      qth04_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L1B", Numbers::sEE(ism).c_str());
      qth05_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L2B", Numbers::sEE(ism).c_str());
      qth06_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L3B", Numbers::sEE(ism).c_str());
      qth07_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser quality %s L4B", Numbers::sEE(ism).c_str());
      qth08_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L1 G01", Numbers::sEE(ism).c_str());
      qth09_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L2 G01", Numbers::sEE(ism).c_str());
      qth10_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L3 G01", Numbers::sEE(ism).c_str());
      qth11_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L4 G01", Numbers::sEE(ism).c_str());
      qth12_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L1 G01", Numbers::sEE(ism).c_str());
      qth13_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L2 G01", Numbers::sEE(ism).c_str());
      qth14_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L3 G01", Numbers::sEE(ism).c_str());
      qth15_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L4 G01", Numbers::sEE(ism).c_str());
      qth16_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L1 G16", Numbers::sEE(ism).c_str());
      qth17_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L2 G16", Numbers::sEE(ism).c_str());
      qth18_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L3 G16", Numbers::sEE(ism).c_str());
      qth19_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser amplitude quality PNs %s L4 G16", Numbers::sEE(ism).c_str());
      qth20_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L1 G16", Numbers::sEE(ism).c_str());
      qth21_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L2 G16", Numbers::sEE(ism).c_str());
      qth22_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L3 G16", Numbers::sEE(ism).c_str());
      qth23_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT laser pedestal quality PNs %s L4 G16", Numbers::sEE(ism).c_str());
      qth24_[ism-1] = dynamic_cast<MEContentsProf2DWithinRangeROOT*> (dbe_->createQTest(ContentsProf2DWithinRangeROOT::getAlgoName(), qtname));

      qth01_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth02_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth03_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth04_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth05_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth06_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth07_[ism-1]->setMeanRange(100.0, 4096.0*12.);
      qth08_[ism-1]->setMeanRange(100.0, 4096.0*12.);

      qth09_[ism-1]->setMeanRange(amplitudeThresholdPnG01_, 4096.0);
      qth10_[ism-1]->setMeanRange(amplitudeThresholdPnG01_, 4096.0);
      qth11_[ism-1]->setMeanRange(amplitudeThresholdPnG01_, 4096.0);
      qth12_[ism-1]->setMeanRange(amplitudeThresholdPnG01_, 4096.0);
      qth13_[ism-1]->setMeanRange(pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0],
				  pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0]);
      qth14_[ism-1]->setMeanRange(pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0],
				  pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0]);
      qth15_[ism-1]->setMeanRange(pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0],
				  pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0]);
      qth16_[ism-1]->setMeanRange(pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0],
				  pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0]);
      qth17_[ism-1]->setMeanRange(amplitudeThresholdPnG16_, 4096.0);
      qth18_[ism-1]->setMeanRange(amplitudeThresholdPnG16_, 4096.0);
      qth19_[ism-1]->setMeanRange(amplitudeThresholdPnG16_, 4096.0);
      qth20_[ism-1]->setMeanRange(amplitudeThresholdPnG16_, 4096.0);
      qth21_[ism-1]->setMeanRange(pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1],
				  pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1]);
      qth22_[ism-1]->setMeanRange(pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1],
				  pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1]);
      qth23_[ism-1]->setMeanRange(pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1],
				  pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1]);
      qth24_[ism-1]->setMeanRange(pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1],
				  pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1]);
      qth01_[ism-1]->setMeanTolerance(percentVariation_);
      qth02_[ism-1]->setMeanTolerance(percentVariation_);
      qth03_[ism-1]->setMeanTolerance(percentVariation_);
      qth04_[ism-1]->setMeanTolerance(percentVariation_);
      qth05_[ism-1]->setMeanTolerance(percentVariation_);
      qth06_[ism-1]->setMeanTolerance(percentVariation_);
      qth07_[ism-1]->setMeanTolerance(percentVariation_);
      qth08_[ism-1]->setMeanTolerance(percentVariation_);

      qth09_[ism-1]->setRMSRange(0.0, 4096.0);
      qth10_[ism-1]->setRMSRange(0.0, 4096.0);
      qth11_[ism-1]->setRMSRange(0.0, 4096.0);
      qth12_[ism-1]->setRMSRange(0.0, 4096.0);
      qth13_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[0]);
      qth14_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[0]);
      qth15_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[0]);
      qth16_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[0]);
      qth17_[ism-1]->setRMSRange(0.0, 4096.0);
      qth18_[ism-1]->setRMSRange(0.0, 4096.0);
      qth19_[ism-1]->setRMSRange(0.0, 4096.0);
      qth20_[ism-1]->setRMSRange(0.0, 4096.0);
      qth21_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[1]);
      qth22_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[1]);
      qth23_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[1]);
      qth24_[ism-1]->setRMSRange(0.0, pedPnRMSThreshold_[1]);
      qth01_[ism-1]->setMinimumEntries(10*1700);
      qth02_[ism-1]->setMinimumEntries(10*1700);
      qth03_[ism-1]->setMinimumEntries(10*1700);
      qth04_[ism-1]->setMinimumEntries(10*1700);
      qth05_[ism-1]->setMinimumEntries(10*1700);
      qth06_[ism-1]->setMinimumEntries(10*1700);
      qth07_[ism-1]->setMinimumEntries(10*1700);
      qth08_[ism-1]->setMinimumEntries(10*1700);

      qth09_[ism-1]->setMinimumEntries(10*10);
      qth10_[ism-1]->setMinimumEntries(10*10);
      qth11_[ism-1]->setMinimumEntries(10*10);
      qth12_[ism-1]->setMinimumEntries(10*10);
      qth13_[ism-1]->setMinimumEntries(10*10);
      qth14_[ism-1]->setMinimumEntries(10*10);
      qth15_[ism-1]->setMinimumEntries(10*10);
      qth16_[ism-1]->setMinimumEntries(10*10);
      qth17_[ism-1]->setMinimumEntries(10*10);
      qth18_[ism-1]->setMinimumEntries(10*10);
      qth19_[ism-1]->setMinimumEntries(10*10);
      qth20_[ism-1]->setMinimumEntries(10*10);
      qth21_[ism-1]->setMinimumEntries(10*10);
      qth22_[ism-1]->setMinimumEntries(10*10);
      qth23_[ism-1]->setMinimumEntries(10*10);
      qth24_[ism-1]->setMinimumEntries(10*10);

      qth01_[ism-1]->setErrorProb(1.00);
      qth02_[ism-1]->setErrorProb(1.00);
      qth03_[ism-1]->setErrorProb(1.00);
      qth04_[ism-1]->setErrorProb(1.00);
      qth05_[ism-1]->setErrorProb(1.00);
      qth06_[ism-1]->setErrorProb(1.00);
      qth07_[ism-1]->setErrorProb(1.00);
      qth08_[ism-1]->setErrorProb(1.00);

      qth09_[ism-1]->setErrorProb(1.00);
      qth10_[ism-1]->setErrorProb(1.00);
      qth11_[ism-1]->setErrorProb(1.00);
      qth12_[ism-1]->setErrorProb(1.00);
      qth13_[ism-1]->setErrorProb(1.00);
      qth14_[ism-1]->setErrorProb(1.00);
      qth15_[ism-1]->setErrorProb(1.00);
      qth16_[ism-1]->setErrorProb(1.00);
      qth17_[ism-1]->setErrorProb(1.00);
      qth18_[ism-1]->setErrorProb(1.00);
      qth19_[ism-1]->setErrorProb(1.00);
      qth20_[ism-1]->setErrorProb(1.00);
      qth21_[ism-1]->setErrorProb(1.00);
      qth22_[ism-1]->setErrorProb(1.00);
      qth23_[ism-1]->setErrorProb(1.00);
      qth24_[ism-1]->setErrorProb(1.00);

      sprintf(qtname, "EELT quality test %s L1", Numbers::sEE(ism).c_str());
      qtg01_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test %s L2", Numbers::sEE(ism).c_str());
      qtg02_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test %s L3", Numbers::sEE(ism).c_str());
      qtg03_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test %s L4", Numbers::sEE(ism).c_str());
      qtg04_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      qtg01_[ism-1]->setMeanRange(1., 6.);
      qtg02_[ism-1]->setMeanRange(1., 6.);
      qtg03_[ism-1]->setMeanRange(1., 6.);
      qtg04_[ism-1]->setMeanRange(1., 6.);

      qtg01_[ism-1]->setErrorProb(1.00);
      qtg02_[ism-1]->setErrorProb(1.00);
      qtg03_[ism-1]->setErrorProb(1.00);
      qtg04_[ism-1]->setErrorProb(1.00);

      sprintf(qtname, "EELT quality test PNs %s L1 G01", Numbers::sEE(ism).c_str());
      qtg05_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L2 G01", Numbers::sEE(ism).c_str());
      qtg06_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L3 G01", Numbers::sEE(ism).c_str());
      qtg07_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L4 G01", Numbers::sEE(ism).c_str());
      qtg08_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L1 G16", Numbers::sEE(ism).c_str());
      qtg09_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L2 G16", Numbers::sEE(ism).c_str());
      qtg10_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L3 G16", Numbers::sEE(ism).c_str());
      qtg11_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      sprintf(qtname, "EELT quality test PNs %s L4 G16", Numbers::sEE(ism).c_str());
      qtg12_[ism-1] = dynamic_cast<MEContentsTH2FWithinRangeROOT*> (dbe_->createQTest(ContentsTH2FWithinRangeROOT::getAlgoName(), qtname));

      qtg05_[ism-1]->setMeanRange(1., 6.);
      qtg06_[ism-1]->setMeanRange(1., 6.);
      qtg07_[ism-1]->setMeanRange(1., 6.);
      qtg08_[ism-1]->setMeanRange(1., 6.);
      qtg09_[ism-1]->setMeanRange(1., 6.);
      qtg10_[ism-1]->setMeanRange(1., 6.);
      qtg11_[ism-1]->setMeanRange(1., 6.);
      qtg12_[ism-1]->setMeanRange(1., 6.);

      qtg05_[ism-1]->setErrorProb(1.00);
      qtg06_[ism-1]->setErrorProb(1.00);
      qtg07_[ism-1]->setErrorProb(1.00);
      qtg08_[ism-1]->setErrorProb(1.00);
      qtg09_[ism-1]->setErrorProb(1.00);
      qtg10_[ism-1]->setErrorProb(1.00);
      qtg11_[ism-1]->setErrorProb(1.00);
      qtg12_[ism-1]->setErrorProb(1.00);

    }

  }

}

void EELaserClient::beginRun(void){

  if ( verbose_ ) cout << "EELaserClient: beginRun" << endl;

  jevt_ = 0;

  this->setup();

  this->subscribe();

}

void EELaserClient::endJob(void) {

  if ( verbose_ ) cout << "EELaserClient: endJob, ievt = " << ievt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EELaserClient::endRun(void) {

  if ( verbose_ ) cout << "EELaserClient: endRun, jevt = " << jevt_ << endl;

  this->unsubscribe();

  this->cleanup();

}

void EELaserClient::setup(void) {

  Char_t histo[200];

  dbe_->setCurrentFolder( "EcalEndcap/EELaserClient" );

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( meg01_[ism-1] ) dbe_->removeElement( meg01_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L1 %s", Numbers::sEE(ism).c_str());
    meg01_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    if ( meg02_[ism-1] ) dbe_->removeElement( meg02_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L2 %s", Numbers::sEE(ism).c_str());
    meg02_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    if ( meg03_[ism-1] ) dbe_->removeElement( meg03_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L3 %s", Numbers::sEE(ism).c_str());
    meg03_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);
    if ( meg04_[ism-1] ) dbe_->removeElement( meg04_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L4 %s", Numbers::sEE(ism).c_str());
    meg04_[ism-1] = dbe_->book2D(histo, histo, 50, Numbers::ix0EE(ism)+0., Numbers::ix0EE(ism)+50., 50, Numbers::iy0EE(ism)+0., Numbers::iy0EE(ism)+50.);

    if ( meg05_[ism-1] ) dbe_->removeElement( meg05_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L1 PNs %s G01", Numbers::sEE(ism).c_str());
    meg05_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg06_[ism-1] ) dbe_->removeElement( meg06_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L2 PNs %s G01", Numbers::sEE(ism).c_str());
    meg06_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg07_[ism-1] ) dbe_->removeElement( meg07_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L3 PNs %s G01", Numbers::sEE(ism).c_str());
    meg07_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg08_[ism-1] ) dbe_->removeElement( meg08_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L4 PNs %s G01", Numbers::sEE(ism).c_str());
    meg08_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg09_[ism-1] ) dbe_->removeElement( meg09_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L1 PNs %s G16", Numbers::sEE(ism).c_str());
    meg09_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg10_[ism-1] ) dbe_->removeElement( meg10_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L2 PNs %s G16", Numbers::sEE(ism).c_str());
    meg10_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg11_[ism-1] ) dbe_->removeElement( meg11_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L3 PNs %s G16", Numbers::sEE(ism).c_str());
    meg11_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);
    if ( meg12_[ism-1] ) dbe_->removeElement( meg12_[ism-1]->getName() );
    sprintf(histo, "EELT laser quality L4 PNs %s G16", Numbers::sEE(ism).c_str());
    meg12_[ism-1] = dbe_->book2D(histo, histo, 10, 0., 10., 1, 0., 5.);

    if ( mea01_[ism-1] ) dbe_->removeElement( mea01_[ism-1]->getName() );;
    sprintf(histo, "EELT amplitude L1A %s", Numbers::sEE(ism).c_str());
    mea01_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea02_[ism-1] ) dbe_->removeElement( mea02_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L2A %s", Numbers::sEE(ism).c_str());
    mea02_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea03_[ism-1] ) dbe_->removeElement( mea03_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L3A %s", Numbers::sEE(ism).c_str());
    mea03_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea04_[ism-1] ) dbe_->removeElement( mea04_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L4A %s", Numbers::sEE(ism).c_str());
    mea04_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea05_[ism-1] ) dbe_->removeElement( mea05_[ism-1]->getName() );;
    sprintf(histo, "EELT amplitude L1B %s", Numbers::sEE(ism).c_str());
    mea05_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea06_[ism-1] ) dbe_->removeElement( mea06_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L2B %s", Numbers::sEE(ism).c_str());
    mea06_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea07_[ism-1] ) dbe_->removeElement( mea07_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L3B %s", Numbers::sEE(ism).c_str());
    mea07_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( mea08_[ism-1] ) dbe_->removeElement( mea08_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude L4B %s", Numbers::sEE(ism).c_str());
    mea08_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);

    if ( met01_[ism-1] ) dbe_->removeElement( met01_[ism-1]->getName() );
    sprintf(histo, "EELT timing L1A %s", Numbers::sEE(ism).c_str());
    met01_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met02_[ism-1] ) dbe_->removeElement( met02_[ism-1]->getName() );
    sprintf(histo, "EELT timing L2A %s", Numbers::sEE(ism).c_str());
    met02_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met03_[ism-1] ) dbe_->removeElement( met03_[ism-1]->getName() );
    sprintf(histo, "EELT timing L3A %s", Numbers::sEE(ism).c_str());
    met03_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met04_[ism-1] ) dbe_->removeElement( met04_[ism-1]->getName() );
    sprintf(histo, "EELT timing L4A %s", Numbers::sEE(ism).c_str());
    met04_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met05_[ism-1] ) dbe_->removeElement( met05_[ism-1]->getName() );
    sprintf(histo, "EELT timing L1B %s", Numbers::sEE(ism).c_str());
    met05_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met06_[ism-1] ) dbe_->removeElement( met06_[ism-1]->getName() );
    sprintf(histo, "EELT timing L2B %s", Numbers::sEE(ism).c_str());
    met06_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met07_[ism-1] ) dbe_->removeElement( met07_[ism-1]->getName() );
    sprintf(histo, "EELT timing L3B %s", Numbers::sEE(ism).c_str());
    met07_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( met08_[ism-1] ) dbe_->removeElement( met08_[ism-1]->getName() );
    sprintf(histo, "EELT timing L4B %s", Numbers::sEE(ism).c_str());
    met08_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);

    if ( metav01_[ism-1] ) dbe_->removeElement( metav01_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L1A %s", Numbers::sEE(ism).c_str());
    metav01_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav02_[ism-1] ) dbe_->removeElement( metav02_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L2A %s", Numbers::sEE(ism).c_str());
    metav02_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav03_[ism-1] ) dbe_->removeElement( metav03_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L3A %s", Numbers::sEE(ism).c_str());
    metav03_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav04_[ism-1] ) dbe_->removeElement( metav04_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L4A %s", Numbers::sEE(ism).c_str());
    metav04_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav05_[ism-1] ) dbe_->removeElement( metav05_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L1B %s", Numbers::sEE(ism).c_str());
    metav05_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav06_[ism-1] ) dbe_->removeElement( metav06_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L2B %s", Numbers::sEE(ism).c_str());
    metav06_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav07_[ism-1] ) dbe_->removeElement( metav07_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L3B %s", Numbers::sEE(ism).c_str());
    metav07_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( metav08_[ism-1] ) dbe_->removeElement( metav08_[ism-1]->getName() );
    sprintf(histo, "EELT timing mean L4B %s", Numbers::sEE(ism).c_str());
    metav08_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);

    if ( metrms01_[ism-1] ) dbe_->removeElement( metrms01_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L1A %s", Numbers::sEE(ism).c_str());
    metrms01_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms02_[ism-1] ) dbe_->removeElement( metrms02_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L2A %s", Numbers::sEE(ism).c_str());
    metrms02_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms03_[ism-1] ) dbe_->removeElement( metrms03_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L3A %s", Numbers::sEE(ism).c_str());
    metrms03_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms04_[ism-1] ) dbe_->removeElement( metrms04_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L4A %s", Numbers::sEE(ism).c_str());
    metrms04_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms05_[ism-1] ) dbe_->removeElement( metrms05_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L1B %s", Numbers::sEE(ism).c_str());
    metrms05_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms06_[ism-1] ) dbe_->removeElement( metrms06_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L2B %s", Numbers::sEE(ism).c_str());
    metrms06_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms07_[ism-1] ) dbe_->removeElement( metrms07_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L3B %s", Numbers::sEE(ism).c_str());
    metrms07_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);
    if ( metrms08_[ism-1] ) dbe_->removeElement( metrms08_[ism-1]->getName() );
    sprintf(histo, "EELT timing rms L4B %s", Numbers::sEE(ism).c_str());
    metrms08_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 0.5);

    if ( meaopn01_[ism-1] ) dbe_->removeElement( meaopn01_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L1A %s", Numbers::sEE(ism).c_str());
    meaopn01_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn02_[ism-1] ) dbe_->removeElement( meaopn02_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L2A %s", Numbers::sEE(ism).c_str());
    meaopn02_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn03_[ism-1] ) dbe_->removeElement( meaopn03_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L3A %s", Numbers::sEE(ism).c_str());
    meaopn03_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn04_[ism-1] ) dbe_->removeElement( meaopn04_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L4A %s", Numbers::sEE(ism).c_str());
    meaopn04_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn05_[ism-1] ) dbe_->removeElement( meaopn05_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L1B %s", Numbers::sEE(ism).c_str());
    meaopn05_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn06_[ism-1] ) dbe_->removeElement( meaopn06_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L2B %s", Numbers::sEE(ism).c_str());
    meaopn06_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn07_[ism-1] ) dbe_->removeElement( meaopn07_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L3B %s", Numbers::sEE(ism).c_str());
    meaopn07_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    if ( meaopn08_[ism-1] ) dbe_->removeElement( meaopn08_[ism-1]->getName() );
    sprintf(histo, "EELT amplitude over PN L4B %s", Numbers::sEE(ism).c_str());
    meaopn08_[ism-1] = dbe_->book1D(histo, histo, 1700, 0., 1700.);
    
    if ( mepnprms01_[ism-1] ) dbe_->removeElement( mepnprms01_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G01 L1", Numbers::sEE(ism).c_str());
    mepnprms01_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms02_[ism-1] ) dbe_->removeElement( mepnprms02_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G01 L2", Numbers::sEE(ism).c_str());
    mepnprms02_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms03_[ism-1] ) dbe_->removeElement( mepnprms03_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G01 L3", Numbers::sEE(ism).c_str());
    mepnprms03_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms04_[ism-1] ) dbe_->removeElement( mepnprms04_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G01 L4", Numbers::sEE(ism).c_str());
    mepnprms04_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms05_[ism-1] ) dbe_->removeElement( mepnprms05_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G16 L1", Numbers::sEE(ism).c_str());
    mepnprms05_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms06_[ism-1] ) dbe_->removeElement( mepnprms06_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G16 L2", Numbers::sEE(ism).c_str());
    mepnprms06_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms07_[ism-1] ) dbe_->removeElement( mepnprms07_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G16 L3", Numbers::sEE(ism).c_str());
    mepnprms07_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);
    if ( mepnprms08_[ism-1] ) dbe_->removeElement( mepnprms08_[ism-1]->getName() );
    sprintf(histo, "EEPDT PNs pedestal rms %s G16 L4", Numbers::sEE(ism).c_str());
    mepnprms08_[ism-1] = dbe_->book1D(histo, histo, 100, 0., 10.);

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    UtilsClient::resetHisto( meg01_[ism-1] );
    UtilsClient::resetHisto( meg02_[ism-1] );
    UtilsClient::resetHisto( meg03_[ism-1] );
    UtilsClient::resetHisto( meg04_[ism-1] );

    UtilsClient::resetHisto( meg05_[ism-1] );
    UtilsClient::resetHisto( meg06_[ism-1] );
    UtilsClient::resetHisto( meg07_[ism-1] );
    UtilsClient::resetHisto( meg08_[ism-1] );
    UtilsClient::resetHisto( meg09_[ism-1] );
    UtilsClient::resetHisto( meg10_[ism-1] );
    UtilsClient::resetHisto( meg11_[ism-1] );
    UtilsClient::resetHisto( meg12_[ism-1] );

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        meg01_[ism-1]->setBinContent( ix, iy, -1. );
        meg02_[ism-1]->setBinContent( ix, iy, -1. );
        meg03_[ism-1]->setBinContent( ix, iy, -1. );
        meg04_[ism-1]->setBinContent( ix, iy, -1. );

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          meg01_[ism-1]->setBinContent( ix, iy, 2. );
          meg02_[ism-1]->setBinContent( ix, iy, 2. );
          meg03_[ism-1]->setBinContent( ix, iy, 2. );
          meg04_[ism-1]->setBinContent( ix, iy, 2. );
        }

      }
    }

    for ( int i = 1; i <= 10; i++ ) {

        meg05_[ism-1]->setBinContent( i, 1, 2. );
        meg06_[ism-1]->setBinContent( i, 1, 2. );
        meg07_[ism-1]->setBinContent( i, 1, 2. );
        meg08_[ism-1]->setBinContent( i, 1, 2. );
        meg09_[ism-1]->setBinContent( i, 1, 2. );
        meg10_[ism-1]->setBinContent( i, 1, 2. );
        meg11_[ism-1]->setBinContent( i, 1, 2. );
        meg12_[ism-1]->setBinContent( i, 1, 2. );

    }

    UtilsClient::resetHisto( mea01_[ism-1] );
    UtilsClient::resetHisto( mea02_[ism-1] );
    UtilsClient::resetHisto( mea03_[ism-1] );
    UtilsClient::resetHisto( mea04_[ism-1] );
    UtilsClient::resetHisto( mea05_[ism-1] );
    UtilsClient::resetHisto( mea06_[ism-1] );
    UtilsClient::resetHisto( mea07_[ism-1] );
    UtilsClient::resetHisto( mea08_[ism-1] );

    UtilsClient::resetHisto( met01_[ism-1] );
    UtilsClient::resetHisto( met02_[ism-1] );
    UtilsClient::resetHisto( met03_[ism-1] );
    UtilsClient::resetHisto( met04_[ism-1] );
    UtilsClient::resetHisto( met05_[ism-1] );
    UtilsClient::resetHisto( met06_[ism-1] );
    UtilsClient::resetHisto( met07_[ism-1] );
    UtilsClient::resetHisto( met08_[ism-1] );

    UtilsClient::resetHisto( metav01_[ism-1] );
    UtilsClient::resetHisto( metav02_[ism-1] );
    UtilsClient::resetHisto( metav03_[ism-1] );
    UtilsClient::resetHisto( metav04_[ism-1] );
    UtilsClient::resetHisto( metav05_[ism-1] );
    UtilsClient::resetHisto( metav06_[ism-1] );
    UtilsClient::resetHisto( metav07_[ism-1] );
    UtilsClient::resetHisto( metav08_[ism-1] );

    UtilsClient::resetHisto( metrms01_[ism-1] );
    UtilsClient::resetHisto( metrms02_[ism-1] );
    UtilsClient::resetHisto( metrms03_[ism-1] );
    UtilsClient::resetHisto( metrms04_[ism-1] );
    UtilsClient::resetHisto( metrms05_[ism-1] );
    UtilsClient::resetHisto( metrms06_[ism-1] );
    UtilsClient::resetHisto( metrms07_[ism-1] );
    UtilsClient::resetHisto( metrms08_[ism-1] );

    UtilsClient::resetHisto( meaopn01_[ism-1] );
    UtilsClient::resetHisto( meaopn02_[ism-1] );
    UtilsClient::resetHisto( meaopn03_[ism-1] );
    UtilsClient::resetHisto( meaopn04_[ism-1] );
    UtilsClient::resetHisto( meaopn05_[ism-1] );
    UtilsClient::resetHisto( meaopn06_[ism-1] );
    UtilsClient::resetHisto( meaopn07_[ism-1] );
    UtilsClient::resetHisto( meaopn08_[ism-1] );
    
    UtilsClient::resetHisto( mepnprms01_[ism-1] );
    UtilsClient::resetHisto( mepnprms02_[ism-1] );
    UtilsClient::resetHisto( mepnprms03_[ism-1] );
    UtilsClient::resetHisto( mepnprms04_[ism-1] );
    UtilsClient::resetHisto( mepnprms05_[ism-1] );
    UtilsClient::resetHisto( mepnprms06_[ism-1] );
    UtilsClient::resetHisto( mepnprms07_[ism-1] );
    UtilsClient::resetHisto( mepnprms08_[ism-1] );
    
  }

}

void EELaserClient::cleanup(void) {

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( cloneME_ ) {
      if ( h01_[ism-1] ) delete h01_[ism-1];
      if ( h02_[ism-1] ) delete h02_[ism-1];
      if ( h03_[ism-1] ) delete h03_[ism-1];
      if ( h04_[ism-1] ) delete h04_[ism-1];
      if ( h05_[ism-1] ) delete h05_[ism-1];
      if ( h06_[ism-1] ) delete h06_[ism-1];
      if ( h07_[ism-1] ) delete h07_[ism-1];
      if ( h08_[ism-1] ) delete h08_[ism-1];

      if ( h09_[ism-1] ) delete h09_[ism-1];
      if ( h10_[ism-1] ) delete h10_[ism-1];
      if ( h11_[ism-1] ) delete h11_[ism-1];
      if ( h12_[ism-1] ) delete h12_[ism-1];

      if ( h13_[ism-1] ) delete h13_[ism-1];
      if ( h14_[ism-1] ) delete h14_[ism-1];
      if ( h15_[ism-1] ) delete h15_[ism-1];
      if ( h16_[ism-1] ) delete h16_[ism-1];
      if ( h17_[ism-1] ) delete h17_[ism-1];
      if ( h18_[ism-1] ) delete h18_[ism-1];
      if ( h19_[ism-1] ) delete h19_[ism-1];
      if ( h20_[ism-1] ) delete h20_[ism-1];

      if ( h21_[ism-1] ) delete h21_[ism-1];
      if ( h22_[ism-1] ) delete h22_[ism-1];
      if ( h23_[ism-1] ) delete h23_[ism-1];
      if ( h24_[ism-1] ) delete h24_[ism-1];

      if ( hs01_[ism-1] ) delete hs01_[ism-1];
      if ( hs02_[ism-1] ) delete hs02_[ism-1];
      if ( hs03_[ism-1] ) delete hs03_[ism-1];
      if ( hs04_[ism-1] ) delete hs04_[ism-1];

      if ( hs05_[ism-1] ) delete hs05_[ism-1];
      if ( hs06_[ism-1] ) delete hs06_[ism-1];
      if ( hs07_[ism-1] ) delete hs07_[ism-1];
      if ( hs08_[ism-1] ) delete hs08_[ism-1];

      if ( i01_[ism-1] ) delete i01_[ism-1];
      if ( i02_[ism-1] ) delete i02_[ism-1];
      if ( i03_[ism-1] ) delete i03_[ism-1];
      if ( i04_[ism-1] ) delete i04_[ism-1];
      if ( i05_[ism-1] ) delete i05_[ism-1];
      if ( i06_[ism-1] ) delete i06_[ism-1];
      if ( i07_[ism-1] ) delete i07_[ism-1];
      if ( i08_[ism-1] ) delete i08_[ism-1];

      if ( i09_[ism-1] ) delete i09_[ism-1];
      if ( i10_[ism-1] ) delete i10_[ism-1];
      if ( i11_[ism-1] ) delete i11_[ism-1];
      if ( i12_[ism-1] ) delete i12_[ism-1];
      if ( i13_[ism-1] ) delete i13_[ism-1];
      if ( i14_[ism-1] ) delete i14_[ism-1];
      if ( i15_[ism-1] ) delete i15_[ism-1];
      if ( i16_[ism-1] ) delete i16_[ism-1];
    }

    h01_[ism-1] = 0;
    h02_[ism-1] = 0;
    h03_[ism-1] = 0;
    h04_[ism-1] = 0;
    h05_[ism-1] = 0;
    h06_[ism-1] = 0;
    h07_[ism-1] = 0;
    h08_[ism-1] = 0;

    h09_[ism-1] = 0;
    h10_[ism-1] = 0;
    h11_[ism-1] = 0;
    h12_[ism-1] = 0;

    h13_[ism-1] = 0;
    h14_[ism-1] = 0;
    h15_[ism-1] = 0;
    h16_[ism-1] = 0;
    h17_[ism-1] = 0;
    h18_[ism-1] = 0;
    h19_[ism-1] = 0;
    h20_[ism-1] = 0;

    h21_[ism-1] = 0;
    h22_[ism-1] = 0;
    h23_[ism-1] = 0;
    h24_[ism-1] = 0;

    hs01_[ism-1] = 0;
    hs02_[ism-1] = 0;
    hs03_[ism-1] = 0;
    hs04_[ism-1] = 0;

    hs05_[ism-1] = 0;
    hs06_[ism-1] = 0;
    hs07_[ism-1] = 0;
    hs08_[ism-1] = 0;

    i01_[ism-1] = 0;
    i02_[ism-1] = 0;
    i03_[ism-1] = 0;
    i04_[ism-1] = 0;
    i05_[ism-1] = 0;
    i06_[ism-1] = 0;
    i07_[ism-1] = 0;
    i08_[ism-1] = 0;

    i09_[ism-1] = 0;
    i10_[ism-1] = 0;
    i11_[ism-1] = 0;
    i12_[ism-1] = 0;
    i13_[ism-1] = 0;
    i14_[ism-1] = 0;
    i15_[ism-1] = 0;
    i16_[ism-1] = 0;

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    dbe_->setCurrentFolder( "EcalEndcap/EELaserClient" );

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
    if ( meg06_[ism-1] ) dbe_->removeElement( meg06_[ism-1]->getName() );
    meg06_[ism-1] = 0;
    if ( meg07_[ism-1] ) dbe_->removeElement( meg07_[ism-1]->getName() );
    meg07_[ism-1] = 0;
    if ( meg08_[ism-1] ) dbe_->removeElement( meg08_[ism-1]->getName() );
    meg08_[ism-1] = 0;
    if ( meg09_[ism-1] ) dbe_->removeElement( meg09_[ism-1]->getName() );
    meg09_[ism-1] = 0;
    if ( meg10_[ism-1] ) dbe_->removeElement( meg10_[ism-1]->getName() );
    meg10_[ism-1] = 0;
    if ( meg11_[ism-1] ) dbe_->removeElement( meg11_[ism-1]->getName() );
    meg11_[ism-1] = 0;
    if ( meg12_[ism-1] ) dbe_->removeElement( meg12_[ism-1]->getName() );
    meg12_[ism-1] = 0;

    if ( mea01_[ism-1] ) dbe_->removeElement( mea01_[ism-1]->getName() );
    mea01_[ism-1] = 0;
    if ( mea02_[ism-1] ) dbe_->removeElement( mea02_[ism-1]->getName() );
    mea02_[ism-1] = 0;
    if ( mea03_[ism-1] ) dbe_->removeElement( mea03_[ism-1]->getName() );
    mea03_[ism-1] = 0;
    if ( mea04_[ism-1] ) dbe_->removeElement( mea04_[ism-1]->getName() );
    mea04_[ism-1] = 0;
    if ( mea05_[ism-1] ) dbe_->removeElement( mea05_[ism-1]->getName() );
    mea05_[ism-1] = 0;
    if ( mea06_[ism-1] ) dbe_->removeElement( mea06_[ism-1]->getName() );
    mea06_[ism-1] = 0;
    if ( mea07_[ism-1] ) dbe_->removeElement( mea07_[ism-1]->getName() );
    mea07_[ism-1] = 0;
    if ( mea08_[ism-1] ) dbe_->removeElement( mea08_[ism-1]->getName() );
    mea08_[ism-1] = 0;

    if ( met01_[ism-1] ) dbe_->removeElement( met01_[ism-1]->getName() );
    met01_[ism-1] = 0;
    if ( met02_[ism-1] ) dbe_->removeElement( met02_[ism-1]->getName() );
    met02_[ism-1] = 0;
    if ( met03_[ism-1] ) dbe_->removeElement( met03_[ism-1]->getName() );
    met03_[ism-1] = 0;
    if ( met04_[ism-1] ) dbe_->removeElement( met04_[ism-1]->getName() );
    met04_[ism-1] = 0;
    if ( met05_[ism-1] ) dbe_->removeElement( met05_[ism-1]->getName() );
    met05_[ism-1] = 0;
    if ( met06_[ism-1] ) dbe_->removeElement( met06_[ism-1]->getName() );
    met06_[ism-1] = 0;
    if ( met07_[ism-1] ) dbe_->removeElement( met07_[ism-1]->getName() );
    met07_[ism-1] = 0;
    if ( met08_[ism-1] ) dbe_->removeElement( met08_[ism-1]->getName() );
    met08_[ism-1] = 0;

    if ( metav01_[ism-1] ) dbe_->removeElement( metav01_[ism-1]->getName() );
    metav01_[ism-1] = 0;
    if ( metav02_[ism-1] ) dbe_->removeElement( metav02_[ism-1]->getName() );
    metav02_[ism-1] = 0;
    if ( metav03_[ism-1] ) dbe_->removeElement( metav03_[ism-1]->getName() );
    metav03_[ism-1] = 0;
    if ( metav04_[ism-1] ) dbe_->removeElement( metav04_[ism-1]->getName() );
    metav04_[ism-1] = 0;
    if ( metav05_[ism-1] ) dbe_->removeElement( metav05_[ism-1]->getName() );
    metav05_[ism-1] = 0;
    if ( metav06_[ism-1] ) dbe_->removeElement( metav06_[ism-1]->getName() );
    metav06_[ism-1] = 0;
    if ( metav07_[ism-1] ) dbe_->removeElement( metav07_[ism-1]->getName() );
    metav07_[ism-1] = 0;
    if ( metav08_[ism-1] ) dbe_->removeElement( metav08_[ism-1]->getName() );
    metav08_[ism-1] = 0;

    if ( metrms01_[ism-1] ) dbe_->removeElement( metrms01_[ism-1]->getName() );
    metrms01_[ism-1] = 0;
    if ( metrms02_[ism-1] ) dbe_->removeElement( metrms02_[ism-1]->getName() );
    metrms02_[ism-1] = 0;
    if ( metrms03_[ism-1] ) dbe_->removeElement( metrms03_[ism-1]->getName() );
    metrms03_[ism-1] = 0;
    if ( metrms04_[ism-1] ) dbe_->removeElement( metrms04_[ism-1]->getName() );
    metrms04_[ism-1] = 0;
    if ( metrms05_[ism-1] ) dbe_->removeElement( metrms05_[ism-1]->getName() );
    metrms05_[ism-1] = 0;
    if ( metrms06_[ism-1] ) dbe_->removeElement( metrms06_[ism-1]->getName() );
    metrms06_[ism-1] = 0;
    if ( metrms07_[ism-1] ) dbe_->removeElement( metrms07_[ism-1]->getName() );
    metrms07_[ism-1] = 0;
    if ( metrms08_[ism-1] ) dbe_->removeElement( metrms08_[ism-1]->getName() );
    metrms08_[ism-1] = 0;

    if ( meaopn01_[ism-1] ) dbe_->removeElement( meaopn01_[ism-1]->getName() );
    meaopn01_[ism-1] = 0;
    if ( meaopn02_[ism-1] ) dbe_->removeElement( meaopn02_[ism-1]->getName() );
    meaopn02_[ism-1] = 0;
    if ( meaopn03_[ism-1] ) dbe_->removeElement( meaopn03_[ism-1]->getName() );
    meaopn03_[ism-1] = 0;
    if ( meaopn04_[ism-1] ) dbe_->removeElement( meaopn04_[ism-1]->getName() );
    meaopn04_[ism-1] = 0;
    if ( meaopn05_[ism-1] ) dbe_->removeElement( meaopn05_[ism-1]->getName() );
    meaopn05_[ism-1] = 0;
    if ( meaopn06_[ism-1] ) dbe_->removeElement( meaopn06_[ism-1]->getName() );
    meaopn06_[ism-1] = 0;
    if ( meaopn07_[ism-1] ) dbe_->removeElement( meaopn07_[ism-1]->getName() );
    meaopn07_[ism-1] = 0;
    if ( meaopn08_[ism-1] ) dbe_->removeElement( meaopn08_[ism-1]->getName() );
    meaopn08_[ism-1] = 0;

    if ( mepnprms01_[ism-1] ) dbe_->removeElement( mepnprms01_[ism-1]->getName() );
    mepnprms01_[ism-1] = 0;
    if ( mepnprms02_[ism-1] ) dbe_->removeElement( mepnprms02_[ism-1]->getName() );
    mepnprms02_[ism-1] = 0;
    if ( mepnprms03_[ism-1] ) dbe_->removeElement( mepnprms03_[ism-1]->getName() );
    mepnprms03_[ism-1] = 0;
    if ( mepnprms04_[ism-1] ) dbe_->removeElement( mepnprms04_[ism-1]->getName() );
    mepnprms04_[ism-1] = 0;
    if ( mepnprms05_[ism-1] ) dbe_->removeElement( mepnprms05_[ism-1]->getName() );
    mepnprms05_[ism-1] = 0;
    if ( mepnprms06_[ism-1] ) dbe_->removeElement( mepnprms06_[ism-1]->getName() );
    mepnprms06_[ism-1] = 0;
    if ( mepnprms07_[ism-1] ) dbe_->removeElement( mepnprms07_[ism-1]->getName() );
    mepnprms07_[ism-1] = 0;
    if ( mepnprms08_[ism-1] ) dbe_->removeElement( mepnprms08_[ism-1]->getName() );
    mepnprms08_[ism-1] = 0;

  }

}

bool EELaserClient::writeDb(EcalCondDBInterface* econn, RunIOV* runiov, MonRunIOV* moniov) {

  bool status = true;

  EcalLogicID ecid;

  MonLaserBlueDat apd_bl;
  map<EcalLogicID, MonLaserBlueDat> dataset1_bl;
  MonLaserGreenDat apd_gr;
  map<EcalLogicID, MonLaserGreenDat> dataset1_gr;
  MonLaserIRedDat apd_ir;
  map<EcalLogicID, MonLaserIRedDat> dataset1_ir;
  MonLaserRedDat apd_rd;
  map<EcalLogicID, MonLaserRedDat> dataset1_rd;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    cout << " SM=" << ism << endl;

    UtilsClient::printBadChannels(qth01_[ism-1]);
    UtilsClient::printBadChannels(qth05_[ism-1]);
    UtilsClient::printBadChannels(qth02_[ism-1]);
    UtilsClient::printBadChannels(qth06_[ism-1]);
    UtilsClient::printBadChannels(qth03_[ism-1]);
    UtilsClient::printBadChannels(qth07_[ism-1]);
    UtilsClient::printBadChannels(qth04_[ism-1]);
    UtilsClient::printBadChannels(qth08_[ism-1]);

//    UtilsClient::printBadChannels(qtg01_[ism-1]);
//    UtilsClient::printBadChannels(qtg02_[ism-1]);
//    UtilsClient::printBadChannels(qtg03_[ism-1]);
//    UtilsClient::printBadChannels(qtg04_[ism-1]);

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism);

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( ! Numbers::validEE(ism, jx, jy) ) continue;

        bool update01;
        bool update02;
        bool update03;
        bool update04;
        bool update05;
        bool update06;
        bool update07;
        bool update08;

        bool update09;
        bool update10;
        bool update11;
        bool update12;
        bool update13;
        bool update14;
        bool update15;
        bool update16;

        float num01, num02, num03, num04, num05, num06, num07, num08;
        float mean01, mean02, mean03, mean04, mean05, mean06, mean07, mean08;
        float rms01, rms02, rms03, rms04, rms05, rms06, rms07, rms08;

        float num09, num10, num11, num12, num13, num14, num15, num16;
        float mean09, mean10, mean11, mean12, mean13, mean14, mean15, mean16;
        float rms09, rms10, rms11, rms12, rms13, rms14, rms15, rms16;

        update01 = UtilsClient::getBinStats(h01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(h02_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(h03_[ism-1], ix, iy, num03, mean03, rms03);
        update04 = UtilsClient::getBinStats(h04_[ism-1], ix, iy, num04, mean04, rms04);
        update05 = UtilsClient::getBinStats(h05_[ism-1], ix, iy, num05, mean05, rms05);
        update06 = UtilsClient::getBinStats(h06_[ism-1], ix, iy, num06, mean06, rms06);
        update07 = UtilsClient::getBinStats(h07_[ism-1], ix, iy, num07, mean07, rms07);
        update08 = UtilsClient::getBinStats(h08_[ism-1], ix, iy, num08, mean08, rms08);

        update09 = UtilsClient::getBinStats(h13_[ism-1], ix, iy, num09, mean09, rms09);
        update10 = UtilsClient::getBinStats(h14_[ism-1], ix, iy, num10, mean10, rms10);
        update11 = UtilsClient::getBinStats(h15_[ism-1], ix, iy, num11, mean11, rms11);
        update12 = UtilsClient::getBinStats(h16_[ism-1], ix, iy, num12, mean12, rms12);
        update13 = UtilsClient::getBinStats(h17_[ism-1], ix, iy, num13, mean13, rms13);
        update14 = UtilsClient::getBinStats(h18_[ism-1], ix, iy, num14, mean14, rms14);
        update15 = UtilsClient::getBinStats(h19_[ism-1], ix, iy, num15, mean15, rms15);
        update16 = UtilsClient::getBinStats(h20_[ism-1], ix, iy, num16, mean16, rms16);

        if ( update01 || update02 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L1A (" << ix << "," << iy << ") " << num01 << " " << mean01 << " " << rms01 << endl;

            cout << endl;

          }

          apd_bl.setAPDMean(mean01);
          apd_bl.setAPDRMS(rms01);

          apd_bl.setAPDOverPNMean(mean02);
          apd_bl.setAPDOverPNRMS(rms02);

          if ( meg01_[ism-1] && int(meg01_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_bl.setTaskStatus(true);
          } else {
            apd_bl.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg01_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
              dataset1_bl[ecid] = apd_bl;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update09 || update10 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L1B (" << ix << "," << iy << ") " << num09 << " " << mean09 << " " << rms09 << endl;

            cout << endl;

          }

          apd_bl.setAPDMean(mean09);
          apd_bl.setAPDRMS(rms09);

          apd_bl.setAPDOverPNMean(mean10);
          apd_bl.setAPDOverPNRMS(rms10);

          if ( meg01_[ism-1] && int(meg01_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_bl.setTaskStatus(true);
          } else {
            apd_bl.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg01_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalBarrel), ic);
              dataset1_bl[ecid] = apd_bl;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update03 || update04 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L2A (" << ix << "," << iy << ") " << num03 << " " << mean03 << " " << rms03 << endl;

            cout << endl;

          }

          apd_ir.setAPDMean(mean03);
          apd_ir.setAPDRMS(rms03);

          apd_ir.setAPDOverPNMean(mean04);
          apd_ir.setAPDOverPNRMS(rms04);

          if ( meg02_[ism-1] && int(meg02_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_ir.setTaskStatus(true);
          } else {
            apd_ir.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg02_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
              dataset1_ir[ecid] = apd_ir;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update11 || update12 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L2B (" << ix << "," << iy << ") " << num11 << " " << mean11 << " " << rms11 << endl;

            cout << endl;

          }

          apd_ir.setAPDMean(mean11);
          apd_ir.setAPDRMS(rms11);

          apd_ir.setAPDOverPNMean(mean12); 
          apd_ir.setAPDOverPNRMS(rms12);

          if ( meg02_[ism-1] && int(meg02_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_ir.setTaskStatus(true);
          } else {
            apd_ir.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg02_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalBarrel), ic);
              dataset1_ir[ecid] = apd_ir;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update05 || update06 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L3A (" << ix << "," << iy << ") " << num05 << " " << mean05 << " " << rms05 << endl;

            cout << endl;

          }

          apd_gr.setAPDMean(mean05);
          apd_gr.setAPDRMS(rms05);

          apd_gr.setAPDOverPNMean(mean06);
          apd_gr.setAPDOverPNRMS(rms06);

          if ( meg03_[ism-1] && int(meg03_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_gr.setTaskStatus(true);
          } else {
            apd_gr.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg03_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
              dataset1_gr[ecid] = apd_gr;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update13 || update14 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L3B (" << ix << "," << iy << ") " << num13 << " " << mean13 << " " << rms13 << endl;

            cout << endl;

          }

          apd_gr.setAPDMean(mean13);
          apd_gr.setAPDRMS(rms13);

          apd_gr.setAPDOverPNMean(mean14);
          apd_gr.setAPDOverPNRMS(rms14);

          if ( meg03_[ism-1] && int(meg03_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_gr.setTaskStatus(true);
          } else {
            apd_gr.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg03_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalBarrel), ic);
              dataset1_gr[ecid] = apd_gr;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update07 || update08 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L4A (" << ix << "," << iy << ") " << num07 << " " << mean07 << " " << rms07 << endl;

            cout << endl;

          }

          apd_rd.setAPDMean(mean07);
          apd_rd.setAPDRMS(rms07);

          apd_rd.setAPDOverPNMean(mean08);
          apd_rd.setAPDOverPNRMS(rms08);

          if ( meg04_[ism-1] && int(meg04_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_rd.setTaskStatus(true);
          } else {
            apd_rd.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg04_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic == -1 ) continue;

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalEndcap), ic);
              dataset1_rd[ecid] = apd_rd;
            } catch (runtime_error &e) {
              cerr << e.what() << endl;
            }
          }

        }

        if ( update15 || update16 ) {

          if ( ix == 1 && iy == 1 ) {

            cout << "Preparing dataset for SM=" << ism << endl;

            cout << "L4B (" << ix << "," << iy << ") " << num15 << " " << mean15 << " " << rms15 << endl;

            cout << endl;

          }

          apd_rd.setAPDMean(mean15);
          apd_rd.setAPDRMS(rms15);

          apd_rd.setAPDOverPNMean(mean16);
          apd_rd.setAPDOverPNRMS(rms16);

          if ( meg04_[ism-1] && int(meg04_[ism-1]->getBinContent( ix, iy )) % 3 == 1. ) {
            apd_rd.setTaskStatus(true);
          } else {
            apd_rd.setTaskStatus(false);
          }

          status = status && UtilsClient::getBinQual(meg04_[ism-1], ix, iy);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( econn ) {
            try {
              ecid = LogicID::getEcalLogicID("EB_crystal_number", Numbers::iSM(ism, EcalBarrel), ic);
              dataset1_rd[ecid] = apd_rd;
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
      cout << "Inserting MonLaserDat ... " << flush;
      if ( dataset1_bl.size() != 0 ) econn->insertDataArraySet(&dataset1_bl, moniov);
      if ( dataset1_ir.size() != 0 ) econn->insertDataArraySet(&dataset1_ir, moniov);
      if ( dataset1_gr.size() != 0 ) econn->insertDataArraySet(&dataset1_gr, moniov);
      if ( dataset1_rd.size() != 0 ) econn->insertDataArraySet(&dataset1_rd, moniov);
      cout << "done." << endl;
    } catch (runtime_error &e) {
      cerr << e.what() << endl;
    }
  }

  cout << endl;

  MonPNBlueDat pn_bl;
  map<EcalLogicID, MonPNBlueDat> dataset2_bl;
  MonPNGreenDat pn_gr;
  map<EcalLogicID, MonPNGreenDat> dataset2_gr;
  MonPNIRedDat pn_ir;
  map<EcalLogicID, MonPNIRedDat> dataset2_ir;
  MonPNRedDat pn_rd;
  map<EcalLogicID, MonPNRedDat> dataset2_rd;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    cout << " SM=" << ism << endl;

    UtilsClient::printBadChannels(qth09_[ism-1]);
    UtilsClient::printBadChannels(qth10_[ism-1]);
    UtilsClient::printBadChannels(qth11_[ism-1]);
    UtilsClient::printBadChannels(qth12_[ism-1]);
    UtilsClient::printBadChannels(qth13_[ism-1]);
    UtilsClient::printBadChannels(qth14_[ism-1]);
    UtilsClient::printBadChannels(qth15_[ism-1]);
    UtilsClient::printBadChannels(qth16_[ism-1]);
    UtilsClient::printBadChannels(qth17_[ism-1]);
    UtilsClient::printBadChannels(qth18_[ism-1]);
    UtilsClient::printBadChannels(qth19_[ism-1]);
    UtilsClient::printBadChannels(qth20_[ism-1]);
    UtilsClient::printBadChannels(qth21_[ism-1]);
    UtilsClient::printBadChannels(qth22_[ism-1]);
    UtilsClient::printBadChannels(qth23_[ism-1]);
    UtilsClient::printBadChannels(qth24_[ism-1]);

//    UtilsClient::printBadChannels(qtg05_[ism-1]);
//    UtilsClient::printBadChannels(qtg06_[ism-1]);
//    UtilsClient::printBadChannels(qtg07_[ism-1]);
//    UtilsClient::printBadChannels(qtg08_[ism-1]);
//    UtilsClient::printBadChannels(qtg09_[ism-1]);
//    UtilsClient::printBadChannels(qtg10_[ism-1]);
//    UtilsClient::printBadChannels(qtg11_[ism-1]);
//    UtilsClient::printBadChannels(qtg12_[ism-1]);

    for ( int i = 1; i <= 10; i++ ) {

      bool update01;
      bool update02;
      bool update03;
      bool update04;
      bool update05;
      bool update06;
      bool update07;
      bool update08;
      bool update09;
      bool update10;
      bool update11;
      bool update12;
      bool update13;
      bool update14;
      bool update15;
      bool update16;

      float num01, num02, num03, num04, num05, num06, num07, num08;
      float num09, num10, num11, num12, num13, num14, num15, num16;
      float mean01, mean02, mean03, mean04, mean05, mean06, mean07, mean08;
      float mean09, mean10, mean11, mean12, mean13, mean14, mean15, mean16;
      float rms01, rms02, rms03, rms04, rms05, rms06, rms07, rms08;
      float rms09, rms10, rms11, rms12, rms13, rms14, rms15, rms16;

      update01 = UtilsClient::getBinStats(i01_[ism-1], 1, i, num01, mean01, rms01);
      update02 = UtilsClient::getBinStats(i02_[ism-1], 1, i, num02, mean02, rms02);
      update03 = UtilsClient::getBinStats(i03_[ism-1], 1, i, num03, mean03, rms03);
      update04 = UtilsClient::getBinStats(i04_[ism-1], 1, i, num04, mean04, rms04);
      update05 = UtilsClient::getBinStats(i05_[ism-1], 1, i, num05, mean05, rms05);
      update06 = UtilsClient::getBinStats(i06_[ism-1], 1, i, num06, mean06, rms06);
      update07 = UtilsClient::getBinStats(i07_[ism-1], 1, i, num07, mean07, rms07);
      update08 = UtilsClient::getBinStats(i08_[ism-1], 1, i, num08, mean08, rms08);
      update09 = UtilsClient::getBinStats(i09_[ism-1], 1, i, num09, mean09, rms09);
      update10 = UtilsClient::getBinStats(i10_[ism-1], 1, i, num10, mean10, rms10);
      update11 = UtilsClient::getBinStats(i11_[ism-1], 1, i, num11, mean11, rms11);
      update12 = UtilsClient::getBinStats(i12_[ism-1], 1, i, num12, mean12, rms12);
      update13 = UtilsClient::getBinStats(i13_[ism-1], 1, i, num13, mean13, rms13);
      update14 = UtilsClient::getBinStats(i14_[ism-1], 1, i, num14, mean14, rms14);
      update15 = UtilsClient::getBinStats(i15_[ism-1], 1, i, num15, mean15, rms15);
      update16 = UtilsClient::getBinStats(i16_[ism-1], 1, i, num16, mean16, rms16);

      if ( update01 || update05 || update09 || update13 ) {

        if ( i == 1 ) {

          cout << "Preparing dataset for SM=" << ism << endl;

          cout << "PNs (" << i << ") L1 G01 " << num01  << " " << mean01 << " " << rms01  << endl;
          cout << "PNs (" << i << ") L1 G16 " << num09  << " " << mean09 << " " << rms09  << endl;

          cout << endl;

        }

        pn_bl.setADCMeanG1(mean01);
        pn_bl.setADCRMSG1(rms01);

        pn_bl.setPedMeanG1(mean05);
        pn_bl.setPedRMSG1(rms05);

        pn_bl.setADCMeanG16(mean09);
        pn_bl.setADCRMSG16(rms09);

        pn_bl.setPedMeanG16(mean13);
        pn_bl.setPedRMSG16(rms13);

        if ( meg05_[ism-1] && int(meg05_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ||
             meg09_[ism-1] && int(meg09_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ) {
          pn_bl.setTaskStatus(true);
        } else {
          pn_bl.setTaskStatus(false);
        }

        status = status && ( UtilsClient::getBinQual(meg05_[ism-1], i, 1) ||
                             UtilsClient::getBinQual(meg09_[ism-1], i, 1) );

        if ( econn ) {
          try {
            ecid = LogicID::getEcalLogicID("EB_LM_PN", Numbers::iSM(ism, EcalEndcap), i-1);
            dataset2_bl[ecid] = pn_bl;
          } catch (runtime_error &e) {
            cerr << e.what() << endl;
          }
        }

      }

      if ( update02 || update06 || update10 || update14 ) {

        if ( i == 1 ) {

          cout << "Preparing dataset for SM=" << ism << endl;

          cout << "PNs (" << i << ") L2 G01 " << num02  << " " << mean02 << " " << rms02  << endl;
          cout << "PNs (" << i << ") L2 G16 " << num10  << " " << mean10 << " " << rms10  << endl;

          cout << endl;

        }

        pn_ir.setADCMeanG1(mean02);
        pn_ir.setADCRMSG1(rms02);

        pn_ir.setPedMeanG1(mean06);
        pn_ir.setPedRMSG1(rms06);

        pn_ir.setADCMeanG16(mean10);
        pn_ir.setADCRMSG16(rms10);

        pn_ir.setPedMeanG16(mean14);
        pn_ir.setPedRMSG16(rms14);

        if ( meg06_[ism-1] && int(meg06_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ||
             meg10_[ism-1] && int(meg10_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ) {
          pn_ir.setTaskStatus(true);
        } else {
          pn_ir.setTaskStatus(false);
        }

        status = status && ( UtilsClient::getBinQual(meg06_[ism-1], i, 1) ||
                             UtilsClient::getBinQual(meg10_[ism-1], i, 1) );

        if ( econn ) {
          try {
            ecid = LogicID::getEcalLogicID("EB_LM_PN", Numbers::iSM(ism, EcalEndcap), i-1);
            dataset2_ir[ecid] = pn_ir;
          } catch (runtime_error &e) {
            cerr << e.what() << endl;
          }
        }

      }

      if ( update03 || update07 || update11 || update15 ) {

        if ( i == 1 ) {

          cout << "Preparing dataset for SM=" << ism << endl;

          cout << "PNs (" << i << ") L3 G01 " << num03  << " " << mean03 << " " << rms03  << endl;
          cout << "PNs (" << i << ") L3 G16 " << num11  << " " << mean11 << " " << rms11  << endl;

          cout << endl;

        }

        pn_gr.setADCMeanG1(mean03);
        pn_gr.setADCRMSG1(rms03);

        pn_gr.setPedMeanG1(mean07);
        pn_gr.setPedRMSG1(rms07);

        pn_gr.setADCMeanG16(mean11);
        pn_gr.setADCRMSG16(rms11);

        pn_gr.setPedMeanG16(mean15);
        pn_gr.setPedRMSG16(rms15);

        if ( meg07_[ism-1] && int(meg07_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ||
             meg11_[ism-1] && int(meg11_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ) {
          pn_gr.setTaskStatus(true);
        } else {
          pn_gr.setTaskStatus(false);
        }

        status = status && ( UtilsClient::getBinQual(meg07_[ism-1], i, 1) ||
                             UtilsClient::getBinQual(meg11_[ism-1], i, 1) );

        if ( econn ) {
          try {
            ecid = LogicID::getEcalLogicID("EB_LM_PN", Numbers::iSM(ism, EcalEndcap), i-1);
            dataset2_gr[ecid] = pn_gr;
          } catch (runtime_error &e) {
            cerr << e.what() << endl;
          }
        }

      }

      if ( update04 || update08 || update12 || update16 ) {

        if ( i == 1 ) {

          cout << "Preparing dataset for SM=" << ism << endl;

          cout << "PNs (" << i << ") L4 G01 " << num04  << " " << mean04 << " " << rms04  << endl;
          cout << "PNs (" << i << ") L4 G16 " << num12  << " " << mean12 << " " << rms12  << endl;

          cout << endl;

        }

        pn_rd.setADCMeanG1(mean04);
        pn_rd.setADCRMSG1(rms04);

        pn_rd.setPedMeanG1(mean08);
        pn_rd.setPedRMSG1(mean08);

        pn_rd.setADCMeanG16(mean12);
        pn_rd.setADCRMSG16(rms12);

        pn_rd.setPedMeanG16(mean16);
        pn_rd.setPedRMSG16(rms16);

        if ( meg08_[ism-1] && int(meg08_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ||
             meg12_[ism-1] && int(meg12_[ism-1]->getBinContent( i, 1 )) % 3 == 1. ) {
          pn_rd.setTaskStatus(true);
        } else {
          pn_rd.setTaskStatus(false);
        }

        status = status && ( UtilsClient::getBinQual(meg08_[ism-1], i, 1) ||
                             UtilsClient::getBinQual(meg12_[ism-1], i, 1) );

        if ( econn ) {
          try {
            ecid = LogicID::getEcalLogicID("EB_LM_PN", Numbers::iSM(ism, EcalEndcap), i-1);
            dataset2_rd[ecid] = pn_rd;
          } catch (runtime_error &e) {
            cerr << e.what() << endl;
          }
        }

      }

    }

  }

  if ( econn ) {
    try {
      cout << "Inserting MonPnDat ... " << flush;
      if ( dataset2_bl.size() != 0 ) econn->insertDataArraySet(&dataset2_bl, moniov);
      if ( dataset2_ir.size() != 0 ) econn->insertDataArraySet(&dataset2_ir, moniov);
      if ( dataset2_gr.size() != 0 ) econn->insertDataArraySet(&dataset2_gr, moniov);
      if ( dataset2_rd.size() != 0 ) econn->insertDataArraySet(&dataset2_rd, moniov);
      cout << "done." << endl;
    } catch (runtime_error &e) {
      cerr << e.what() << endl;
    }
  }

  return status;

}

void EELaserClient::subscribe(void){

  if ( verbose_ ) cout << "EELaserClient: subscribe" << endl;

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->subscribe(histo, ism);

  }

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    if ( enableMonitorDaemon_ ) {
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A", Numbers::sEE(ism).c_str());
      if ( qth01_[ism-1] ) dbe_->useQTest(histo, qth01_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A", Numbers::sEE(ism).c_str());
      if ( qth02_[ism-1] ) dbe_->useQTest(histo, qth02_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A", Numbers::sEE(ism).c_str());
      if ( qth03_[ism-1] ) dbe_->useQTest(histo, qth03_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A", Numbers::sEE(ism).c_str());
      if ( qth04_[ism-1] ) dbe_->useQTest(histo, qth04_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B", Numbers::sEE(ism).c_str());
      if ( qth05_[ism-1] ) dbe_->useQTest(histo, qth05_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B", Numbers::sEE(ism).c_str());
      if ( qth06_[ism-1] ) dbe_->useQTest(histo, qth06_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B", Numbers::sEE(ism).c_str());
      if ( qth07_[ism-1] ) dbe_->useQTest(histo, qth07_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B", Numbers::sEE(ism).c_str());
      if ( qth08_[ism-1] ) dbe_->useQTest(histo, qth08_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1", Numbers::sEE(ism).c_str());
      if ( qth09_[ism-1] ) dbe_->useQTest(histo, qth09_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2", Numbers::sEE(ism).c_str());
      if ( qth10_[ism-1] ) dbe_->useQTest(histo, qth10_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3", Numbers::sEE(ism).c_str());
      if ( qth11_[ism-1] ) dbe_->useQTest(histo, qth11_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4", Numbers::sEE(ism).c_str());
      if ( qth12_[ism-1] ) dbe_->useQTest(histo, qth12_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1", Numbers::sEE(ism).c_str());
      if ( qth13_[ism-1] ) dbe_->useQTest(histo, qth13_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2", Numbers::sEE(ism).c_str());
      if ( qth14_[ism-1] ) dbe_->useQTest(histo, qth14_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3", Numbers::sEE(ism).c_str());
      if ( qth15_[ism-1] ) dbe_->useQTest(histo, qth15_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4", Numbers::sEE(ism).c_str());
      if ( qth16_[ism-1] ) dbe_->useQTest(histo, qth16_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1", Numbers::sEE(ism).c_str());
      if ( qth17_[ism-1] ) dbe_->useQTest(histo, qth17_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2", Numbers::sEE(ism).c_str());
      if ( qth18_[ism-1] ) dbe_->useQTest(histo, qth18_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3", Numbers::sEE(ism).c_str());
      if ( qth19_[ism-1] ) dbe_->useQTest(histo, qth19_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4", Numbers::sEE(ism).c_str());
      if ( qth20_[ism-1] ) dbe_->useQTest(histo, qth20_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1", Numbers::sEE(ism).c_str());
      if ( qth21_[ism-1] ) dbe_->useQTest(histo, qth21_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2", Numbers::sEE(ism).c_str());
      if ( qth22_[ism-1] ) dbe_->useQTest(histo, qth22_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3", Numbers::sEE(ism).c_str());
      if ( qth23_[ism-1] ) dbe_->useQTest(histo, qth23_[ism-1]->getName());
      sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4", Numbers::sEE(ism).c_str());
      if ( qth24_[ism-1] ) dbe_->useQTest(histo, qth24_[ism-1]->getName());
    } else {
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A", Numbers::sEE(ism).c_str());
      if ( qth01_[ism-1] ) dbe_->useQTest(histo, qth01_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A", Numbers::sEE(ism).c_str());
      if ( qth02_[ism-1] ) dbe_->useQTest(histo, qth02_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A", Numbers::sEE(ism).c_str());
      if ( qth03_[ism-1] ) dbe_->useQTest(histo, qth03_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A", Numbers::sEE(ism).c_str());
      if ( qth04_[ism-1] ) dbe_->useQTest(histo, qth04_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B", Numbers::sEE(ism).c_str());
      if ( qth05_[ism-1] ) dbe_->useQTest(histo, qth05_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B", Numbers::sEE(ism).c_str());
      if ( qth06_[ism-1] ) dbe_->useQTest(histo, qth06_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B", Numbers::sEE(ism).c_str());
      if ( qth07_[ism-1] ) dbe_->useQTest(histo, qth07_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B", Numbers::sEE(ism).c_str());
      if ( qth08_[ism-1] ) dbe_->useQTest(histo, qth08_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1", Numbers::sEE(ism).c_str());
      if ( qth09_[ism-1] ) dbe_->useQTest(histo, qth09_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2", Numbers::sEE(ism).c_str());
      if ( qth10_[ism-1] ) dbe_->useQTest(histo, qth10_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3", Numbers::sEE(ism).c_str());
      if ( qth11_[ism-1] ) dbe_->useQTest(histo, qth11_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4", Numbers::sEE(ism).c_str());
      if ( qth12_[ism-1] ) dbe_->useQTest(histo, qth12_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1", Numbers::sEE(ism).c_str());
      if ( qth13_[ism-1] ) dbe_->useQTest(histo, qth13_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2", Numbers::sEE(ism).c_str());
      if ( qth14_[ism-1] ) dbe_->useQTest(histo, qth14_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3", Numbers::sEE(ism).c_str());
      if ( qth15_[ism-1] ) dbe_->useQTest(histo, qth15_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4", Numbers::sEE(ism).c_str());
      if ( qth16_[ism-1] ) dbe_->useQTest(histo, qth16_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1", Numbers::sEE(ism).c_str());
      if ( qth17_[ism-1] ) dbe_->useQTest(histo, qth17_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2", Numbers::sEE(ism).c_str());
      if ( qth18_[ism-1] ) dbe_->useQTest(histo, qth18_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3", Numbers::sEE(ism).c_str());
      if ( qth19_[ism-1] ) dbe_->useQTest(histo, qth19_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4", Numbers::sEE(ism).c_str());
      if ( qth20_[ism-1] ) dbe_->useQTest(histo, qth20_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1", Numbers::sEE(ism).c_str());
      if ( qth21_[ism-1] ) dbe_->useQTest(histo, qth21_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2", Numbers::sEE(ism).c_str());
      if ( qth22_[ism-1] ) dbe_->useQTest(histo, qth22_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3", Numbers::sEE(ism).c_str());
      if ( qth23_[ism-1] ) dbe_->useQTest(histo, qth23_[ism-1]->getName());
      sprintf(histo, "EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4", Numbers::sEE(ism).c_str());
      if ( qth24_[ism-1] ) dbe_->useQTest(histo, qth24_[ism-1]->getName());
    }

    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L1 %s", Numbers::sEE(ism).c_str());
    if ( qtg01_[ism-1] ) dbe_->useQTest(histo, qtg01_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L2 %s", Numbers::sEE(ism).c_str());
    if ( qtg02_[ism-1] ) dbe_->useQTest(histo, qtg02_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L3 %s", Numbers::sEE(ism).c_str());
    if ( qtg03_[ism-1] ) dbe_->useQTest(histo, qtg03_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L4 %s", Numbers::sEE(ism).c_str());
    if ( qtg04_[ism-1] ) dbe_->useQTest(histo, qtg04_[ism-1]->getName());

    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L1 PNs %s G01", Numbers::sEE(ism).c_str());
    if ( qtg05_[ism-1] ) dbe_->useQTest(histo, qtg05_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L2 PNs %s G01", Numbers::sEE(ism).c_str());
    if ( qtg06_[ism-1] ) dbe_->useQTest(histo, qtg06_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L3 PNs %s G01", Numbers::sEE(ism).c_str());
    if ( qtg07_[ism-1] ) dbe_->useQTest(histo, qtg07_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L4 PNs %s G01", Numbers::sEE(ism).c_str());
    if ( qtg08_[ism-1] ) dbe_->useQTest(histo, qtg08_[ism-1]->getName());

    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L1 PNs %s G16", Numbers::sEE(ism).c_str());
    if ( qtg09_[ism-1] ) dbe_->useQTest(histo, qtg09_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L2 PNs %s G16", Numbers::sEE(ism).c_str());
    if ( qtg10_[ism-1] ) dbe_->useQTest(histo, qtg10_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L3 PNs %s G16", Numbers::sEE(ism).c_str());
    if ( qtg11_[ism-1] ) dbe_->useQTest(histo, qtg11_[ism-1]->getName());
    sprintf(histo, "EcalEndcap/EELaserTask/EELT laser quality L4 PNs %s G16", Numbers::sEE(ism).c_str());
    if ( qtg12_[ism-1] ) dbe_->useQTest(histo, qtg12_[ism-1]->getName());

  }

}

void EELaserClient::subscribeNew(void){

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4A", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4B", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->subscribeNew(histo, ism);

  }

}

void EELaserClient::unsubscribe(void){

  if ( verbose_ ) cout << "EELaserClient: unsubscribe" << endl;

  Char_t histo[200];

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    unsigned int ism = superModules_[i];

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT timing %s L1B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT timing %s L2B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT timing %s L3B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT timing %s L4B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4A", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/EELT shape %s L1B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/EELT shape %s L2B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/EELT shape %s L3B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/EELT shape %s L4B", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);
    sprintf(histo, "*/EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4", Numbers::sEE(ism).c_str());
    mui_->unsubscribe(histo, ism);

  }

}

void EELaserClient::softReset(void){

}

void EELaserClient::analyze(void){

  ievt_++;
  jevt_++;
  if ( ievt_ % 10 == 0 ) {
    if ( verbose_ ) cout << "EELaserClient: ievt/jevt = " << ievt_ << "/" << jevt_ << endl;
  }

  uint64_t bits01 = 0;
  bits01 |= EcalErrorDictionary::getMask("LASER_MEAN_WARNING");
  bits01 |= EcalErrorDictionary::getMask("LASER_RMS_WARNING");
  bits01 |= EcalErrorDictionary::getMask("LASER_MEAN_OVER_PN_WARNING");
  bits01 |= EcalErrorDictionary::getMask("LASER_RMS_OVER_PN_WARNING");

  uint64_t bits02 = 0;
  bits02 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_MEAN_WARNING");
  bits02 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_RMS_WARNING");
  bits02 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_MEAN_ERROR");
  bits02 |= EcalErrorDictionary::getMask("PEDESTAL_LOW_GAIN_RMS_ERROR");

  uint64_t bits03 = 0;
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_MEAN_WARNING");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_RMS_WARNING");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_MEAN_ERROR");
  bits03 |= EcalErrorDictionary::getMask("PEDESTAL_MIDDLE_GAIN_RMS_ERROR");

  uint64_t bits04 = 0;
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_MEAN_WARNING");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_RMS_WARNING");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_MEAN_ERROR");
  bits04 |= EcalErrorDictionary::getMask("PEDESTAL_HIGH_GAIN_RMS_ERROR");

  map<EcalLogicID, RunCrystalErrorsDat> mask1;
  map<EcalLogicID, RunPNErrorsDat> mask2;

  EcalErrorMask::fetchDataSet(&mask1);
  EcalErrorMask::fetchDataSet(&mask2);

  Char_t histo[200];

  MonitorElement* me;

  for ( unsigned int i=0; i<superModules_.size(); i++ ) {

    int ism = superModules_[i];

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h02_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h04_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h04_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h05_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h05_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h06_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h06_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h07_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h07_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h08_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h08_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT timing %s L1A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h09_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h09_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT timing %s L2A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h10_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h10_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT timing %s L3A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h11_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h11_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT timing %s L4A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h12_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h12_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT amplitude %s L1B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h13_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h13_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT amplitude over PN %s L1B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h14_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h14_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT amplitude %s L2B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h15_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h15_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT amplitude over PN %s L2B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h16_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h16_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT amplitude %s L3B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h17_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h17_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT amplitude over PN %s L3B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h18_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h18_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT amplitude %s L4B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h19_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h19_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT amplitude over PN %s L4B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h20_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h20_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT timing %s L1B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h21_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h21_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT timing %s L2B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h22_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h22_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT timing %s L3B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h23_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h23_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT timing %s L4B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    h24_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, h24_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT shape %s L1A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT shape %s L2A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs02_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT shape %s L3A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT shape %s L4A").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs04_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs04_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/EELT shape %s L1B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs05_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs05_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/EELT shape %s L2B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs06_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs06_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/EELT shape %s L3B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs07_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs07_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/EELT shape %s L4B").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    hs08_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, hs08_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs amplitude %s G01 L1").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i01_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i01_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs amplitude %s G01 L2").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i02_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i02_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs amplitude %s G01 L3").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i03_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i03_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs amplitude %s G01 L4").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i04_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i04_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/PN/Gain01/EEPDT PNs pedestal %s G01 L1").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i05_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i05_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/PN/Gain01/EEPDT PNs pedestal %s G01 L2").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i06_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i06_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/PN/Gain01/EEPDT PNs pedestal %s G01 L3").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i07_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i07_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/PN/Gain01/EEPDT PNs pedestal %s G01 L4").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i08_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i08_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs amplitude %s G16 L1").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i09_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i09_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs amplitude %s G16 L2").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i10_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i10_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs amplitude %s G16 L3").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i11_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i11_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs amplitude %s G16 L4").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i12_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i12_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser1/PN/Gain16/EEPDT PNs pedestal %s G16 L1").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i13_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i13_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser2/PN/Gain16/EEPDT PNs pedestal %s G16 L2").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i14_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i14_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser3/PN/Gain16/EEPDT PNs pedestal %s G16 L3").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i15_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i15_[ism-1] );

    sprintf(histo, (prefixME_+"EcalEndcap/EELaserTask/Laser4/PN/Gain16/EEPDT PNs pedestal %s G16 L4").c_str(), Numbers::sEE(ism).c_str());
    me = dbe_->get(histo);
    i16_[ism-1] = UtilsClient::getHisto<TProfile2D*>( me, cloneME_, i16_[ism-1] );

    UtilsClient::resetHisto( meg01_[ism-1] );
    UtilsClient::resetHisto( meg02_[ism-1] );
    UtilsClient::resetHisto( meg03_[ism-1] );
    UtilsClient::resetHisto( meg04_[ism-1] );

    UtilsClient::resetHisto( meg05_[ism-1] );
    UtilsClient::resetHisto( meg06_[ism-1] );
    UtilsClient::resetHisto( meg07_[ism-1] );
    UtilsClient::resetHisto( meg08_[ism-1] );
    UtilsClient::resetHisto( meg09_[ism-1] );
    UtilsClient::resetHisto( meg10_[ism-1] );
    UtilsClient::resetHisto( meg11_[ism-1] );
    UtilsClient::resetHisto( meg12_[ism-1] );

    UtilsClient::resetHisto( mea01_[ism-1] );
    UtilsClient::resetHisto( mea02_[ism-1] );
    UtilsClient::resetHisto( mea03_[ism-1] );
    UtilsClient::resetHisto( mea04_[ism-1] );
    UtilsClient::resetHisto( mea05_[ism-1] );
    UtilsClient::resetHisto( mea06_[ism-1] );
    UtilsClient::resetHisto( mea07_[ism-1] );
    UtilsClient::resetHisto( mea08_[ism-1] );

    UtilsClient::resetHisto( met01_[ism-1] );
    UtilsClient::resetHisto( met02_[ism-1] );
    UtilsClient::resetHisto( met03_[ism-1] );
    UtilsClient::resetHisto( met04_[ism-1] );
    UtilsClient::resetHisto( met05_[ism-1] );
    UtilsClient::resetHisto( met06_[ism-1] );
    UtilsClient::resetHisto( met07_[ism-1] );
    UtilsClient::resetHisto( met08_[ism-1] );

    UtilsClient::resetHisto( metav01_[ism-1] );
    UtilsClient::resetHisto( metav02_[ism-1] );
    UtilsClient::resetHisto( metav03_[ism-1] );
    UtilsClient::resetHisto( metav04_[ism-1] );
    UtilsClient::resetHisto( metav05_[ism-1] );
    UtilsClient::resetHisto( metav06_[ism-1] );
    UtilsClient::resetHisto( metav07_[ism-1] );
    UtilsClient::resetHisto( metav08_[ism-1] );

    UtilsClient::resetHisto( metrms01_[ism-1] );
    UtilsClient::resetHisto( metrms02_[ism-1] );
    UtilsClient::resetHisto( metrms03_[ism-1] );
    UtilsClient::resetHisto( metrms04_[ism-1] );
    UtilsClient::resetHisto( metrms05_[ism-1] );
    UtilsClient::resetHisto( metrms06_[ism-1] );
    UtilsClient::resetHisto( metrms07_[ism-1] );
    UtilsClient::resetHisto( metrms08_[ism-1] );

    UtilsClient::resetHisto( meaopn01_[ism-1] );
    UtilsClient::resetHisto( meaopn02_[ism-1] );
    UtilsClient::resetHisto( meaopn03_[ism-1] );
    UtilsClient::resetHisto( meaopn04_[ism-1] );
    UtilsClient::resetHisto( meaopn05_[ism-1] );
    UtilsClient::resetHisto( meaopn06_[ism-1] );
    UtilsClient::resetHisto( meaopn07_[ism-1] );
    UtilsClient::resetHisto( meaopn08_[ism-1] );

    UtilsClient::resetHisto( mepnprms01_[ism-1] );
    UtilsClient::resetHisto( mepnprms02_[ism-1] );
    UtilsClient::resetHisto( mepnprms03_[ism-1] );
    UtilsClient::resetHisto( mepnprms04_[ism-1] );
    UtilsClient::resetHisto( mepnprms05_[ism-1] );
    UtilsClient::resetHisto( mepnprms06_[ism-1] );
    UtilsClient::resetHisto( mepnprms07_[ism-1] );
    UtilsClient::resetHisto( mepnprms08_[ism-1] );
    
    float meanAmplL1A, meanAmplL2A, meanAmplL3A, meanAmplL4A;
    float meanAmplL1B, meanAmplL2B, meanAmplL3B, meanAmplL4B;

    int nCryL1A, nCryL2A, nCryL3A, nCryL4A;
    int nCryL1B, nCryL2B, nCryL3B, nCryL4B;

    meanAmplL1A = meanAmplL2A = meanAmplL3A = meanAmplL4A = 0.;
    meanAmplL1B = meanAmplL2B = meanAmplL3B = meanAmplL4B = 0.;

    nCryL1A = nCryL2A = nCryL3A = nCryL4A = 0;
    nCryL1B = nCryL2B = nCryL3B = nCryL4B = 0;

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        bool update01;
        bool update02;
        bool update03;
        bool update04;
        bool update05;
        bool update06;
        bool update07;
        bool update08;

        float num01, num02, num03, num04, num05, num06, num07, num08;
        float mean01, mean02, mean03, mean04, mean05, mean06, mean07, mean08;
        float rms01, rms02, rms03, rms04, rms05, rms06, rms07, rms08;

        update01 = UtilsClient::getBinStats(h01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(h03_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(h05_[ism-1], ix, iy, num03, mean03, rms03);
        update04 = UtilsClient::getBinStats(h07_[ism-1], ix, iy, num04, mean04, rms04);
        update05 = UtilsClient::getBinStats(h13_[ism-1], ix, iy, num05, mean05, rms05);
        update06 = UtilsClient::getBinStats(h15_[ism-1], ix, iy, num06, mean06, rms06);
        update07 = UtilsClient::getBinStats(h17_[ism-1], ix, iy, num07, mean07, rms07);
        update08 = UtilsClient::getBinStats(h19_[ism-1], ix, iy, num08, mean08, rms08);

        if ( update01 ) {
          meanAmplL1A += mean01;
          nCryL1A++;
        }

        if ( update02 ) {
          meanAmplL2A += mean02;
          nCryL2A++;
        }

        if ( update03 ) {
          meanAmplL3A += mean03;
          nCryL3A++;
        }

        if ( update04 ) {
          meanAmplL4A += mean04;
          nCryL4A++;
        }

        if ( update05 ) {
          meanAmplL1B += mean05;
          nCryL1B++;
        }

        if ( update06 ) {
          meanAmplL2B += mean06;
          nCryL2B++;
        }

        if ( update07 ) {
          meanAmplL3B += mean07;
          nCryL3B++;
        }

        if ( update08 ) {
          meanAmplL4B += mean08;
          nCryL4B++;
        }

      }
    }

    if ( nCryL1A > 0 ) meanAmplL1A /= float (nCryL1A);
    if ( nCryL2A > 0 ) meanAmplL2A /= float (nCryL2A);
    if ( nCryL3A > 0 ) meanAmplL3A /= float (nCryL3A);
    if ( nCryL4A > 0 ) meanAmplL4A /= float (nCryL4A);
    if ( nCryL1B > 0 ) meanAmplL1B /= float (nCryL1B);
    if ( nCryL2B > 0 ) meanAmplL2B /= float (nCryL2B);
    if ( nCryL3B > 0 ) meanAmplL3B /= float (nCryL3B);
    if ( nCryL4B > 0 ) meanAmplL4B /= float (nCryL4B);

    for ( int ix = 1; ix <= 50; ix++ ) {
      for ( int iy = 1; iy <= 50; iy++ ) {

        if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, -1.);
        if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, -1.);
        if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, -1.);
        if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent( ix, iy, -1.);

        int jx = ix + Numbers::ix0EE(ism);
        int jy = iy + Numbers::iy0EE(ism); 

        if ( ism >= 1 && ism <= 9 ) jx = 101 - jx;

        if ( Numbers::validEE(ism, jx, jy) ) {
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, 2.);
          if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, 2.);
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, 2.);
          if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent( ix, iy, 2.);
        }

        bool update01;
        bool update02;
        bool update03;
        bool update04;
        bool update05;
        bool update06;
        bool update07;
        bool update08;
        bool update09;
        bool update10;
        bool update11;
        bool update12;

        bool update13;
        bool update14;
        bool update15;
        bool update16;
        bool update17;
        bool update18;
        bool update19;
        bool update20;
        bool update21;
        bool update22;
        bool update23;
        bool update24;

        float num01, num02, num03, num04, num05, num06, num07, num08;
        float num09, num10, num11, num12;
        float mean01, mean02, mean03, mean04, mean05, mean06, mean07, mean08;
        float mean09, mean10, mean11, mean12;
        float rms01, rms02, rms03, rms04, rms05, rms06, rms07, rms08;
        float rms09, rms10, rms11, rms12;

        float num13, num14, num15, num16, num17, num18, num19, num20;
        float num21, num22, num23, num24;
        float mean13, mean14, mean15, mean16, mean17, mean18, mean19, mean20;
        float mean21, mean22, mean23, mean24;
        float rms13, rms14, rms15, rms16, rms17, rms18, rms19, rms20;
        float rms21, rms22, rms23, rms24;

        update01 = UtilsClient::getBinStats(h01_[ism-1], ix, iy, num01, mean01, rms01);
        update02 = UtilsClient::getBinStats(h02_[ism-1], ix, iy, num02, mean02, rms02);
        update03 = UtilsClient::getBinStats(h03_[ism-1], ix, iy, num03, mean03, rms03);
        update04 = UtilsClient::getBinStats(h04_[ism-1], ix, iy, num04, mean04, rms04);
        update05 = UtilsClient::getBinStats(h05_[ism-1], ix, iy, num05, mean05, rms05);
        update06 = UtilsClient::getBinStats(h06_[ism-1], ix, iy, num06, mean06, rms06);
        update07 = UtilsClient::getBinStats(h07_[ism-1], ix, iy, num07, mean07, rms07);
        update08 = UtilsClient::getBinStats(h08_[ism-1], ix, iy, num08, mean08, rms08);
        update09 = UtilsClient::getBinStats(h09_[ism-1], ix, iy, num09, mean09, rms09);
        update10 = UtilsClient::getBinStats(h10_[ism-1], ix, iy, num10, mean10, rms10);
        update11 = UtilsClient::getBinStats(h11_[ism-1], ix, iy, num11, mean11, rms11);
        update12 = UtilsClient::getBinStats(h12_[ism-1], ix, iy, num12, mean12, rms12);
	
         // other SM half
	
        update13 = UtilsClient::getBinStats(h13_[ism-1], ix, iy, num13, mean13, rms13);
        update14 = UtilsClient::getBinStats(h14_[ism-1], ix, iy, num14, mean14, rms14);
        update15 = UtilsClient::getBinStats(h15_[ism-1], ix, iy, num15, mean15, rms15);
        update16 = UtilsClient::getBinStats(h16_[ism-1], ix, iy, num16, mean16, rms16);
        update17 = UtilsClient::getBinStats(h17_[ism-1], ix, iy, num17, mean17, rms17);
        update18 = UtilsClient::getBinStats(h18_[ism-1], ix, iy, num18, mean18, rms18);
        update19 = UtilsClient::getBinStats(h19_[ism-1], ix, iy, num19, mean19, rms19);
        update20 = UtilsClient::getBinStats(h20_[ism-1], ix, iy, num20, mean20, rms20);
        update21 = UtilsClient::getBinStats(h21_[ism-1], ix, iy, num21, mean21, rms21);
        update22 = UtilsClient::getBinStats(h22_[ism-1], ix, iy, num22, mean22, rms22);
        update23 = UtilsClient::getBinStats(h23_[ism-1], ix, iy, num23, mean23, rms23);
        update24 = UtilsClient::getBinStats(h24_[ism-1], ix, iy, num24, mean24, rms24);

        if ( update01 ) {

          float val;

          val = 1.;
          if ( fabs(mean01 - meanAmplL1A) > fabs(percentVariation_ * meanAmplL1A) )
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

        if ( update13 ) {

          float val;

          val = 1.;
          if ( fabs(mean13 - meanAmplL1B) > fabs(percentVariation_ * meanAmplL1B) )
           val = 0.;
          if ( meg01_[ism-1] ) meg01_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea05_[ism-1] ) {
              if ( mean13 > 0. ) {
                mea05_[ism-1]->setBinContent( ic, mean13 );
                mea05_[ism-1]->setBinError( ic, rms13 );
            } else {
                mea05_[ism-1]->setEntries( 1.+mea05_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update03 ) {

          float val;

          val = 1.;
          if ( fabs(mean03 - meanAmplL2A) > fabs(percentVariation_ * meanAmplL2A) )
            val = 0.;
          if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, val);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea02_[ism-1] ) {
              if ( mean03 > 0. ) {
                mea02_[ism-1]->setBinContent( ic, mean03 );
                mea02_[ism-1]->setBinError( ic, rms03 );
              } else {
                mea02_[ism-1]->setEntries( 1.+mea02_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update15 ) {

          float val;

          val = 1.;
          if ( fabs(mean15 - meanAmplL2B) > fabs(percentVariation_ * meanAmplL2B) )
            val = 0.;
          if ( meg02_[ism-1] ) meg02_[ism-1]->setBinContent( ix, iy, val);

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea06_[ism-1] ) {
              if ( mean15 > 0. ) {
                mea06_[ism-1]->setBinContent( ic, mean15 );
                mea06_[ism-1]->setBinError( ic, rms15 );
              } else {
                mea06_[ism-1]->setEntries( 1.+mea06_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update05 ) {

          float val;

          val = 1.;
          if ( fabs(mean05 - meanAmplL3A) > fabs(percentVariation_ * meanAmplL3A) )
            val = 0.;
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea03_[ism-1] ) {
              if ( mean05 > 0. ) {
                mea03_[ism-1]->setBinContent( ic, mean05 );
                mea03_[ism-1]->setBinError( ic, rms05 );
              } else {
                mea03_[ism-1]->setEntries( 1.+mea03_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update17 ) {

          float val;

          val = 1.;
          if ( fabs(mean17 - meanAmplL3B) > fabs(percentVariation_ * meanAmplL3B) )
            val = 0.;
          if ( meg03_[ism-1] ) meg03_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea07_[ism-1] ) {
              if ( mean17 > 0. ) {
                mea07_[ism-1]->setBinContent( ic, mean17 );
                mea07_[ism-1]->setBinError( ic, rms17 );
              } else {
                mea07_[ism-1]->setEntries( 1.+mea07_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update07 ) {

          float val;

          val = 1.;
          if ( fabs(mean07 - meanAmplL4A) > fabs(percentVariation_ * meanAmplL4A) )
            val = 0.;
          if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea04_[ism-1] ) {
              if ( mean07 > 0. ) {
                mea04_[ism-1]->setBinContent( ic, mean07 );
                mea04_[ism-1]->setBinError( ic, rms07 );
              } else {
                mea04_[ism-1]->setEntries( 1.+mea04_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update19 ) {

          float val;

          val = 1.;
          if ( fabs(mean19 - meanAmplL4B) > fabs(percentVariation_ * meanAmplL4B) )
            val = 0.;
          if ( meg04_[ism-1] ) meg04_[ism-1]->setBinContent( ix, iy, val );

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( mea08_[ism-1] ) {
              if ( mean19 > 0. ) {
                mea08_[ism-1]->setBinContent( ic, mean19 );
                mea08_[ism-1]->setBinError( ic, rms19 );
              } else {
                mea08_[ism-1]->setEntries( 1.+mea08_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update02 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn01_[ism-1] ) {
              if ( mean02 > 0. ) {
                meaopn01_[ism-1]->setBinContent( ic, mean02 );
                meaopn01_[ism-1]->setBinError( ic, rms02 );
              } else {
                meaopn01_[ism-1]->setEntries( 1.+meaopn01_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update14 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn05_[ism-1] ) {
              if ( mean14 > 0. ) {
                meaopn05_[ism-1]->setBinContent( ic, mean14 );
                meaopn05_[ism-1]->setBinError( ic, rms14 );
              } else {
                meaopn05_[ism-1]->setEntries( 1.+meaopn05_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update04 ) {


          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn02_[ism-1] ) {
              if ( mean04 > 0. ) {
                meaopn02_[ism-1]->setBinContent( ic, mean04 );
                meaopn02_[ism-1]->setBinError( ic, rms04 );
              } else {
                meaopn02_[ism-1]->setEntries( 1.+meaopn02_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update16 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn06_[ism-1] ) {
              if ( mean16 > 0. ) {
                meaopn06_[ism-1]->setBinContent( ic, mean16 );
                meaopn06_[ism-1]->setBinError( ic, rms16 );
              } else {
                meaopn06_[ism-1]->setEntries( 1.+meaopn06_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update06 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn03_[ism-1] ) {
              if ( mean06 > 0. ) {
                meaopn03_[ism-1]->setBinContent( ic, mean06 );
                meaopn03_[ism-1]->setBinError( ic, rms06 );
              } else {
                meaopn03_[ism-1]->setEntries( 1.+meaopn03_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update18 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn07_[ism-1] ) {
              if ( mean18 > 0. ) {
                meaopn07_[ism-1]->setBinContent( ic, mean18 );
                meaopn07_[ism-1]->setBinError( ic, rms18 );
              } else {
                meaopn07_[ism-1]->setEntries( 1.+meaopn07_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update08 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn04_[ism-1] ) {
              if ( mean08 > 0. ) {
                meaopn04_[ism-1]->setBinContent( ic, mean08 );
                meaopn04_[ism-1]->setBinError( ic, rms08 );
              } else {
                meaopn04_[ism-1]->setEntries( 1.+meaopn04_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update20 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( meaopn08_[ism-1] ) {
              if ( mean20 > 0. ) {
                meaopn08_[ism-1]->setBinContent( ic, mean20 );
                meaopn08_[ism-1]->setBinError( ic, rms20 );
              } else {
                meaopn08_[ism-1]->setEntries( 1.+meaopn08_[ism-1]->getEntries() );
              }
            }
          }

        }

        if ( update09 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met01_[ism-1] ) {
              if ( mean09 > 0. ) {
                met01_[ism-1]->setBinContent( ic, mean09 );
                met01_[ism-1]->setBinError( ic, rms09 );
              } else {
                met01_[ism-1]->setEntries(1.+met01_[ism-1]->getEntries());
              }
            }

            if ( metav01_[ism-1] ) metav01_[ism-1] ->Fill(mean09);
            if ( metrms01_[ism-1] ) metrms01_[ism-1]->Fill(rms09);
          }

        }

        if ( update21 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met05_[ism-1] ) {
              if ( mean21 > 0. ) {
                met05_[ism-1]->setBinContent( ic, mean21 );
                met05_[ism-1]->setBinError( ic, rms21 );
              } else {
                met05_[ism-1]->setEntries(1.+met05_[ism-1]->getEntries());
              }
            }

            if ( metav05_[ism-1] ) metav05_[ism-1] ->Fill(mean21);
            if ( metrms05_[ism-1] ) metrms05_[ism-1]->Fill(rms21);
          }

        }

        if ( update10 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met02_[ism-1] ) {
              if ( mean10 > 0. ) {
                met02_[ism-1]->setBinContent( ic, mean10 );
                met02_[ism-1]->setBinError( ic, rms10 );
              } else {
                met02_[ism-1]->setEntries(1.+met02_[ism-1]->getEntries());
              }
            }

            if ( metav02_[ism-1] ) metav02_[ism-1] ->Fill(mean10);
            if ( metrms02_[ism-1] ) metrms02_[ism-1]->Fill(rms10);
          }

        }

        if ( update22 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met06_[ism-1] ) {
              if ( mean22 > 0. ) {
                met06_[ism-1]->setBinContent( ic, mean22 );
                met06_[ism-1]->setBinError( ic, rms22 );
              } else {
                met06_[ism-1]->setEntries(1.+met06_[ism-1]->getEntries());
              }
            }

            if ( metav06_[ism-1] ) metav06_[ism-1] ->Fill(mean22);
            if ( metrms06_[ism-1] ) metrms06_[ism-1]->Fill(rms22);
          }

        }

        if ( update11 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met03_[ism-1] ) {
              if ( mean11 > 0. ) {
                met03_[ism-1]->setBinContent( ic, mean11 );
                met03_[ism-1]->setBinError( ic, rms11 );
              } else {
                met03_[ism-1]->setEntries(1.+met03_[ism-1]->getEntries());
              }
            }

            if ( metav03_[ism-1] ) metav03_[ism-1] ->Fill(mean11);
            if ( metrms03_[ism-1] ) metrms03_[ism-1]->Fill(rms11);
          }

        }

        if ( update23 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met07_[ism-1] ) {
              if ( mean23 > 0. ) {
                met07_[ism-1]->setBinContent( ic, mean23 );
                met07_[ism-1]->setBinError( ic, rms23 );
              } else {
                met07_[ism-1]->setEntries(1.+met07_[ism-1]->getEntries());
              }
            }

            if ( metav07_[ism-1] ) metav07_[ism-1] ->Fill(mean23);
            if ( metrms07_[ism-1] ) metrms07_[ism-1]->Fill(rms23);
          }

        }

        if ( update12 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met04_[ism-1] ) {
              if ( mean12 > 0. ) {
                met04_[ism-1]->setBinContent( ic, mean12 );
                met04_[ism-1]->setBinError( ic, rms12 );
              } else {
                met04_[ism-1]->setEntries(1.+met04_[ism-1]->getEntries());
              }
            }

            if ( metav04_[ism-1] ) metav04_[ism-1] ->Fill(mean12);
            if ( metrms04_[ism-1] ) metrms04_[ism-1]->Fill(rms12);
          }

        }

        if ( update24 ) {

          int ic = Numbers::icEE(ism, ix, iy);

          if ( ic != -1 ) {
            if ( met08_[ism-1] ) {
              if ( mean24 > 0. ) {
                met08_[ism-1]->setBinContent( ic, mean24 );
                met08_[ism-1]->setBinError( ic, rms24 );
              } else {
                met08_[ism-1]->setEntries(1.+met08_[ism-1]->getEntries());
              }
            }

            if ( metav08_[ism-1] ) metav08_[ism-1] ->Fill(mean24);
            if ( metrms08_[ism-1] ) metrms08_[ism-1]->Fill(rms24);
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

            int ic = Numbers::icEE(ism, ix, iy);

            if ( ic == -1 ) continue;

            EcalLogicID ecid = m->first;

            if ( ecid.getID1() == Numbers::iSM(ism, EcalEndcap) && ecid.getID2() == ic ) {
              if ( (m->second).getErrorBits() & bits01 ) {
                if ( meg01_[ism-1] ) {
                  float val = int(meg01_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg01_[ism-1]->setBinContent( ix, iy, val+3 );
                }
                if ( meg02_[ism-1] ) {
                  float val = int(meg02_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg02_[ism-1]->setBinContent( ix, iy, val+3 );
                }
                if ( meg03_[ism-1] ) {
                  float val = int(meg03_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg03_[ism-1]->setBinContent( ix, iy, val+3 );
                }
                if ( meg04_[ism-1] ) {
                  float val = int(meg04_[ism-1]->getBinContent(ix, iy)) % 3;
                  meg04_[ism-1]->setBinContent( ix, iy, val+3 );
                }
              }
            }

          }
        }

      }
    }

    for ( int i = 1; i <= 10; i++ ) {

      if ( meg05_[ism-1] ) meg05_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg06_[ism-1] ) meg06_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg07_[ism-1] ) meg07_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg08_[ism-1] ) meg08_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg09_[ism-1] ) meg09_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg10_[ism-1] ) meg10_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg11_[ism-1] ) meg11_[ism-1]->setBinContent( i, 1, 2. );
      if ( meg12_[ism-1] ) meg12_[ism-1]->setBinContent( i, 1, 2. );

      bool update01;
      bool update02;
      bool update03;
      bool update04;
      bool update05;
      bool update06;
      bool update07;
      bool update08;
      bool update09;
      bool update10;
      bool update11;
      bool update12;
      bool update13;
      bool update14;
      bool update15;
      bool update16;

      float num01, num02, num03, num04, num05, num06, num07, num08;
      float num09, num10, num11, num12, num13, num14, num15, num16;
      float mean01, mean02, mean03, mean04, mean05, mean06, mean07, mean08;
      float mean09, mean10, mean11, mean12, mean13, mean14, mean15, mean16;
      float rms01, rms02, rms03, rms04, rms05, rms06, rms07, rms08;
      float rms09, rms10, rms11, rms12, rms13, rms14, rms15, rms16;

      update01 = UtilsClient::getBinStats(i01_[ism-1], 1, i, num01, mean01, rms01);
      update02 = UtilsClient::getBinStats(i02_[ism-1], 1, i, num02, mean02, rms02);
      update03 = UtilsClient::getBinStats(i03_[ism-1], 1, i, num03, mean03, rms03);
      update04 = UtilsClient::getBinStats(i04_[ism-1], 1, i, num04, mean04, rms04);
      update05 = UtilsClient::getBinStats(i05_[ism-1], 1, i, num05, mean05, rms05);
      update06 = UtilsClient::getBinStats(i06_[ism-1], 1, i, num06, mean06, rms06);
      update07 = UtilsClient::getBinStats(i07_[ism-1], 1, i, num07, mean07, rms07);
      update08 = UtilsClient::getBinStats(i08_[ism-1], 1, i, num08, mean08, rms08);
      update09 = UtilsClient::getBinStats(i09_[ism-1], 1, i, num09, mean09, rms09);
      update10 = UtilsClient::getBinStats(i10_[ism-1], 1, i, num10, mean10, rms10);
      update11 = UtilsClient::getBinStats(i11_[ism-1], 1, i, num11, mean11, rms11);
      update12 = UtilsClient::getBinStats(i12_[ism-1], 1, i, num12, mean12, rms12);
      update13 = UtilsClient::getBinStats(i13_[ism-1], 1, i, num13, mean13, rms13);
      update14 = UtilsClient::getBinStats(i14_[ism-1], 1, i, num14, mean14, rms14);
      update15 = UtilsClient::getBinStats(i15_[ism-1], 1, i, num15, mean15, rms15);
      update16 = UtilsClient::getBinStats(i16_[ism-1], 1, i, num16, mean16, rms16);

      if ( update01 && update05 ) {

        float val;

        val = 1.;
        if ( mean01 < amplitudeThresholdPnG01_ )
          val = 0.;
        if ( mean05 <  pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0] ||
	     pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0] < mean05)
          val = 0.;
        if ( rms05 > pedPnRMSThreshold_[0] )
          val = 0.;

        if ( meg05_[ism-1] ) meg05_[ism-1]->setBinContent(i, 1, val);
        if ( mepnprms01_[ism-1] ) mepnprms01_[ism-1]->Fill(rms05);

      }

      if ( update02 && update06 ) {

        float val;

        val = 1.;
        if ( mean02 < amplitudeThresholdPnG01_ )
          val = 0.;
        if ( mean06 <  pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0] ||
	     pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0] < mean06)
          val = 0.;
        if ( rms06 > pedPnRMSThreshold_[0] )
          val = 0.;

        if ( meg06_[ism-1] )           meg06_[ism-1]->setBinContent(i, 1, val);
        if ( mepnprms02_[ism-1] ) mepnprms02_[ism-1]->Fill(rms06);

      }

      if ( update03 && update07 ) {

        float val;

        val = 1.;
        if ( mean03 < amplitudeThresholdPnG01_ )
          val = 0.;
      if ( mean07 <  pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0] ||
	   pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0] < mean07)
          val = 0.;
      if ( rms07 > pedPnRMSThreshold_[0] )
          val = 0.;

      if ( meg07_[ism-1] )           meg07_[ism-1]->setBinContent(i, 1, val);
      if ( mepnprms03_[ism-1] ) mepnprms03_[ism-1]->Fill(rms07);

      }

      if ( update04 && update08 ) {

        float val;

        val = 1.;
        if ( mean04 < amplitudeThresholdPnG01_ )
          val = 0.;
       if ( mean08 <  pedPnExpectedMean_[0] - pedPnDiscrepancyMean_[0] ||
	    pedPnExpectedMean_[0] + pedPnDiscrepancyMean_[0] < mean08)
           val = 0.;
       if ( rms08 > pedPnRMSThreshold_[0] )
           val = 0.;
       
       if ( meg08_[ism-1] )           meg08_[ism-1]->setBinContent(i, 1, val);
       if ( mepnprms04_[ism-1] ) mepnprms04_[ism-1]->Fill(rms08);
       
      }
      
      if ( update09 && update13 ) {

        float val;

        val = 1.;
        if ( mean09 < amplitudeThresholdPnG16_ )
          val = 0.;
        if ( mean13 <  pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1] ||
	     pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1] < mean13)
          val = 0.;
        if ( rms13 > pedPnRMSThreshold_[1] )
          val = 0.;
	
        if ( meg09_[ism-1] )           meg09_[ism-1]->setBinContent(i, 1, val);
        if ( mepnprms05_[ism-1] ) mepnprms05_[ism-1]->Fill(rms13);

      }

      if ( update10 && update14 ) {

        float val;

        val = 1.;
        if ( mean10 < amplitudeThresholdPnG16_ )
          val = 0.;
	//        if ( mean14 < pedestalThresholdPn_ )
       if ( mean14 <  pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1] ||
	    pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1] < mean14)
           val = 0.;
       if ( rms14 > pedPnRMSThreshold_[1] )
          val = 0.;
       
       if ( meg10_[ism-1] )           meg10_[ism-1]->setBinContent(i, 1, val);
       if ( mepnprms06_[ism-1] ) mepnprms06_[ism-1]->Fill(rms14);
       
      }

      if ( update11 && update15 ) {

        float val;

        val = 1.;
        if ( mean11 < amplitudeThresholdPnG16_ )
          val = 0.;
	//        if ( mean15 < pedestalThresholdPn_ )
        if ( mean15 <  pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1] ||
	     pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1] < mean15)
          val = 0.;
        if ( rms15 > pedPnRMSThreshold_[1] )
          val = 0.;

        if ( meg11_[ism-1] ) meg11_[ism-1]->setBinContent(i, 1, val);
        if ( mepnprms07_[ism-1] ) mepnprms07_[ism-1]->Fill(rms15);

      }

      if ( update12 && update16 ) {

        float val;

        val = 1.;
        if ( mean12 < amplitudeThresholdPnG16_ )
          val = 0.;
	//        if ( mean16 < pedestalThresholdPn_ )
        if ( mean16 <  pedPnExpectedMean_[1] - pedPnDiscrepancyMean_[1] ||
  	              pedPnExpectedMean_[1] + pedPnDiscrepancyMean_[1] < mean16)
          val = 0.;
        if ( rms16 > pedPnRMSThreshold_[1] )
          val = 0.;

        if ( meg12_[ism-1] ) meg12_[ism-1]->setBinContent(i, 1, val);
        if ( mepnprms08_[ism-1] ) mepnprms08_[ism-1]->Fill(rms16);

      }

      // masking

      if ( mask2.size() != 0 ) {
        map<EcalLogicID, RunPNErrorsDat>::const_iterator m;
        for (m = mask2.begin(); m != mask2.end(); m++) {

          EcalLogicID ecid = m->first;

          if ( ecid.getID1() == Numbers::iSM(ism, EcalEndcap) && ecid.getID2() == i-1 ) {
            if ( (m->second).getErrorBits() & (bits01|bits02) ) {
              if ( meg05_[ism-1] ) {
                float val = int(meg05_[ism-1]->getBinContent(i, 1)) % 3;
                meg05_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits02) ) {
              if ( meg06_[ism-1] ) {
                float val = int(meg06_[ism-1]->getBinContent(i, 1)) % 3;
                meg06_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits02) ) {
              if ( meg07_[ism-1] ) {
                float val = int(meg07_[ism-1]->getBinContent(i, 1)) % 3;
                meg07_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits02) ) {
              if ( meg08_[ism-1] ) {
                float val = int(meg08_[ism-1]->getBinContent(i, 1)) % 3;
                meg08_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits04) ) {
              if ( meg09_[ism-1] ) {
                float val = int(meg09_[ism-1]->getBinContent(i, 1)) % 3;
                meg09_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits04) ) {
              if ( meg10_[ism-1] ) {
                float val = int(meg10_[ism-1]->getBinContent(i, 1)) % 3;
                meg10_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits04) ) {
              if ( meg11_[ism-1] ) {
                float val = int(meg11_[ism-1]->getBinContent(i, 1)) % 3;
                meg11_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
            if ( (m->second).getErrorBits() & (bits01|bits04) ) {
              if ( meg12_[ism-1] ) {
                float val = int(meg12_[ism-1]->getBinContent(i, 1)) % 3;
                meg12_[ism-1]->setBinContent( i, 1, val+3 );
              }
            }
          }

        }
      }

    }

    vector<dqm::me_util::Channel> badChannels01;
    vector<dqm::me_util::Channel> badChannels02;
    vector<dqm::me_util::Channel> badChannels03;
    vector<dqm::me_util::Channel> badChannels04;
    vector<dqm::me_util::Channel> badChannels05;
    vector<dqm::me_util::Channel> badChannels06;
    vector<dqm::me_util::Channel> badChannels07;
    vector<dqm::me_util::Channel> badChannels08;

    if ( qth01_[ism-1] ) badChannels01 = qth01_[ism-1]->getBadChannels();
    if ( qth02_[ism-1] ) badChannels02 = qth02_[ism-1]->getBadChannels();
    if ( qth03_[ism-1] ) badChannels03 = qth03_[ism-1]->getBadChannels();
    if ( qth04_[ism-1] ) badChannels04 = qth04_[ism-1]->getBadChannels();
    if ( qth05_[ism-1] ) badChannels05 = qth05_[ism-1]->getBadChannels();
    if ( qth06_[ism-1] ) badChannels06 = qth06_[ism-1]->getBadChannels();
    if ( qth07_[ism-1] ) badChannels07 = qth07_[ism-1]->getBadChannels();
    if ( qth08_[ism-1] ) badChannels08 = qth08_[ism-1]->getBadChannels();

    vector<dqm::me_util::Channel> badChannels09;
    vector<dqm::me_util::Channel> badChannels10;
    vector<dqm::me_util::Channel> badChannels11;
    vector<dqm::me_util::Channel> badChannels12;
    vector<dqm::me_util::Channel> badChannels13;
    vector<dqm::me_util::Channel> badChannels14;
    vector<dqm::me_util::Channel> badChannels15;
    vector<dqm::me_util::Channel> badChannels16;
    vector<dqm::me_util::Channel> badChannels17;
    vector<dqm::me_util::Channel> badChannels18;
    vector<dqm::me_util::Channel> badChannels19;
    vector<dqm::me_util::Channel> badChannels20;
    vector<dqm::me_util::Channel> badChannels21;
    vector<dqm::me_util::Channel> badChannels22;
    vector<dqm::me_util::Channel> badChannels23;
    vector<dqm::me_util::Channel> badChannels24;

    if ( qth09_[ism-1] ) badChannels09 = qth09_[ism-1]->getBadChannels();
    if ( qth10_[ism-1] ) badChannels10 = qth10_[ism-1]->getBadChannels();
    if ( qth11_[ism-1] ) badChannels11 = qth11_[ism-1]->getBadChannels();
    if ( qth12_[ism-1] ) badChannels12 = qth12_[ism-1]->getBadChannels();
    if ( qth13_[ism-1] ) badChannels13 = qth13_[ism-1]->getBadChannels();
    if ( qth14_[ism-1] ) badChannels14 = qth14_[ism-1]->getBadChannels();
    if ( qth15_[ism-1] ) badChannels15 = qth15_[ism-1]->getBadChannels();
    if ( qth16_[ism-1] ) badChannels16 = qth16_[ism-1]->getBadChannels();
    if ( qth17_[ism-1] ) badChannels09 = qth17_[ism-1]->getBadChannels();
    if ( qth18_[ism-1] ) badChannels10 = qth18_[ism-1]->getBadChannels();
    if ( qth19_[ism-1] ) badChannels11 = qth19_[ism-1]->getBadChannels();
    if ( qth20_[ism-1] ) badChannels12 = qth20_[ism-1]->getBadChannels();
    if ( qth21_[ism-1] ) badChannels13 = qth21_[ism-1]->getBadChannels();
    if ( qth22_[ism-1] ) badChannels14 = qth22_[ism-1]->getBadChannels();
    if ( qth23_[ism-1] ) badChannels15 = qth23_[ism-1]->getBadChannels();
    if ( qth24_[ism-1] ) badChannels16 = qth24_[ism-1]->getBadChannels();

  }

}

void EELaserClient::htmlOutput(int run, string htmlDir, string htmlName){

  cout << "Preparing EELaserClient html output ..." << endl;

  ofstream htmlFile;

  htmlFile.open((htmlDir + htmlName).c_str());

  // html page header
  htmlFile << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">  " << endl;
  htmlFile << "<html>  " << endl;
  htmlFile << "<head>  " << endl;
  htmlFile << "  <meta content=\"text/html; charset=ISO-8859-1\"  " << endl;
  htmlFile << " http-equiv=\"content-type\">  " << endl;
  htmlFile << "  <title>Monitor:LaserTask output</title> " << endl;
  htmlFile << "</head>  " << endl;
  htmlFile << "<style type=\"text/css\"> td { font-weight: bold } </style>" << endl;
  htmlFile << "<body>  " << endl;
  //htmlFile << "<br>  " << endl;
  htmlFile << "<a name=""top""></a>" << endl;
  htmlFile << "<h2>Run:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;" << endl;
  htmlFile << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">" << run << "</span></h2>" << endl;
  htmlFile << "<h2>Monitoring task:&nbsp;&nbsp;&nbsp;&nbsp; <span " << endl;
  htmlFile << " style=\"color: rgb(0, 0, 153);\">LASER</span></h2> " << endl;
  htmlFile << "<hr>" << endl;
  htmlFile << "<table border=1><tr><td bgcolor=red>channel has problems in this task</td>" << endl;
  htmlFile << "<td bgcolor=lime>channel has NO problems</td>" << endl;
  htmlFile << "<td bgcolor=yellow>channel is missing</td></table>" << endl;
  htmlFile << "<hr>" << endl;
  //   htmlFile << "<table border=1><tr><td>L1 = blu laser</td>" << endl;
  //   htmlFile << "<td>L2 = green laser</td>" << endl;
  //   htmlFile << "<td>L3 = red laser</td>" << endl;
  //   htmlFile << "<td>L4 = infrared laser</td></table>" << endl;
  htmlFile << "<table style=\"width: 600px;\" border=\"0\">" << endl;
  htmlFile << "<tbody>" << endl;
  htmlFile << "<tr>" << endl;
  htmlFile << "<td style=\"text-align: center;\">" << endl;
  htmlFile << "<div style=\"text-align: center;\"> </div>" << endl;
  htmlFile << "<table style=\"width: 482px; height: 35px;\" border=\"1\">" << endl;
  htmlFile << "<tbody>" << endl;
  htmlFile << "<tr>" << endl;
  htmlFile << "<td style=\"text-align: center;\">L1 = blue laser </td>" << endl;
  htmlFile << "<td style=\"vertical-align: top; text-align: center;\">L2 =green laser </td>" << endl;
  htmlFile << "<td style=\"vertical-align: top; text-align: center;\">L3 =red laser </td>" << endl;
  htmlFile << "<td style=\"vertical-align: top; text-align: center;\">L4 =infrared laser </td>" << endl;
  htmlFile << "</tr>" << endl;
  htmlFile << "</tbody>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "</td>" << endl;
  htmlFile << "<td align=\"center\">" << endl;
  htmlFile << "<div style=\"text-align: center;\"> </div>" << endl;
  htmlFile << "<table style=\"width: 255px; height: 35px;\" border=\"1\">" << endl;
  htmlFile << "<tbody>" << endl;
  htmlFile << "<tr>" << endl;
  htmlFile << "<td style=\"text-align: center;\">A=first half</td>" << endl;
  htmlFile << "<td style=\"vertical-align: top; text-align: center;\">B=second half<br>" << endl;
  htmlFile << "</td>" << endl;
  htmlFile << "</tr>" << endl;
  htmlFile << "</tbody>" << endl;
  htmlFile << "</table>" << endl;
  htmlFile << "</td>" << endl;
  htmlFile << "</tr>" << endl;
  htmlFile << "</tbody>" << endl;
  htmlFile << "</table>" << endl;
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

  const double histMax = 1.e15;

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
  
  string imgNameQual[8], imgNameAmp[8], imgNameTim[8], imgNameTimav[8], imgNameTimrms[8], imgNameShape[8], imgNameAmpoPN[8], imgNameMEPnQualG01[8], imgNameMEPnG01[8], imgNameMEPnPedG01[8], imgNameMEPnRmsPedG01[8], imgNameMEPnQualG16[8], imgNameMEPnG16[8], imgNameMEPnPedG16[8], imgNameMEPnRmsPedG16[8], imgName, meName;

  TCanvas* cQual   = new TCanvas("cQual", "Temp", 2*csize, 2*csize);
  TCanvas* cQualPN = new TCanvas("cQualPN", "Temp", 2*csize, csize);
  TCanvas* cAmp    = new TCanvas("cAmp", "Temp", csize, csize);
  TCanvas* cTim    = new TCanvas("cTim", "Temp", csize, csize);
  TCanvas* cTimav  = new TCanvas("cTimav", "Temp", csize, csize);
  TCanvas* cTimrms = new TCanvas("cTimrms", "Temp", csize, csize);
  TCanvas* cShape  = new TCanvas("cShape", "Temp", csize, csize);
  TCanvas* cAmpoPN = new TCanvas("cAmpoPN", "Temp", csize, csize);
  TCanvas* cPed    = new TCanvas("cPed", "Temp", csize, csize);

  TH2F* obj2f;
  TH1F* obj1f;
  TH1D* obj1d;

  // Loop on barrel supermodules

  for ( unsigned int i=0; i<superModules_.size(); i ++ ) {

    int ism = superModules_[i];

    // Loop on wavelength times 2 'sides'

    for ( int iCanvas = 1 ; iCanvas <= 4 * 2 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;
      if ( iCanvas == 4+2 || iCanvas == 4+3 ) continue;

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
        case 4:
          obj2f = UtilsClient::getHisto<TH2F*>( meg04_[ism-1] );
          break;
        case 5:
        case 6:
        case 7:
        case 8:
          obj2f = 0;
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
        obj2f->GetYaxis()->SetLabelSize(0.02);
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
        case 4:
          obj1f = UtilsClient::getHisto<TH1F*>( mea04_[ism-1] );
          break;
        case 5:
          obj1f = UtilsClient::getHisto<TH1F*>( mea05_[ism-1] );
          break;
        case 6:
          obj1f = UtilsClient::getHisto<TH1F*>( mea06_[ism-1] );
          break;
        case 7:
          obj1f = UtilsClient::getHisto<TH1F*>( mea07_[ism-1] );
          break;
        case 8:
          obj1f = UtilsClient::getHisto<TH1F*>( mea08_[ism-1] );
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
        cAmp->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Timing distributions

      imgNameTim[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( met01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( met02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( met03_[ism-1] );
          break;
        case 4:
          obj1f = UtilsClient::getHisto<TH1F*>( met04_[ism-1] );
          break;
        case 5:
          obj1f = UtilsClient::getHisto<TH1F*>( met05_[ism-1] );
          break;
        case 6:
          obj1f = UtilsClient::getHisto<TH1F*>( met06_[ism-1] );
          break;
        case 7:
          obj1f = UtilsClient::getHisto<TH1F*>( met07_[ism-1] );
          break;
        case 8:
          obj1f = UtilsClient::getHisto<TH1F*>( met08_[ism-1] );
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
        imgNameTim[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameTim[iCanvas-1];

        cTim->cd();
        gStyle->SetOptStat("euo");
        obj1f->SetStats(kTRUE);
        obj1f->SetMinimum(0.0);
        obj1f->SetMaximum(10.0);
        obj1f->Draw();
        cTim->Update();
        cTim->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Timing mean distributions

      imgNameTimav[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( metav01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( metav02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( metav03_[ism-1] );
          break;
        case 4:
          obj1f = UtilsClient::getHisto<TH1F*>( metav04_[ism-1] );
          break;
        case 5:
          obj1f = UtilsClient::getHisto<TH1F*>( metav05_[ism-1] );
          break;
        case 6:
          obj1f = UtilsClient::getHisto<TH1F*>( metav06_[ism-1] );
          break;
        case 7:
          obj1f = UtilsClient::getHisto<TH1F*>( metav07_[ism-1] );
          break;
        case 8:
          obj1f = UtilsClient::getHisto<TH1F*>( metav08_[ism-1] );
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
        imgNameTimav[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameTimav[iCanvas-1];

        cTimav->cd();
        gStyle->SetOptStat("euomr");
        obj1f->SetStats(kTRUE);
        if ( obj1f->GetMaximum(histMax) > 0. ) {
          gPad->SetLogy(1);
        } else {
          gPad->SetLogy(0);
        }
        obj1f->Draw();
        cTimav->Update();
        cTimav->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Timing rms distributions

      imgNameTimrms[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms03_[ism-1] );
          break;
        case 4:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms04_[ism-1] );
          break;
        case 5:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms05_[ism-1] );
          break;
        case 6:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms06_[ism-1] );
          break;
        case 7:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms07_[ism-1] );
          break;
        case 8:
          obj1f = UtilsClient::getHisto<TH1F*>( metrms08_[ism-1] );
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
        imgNameTimrms[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameTimrms[iCanvas-1];

        cTimrms->cd();
        gStyle->SetOptStat("euomr");
        obj1f->SetStats(kTRUE);
        if ( obj1f->GetMaximum(histMax) > 0. ) {
          gPad->SetLogy(1);
        } else {
          gPad->SetLogy(0);
        }
        obj1f->Draw();
        cTimrms->Update();
        cTimrms->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Shape distributions

      imgNameShape[iCanvas-1] = "";

      obj1d = 0;
      switch ( iCanvas ) {
        case 1:
          if ( hs01_[ism-1] ) obj1d = hs01_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 2:
          if ( hs02_[ism-1] ) obj1d = hs02_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 3:
          if ( hs03_[ism-1] ) obj1d = hs03_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 4:
          if ( hs04_[ism-1] ) obj1d = hs04_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 5:
          if ( hs05_[ism-1] ) obj1d = hs05_[ism-1]->ProjectionY("_py", 1681, 1681, "e");
          break;
        case 6:
          if ( hs06_[ism-1] ) obj1d = hs06_[ism-1]->ProjectionY("_py", 1681, 1681, "e");
          break;
        case 7:
          if ( hs07_[ism-1] ) obj1d = hs07_[ism-1]->ProjectionY("_py", 1681, 1681, "e");
          break;
        case 8:
          if ( hs08_[ism-1] ) obj1d = hs08_[ism-1]->ProjectionY("_py", 1681, 1681, "e");
          break;
        default:
          break;
      }

      if ( obj1d ) {
        meName = obj1d->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1, "_");
          }
        }
        imgNameShape[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameShape[iCanvas-1];

        cShape->cd();
        gStyle->SetOptStat("euo");
        obj1d->SetStats(kTRUE);
//        if ( obj1d->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1d->Draw();
        cShape->Update();
        cShape->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

        delete obj1d;

      }

      // Amplitude over PN distributions

      imgNameAmpoPN[iCanvas-1] = "";

      obj1f = 0;
      switch ( iCanvas ) {
        case 1:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn01_[ism-1] );
          break;
        case 2:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn02_[ism-1] );
          break;
        case 3:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn03_[ism-1] );
          break;
        case 4:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn04_[ism-1] );
          break;
        case 5:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn05_[ism-1] );
          break;
        case 6:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn06_[ism-1] );
          break;
        case 7:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn07_[ism-1] );
          break;
        case 8:
          obj1f = UtilsClient::getHisto<TH1F*>( meaopn08_[ism-1] );
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
        imgNameAmpoPN[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameAmpoPN[iCanvas-1];

        cAmpoPN->cd();
        gStyle->SetOptStat("euo");
        obj1f->SetStats(kTRUE);
//        if ( obj1f->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1f->SetMinimum(0.0);
        obj1f->SetMaximum(20.0);
        obj1f->Draw();
        cAmpoPN->Update();
        cAmpoPN->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

      }

      // Monitoring elements plots

      imgNameMEPnQualG01[iCanvas-1] = "";

      obj2f = 0;
      switch ( iCanvas ) {
      case 1:
        obj2f = UtilsClient::getHisto<TH2F*>( meg05_[ism-1] );
        break;
      case 2:
        obj2f = UtilsClient::getHisto<TH2F*>( meg06_[ism-1] );
        break;
      case 3:
        obj2f = UtilsClient::getHisto<TH2F*>( meg07_[ism-1] );
        break;
      case 4:
        obj2f = UtilsClient::getHisto<TH2F*>( meg08_[ism-1] );
        break;
      case 5:
      case 6:
      case 7:
      case 8:
        obj2f = 0;
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
        imgNameMEPnQualG01[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnQualG01[iCanvas-1];

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

      imgNameMEPnQualG16[iCanvas-1] = "";

      obj2f = 0;
      switch ( iCanvas ) {
      case 1:
        obj2f = UtilsClient::getHisto<TH2F*>( meg09_[ism-1] );
        break;
      case 2:
        obj2f = UtilsClient::getHisto<TH2F*>( meg10_[ism-1] );
        break;
      case 3:
        obj2f = UtilsClient::getHisto<TH2F*>( meg11_[ism-1] );
        break;
      case 4:
        obj2f = UtilsClient::getHisto<TH2F*>( meg12_[ism-1] );
        break;
      case 5:
      case 6:
      case 7:
      case 8:
        obj2f = 0;
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
        imgNameMEPnQualG16[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnQualG16[iCanvas-1];

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

      imgNameMEPnG01[iCanvas-1] = "";

      obj1d = 0;
      switch ( iCanvas ) {
        case 1:
          if ( i01_[ism-1] ) obj1d = i01_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 2:
          if ( i02_[ism-1] ) obj1d = i02_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 3:
          if ( i03_[ism-1] ) obj1d = i03_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 4:
          if ( i04_[ism-1] ) obj1d = i04_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 5:
        case 6:
        case 7:
        case 8:
          obj2f = 0;
          break;
        default:
          break;
      }

      if ( obj1d ) {

        meName = obj1d->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnG01[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnG01[iCanvas-1];

        cAmp->cd();
        gStyle->SetOptStat("euo");
        obj1d->SetStats(kTRUE);
//        if ( obj1d->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1d->SetMinimum(0.0);
        obj1d->Draw();
        cAmp->Update();
        cAmp->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

        delete obj1d;

      }

      imgNameMEPnG16[iCanvas-1] = "";

      obj1d = 0;
      switch ( iCanvas ) {
        case 1:
          if ( i09_[ism-1] ) obj1d = i09_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 2:
          if ( i10_[ism-1] ) obj1d = i10_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 3:
          if ( i11_[ism-1] ) obj1d = i11_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 4:
          if ( i12_[ism-1] ) obj1d = i12_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 5:
        case 6:
        case 7:
        case 8:
          obj2f = 0;
          break;
        default:
          break;
      }

      if ( obj1d ) {

        meName = obj1d->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnG16[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnG16[iCanvas-1];

        cAmp->cd();
        gStyle->SetOptStat("euo");
        obj1d->SetStats(kTRUE);
//        if ( obj1d->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1d->SetMinimum(0.0);
        obj1d->Draw();
        cAmp->Update();
        cAmp->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

        delete obj1d;

      }

      // Monitoring elements plots

      imgNameMEPnPedG01[iCanvas-1] = "";

      obj1d = 0;
      switch ( iCanvas ) {
        case 1:
          if ( i05_[ism-1] ) obj1d = i05_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 2:
          if ( i06_[ism-1] ) obj1d = i06_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 3:
          if ( i07_[ism-1] ) obj1d = i07_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 4:
          if ( i08_[ism-1] ) obj1d = i08_[ism-1]->ProjectionY("_py", 1, 1, "e");
          break;
        case 5:
        case 6:
        case 7:
        case 8:
          obj2f = 0;
          break;
        default:
          break;
      }

      if ( obj1d ) {

        meName = obj1d->GetName();

        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnPedG01[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnPedG01[iCanvas-1];

        cPed->cd();
        gStyle->SetOptStat("euo");
        obj1d->SetStats(kTRUE);
//        if ( obj1d->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1d->SetMinimum(0.0);
        obj1d->Draw();
        cPed->Update();
        cPed->SaveAs(imgName.c_str());
        gPad->SetLogy(0);

        delete obj1d;

      }
      
      
      imgNameMEPnRmsPedG01[iCanvas-1] = "";
      
      obj1f = 0;
      switch ( iCanvas ) {
      case 1:
	if ( mepnprms01_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms01_[ism-1]);
	break;
      case 2:
	if ( mepnprms02_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms02_[ism-1]);
	break;
      case 3:
	if ( mepnprms03_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms03_[ism-1]);
	break;
      case 4:
	if ( mepnprms04_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms04_[ism-1]);
	break;
      case 5:
      case 6:
      case 7:
      case 8:
	obj2f = 0;
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
	imgNameMEPnRmsPedG01[iCanvas-1] = meName + ".png";
	imgName = htmlDir + imgNameMEPnRmsPedG01[iCanvas-1];
  	
	cPed->cd();
	gStyle->SetOptStat("euomr");
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
      
      
      imgNameMEPnRmsPedG16[iCanvas-1] = "";
      
      obj1f = 0;
      switch ( iCanvas ) {
      case 1:
	if ( mepnprms05_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms05_[ism-1]);
	break;
      case 2:
	if ( mepnprms06_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms06_[ism-1]);
	break;
      case 3:
	if ( mepnprms07_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms07_[ism-1]);
	break;
      case 4:
	if ( mepnprms08_[ism-1] ) obj1f =  UtilsClient::getHisto<TH1F*>(mepnprms08_[ism-1]);
	break;
      case 5:
      case 6:
      case 7:
      case 8:
	obj2f = 0;
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
	imgNameMEPnRmsPedG16[iCanvas-1] = meName + ".png";
	imgName = htmlDir + imgNameMEPnRmsPedG16[iCanvas-1];
  	
	cPed->cd();
	gStyle->SetOptStat("euomr");
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
      
      imgNameMEPnPedG16[iCanvas-1] = "";
      
      obj1d = 0;
      switch ( iCanvas ) {
      case 1:
	if ( i13_[ism-1] ) obj1d = i13_[ism-1]->ProjectionY("_py", 1, 1, "e");
	break;
      case 2:
	if ( i14_[ism-1] ) obj1d = i14_[ism-1]->ProjectionY("_py", 1, 1, "e");
	break;
      case 3:
	if ( i15_[ism-1] ) obj1d = i15_[ism-1]->ProjectionY("_py", 1, 1, "e");
	break;
      case 4:
	if ( i16_[ism-1] ) obj1d = i16_[ism-1]->ProjectionY("_py", 1, 1, "e");
	break;
      case 5:
      case 6:
      case 7:
      case 8:
	obj2f = 0;
	break;
      default:
	break;
      }
      
      if ( obj1d ) {
	
        meName = obj1d->GetName();
	
        for ( unsigned int i = 0; i < meName.size(); i++ ) {
          if ( meName.substr(i, 1) == " " )  {
            meName.replace(i, 1 ,"_" );
          }
        }
        imgNameMEPnPedG16[iCanvas-1] = meName + ".png";
        imgName = htmlDir + imgNameMEPnPedG16[iCanvas-1];
	
        cPed->cd();
        gStyle->SetOptStat("euo");
        obj1d->SetStats(kTRUE);
//        if ( obj1d->GetMaximum(histMax) > 0. ) {
//          gPad->SetLogy(1);
//        } else {
//          gPad->SetLogy(0);
//        }
        obj1d->SetMinimum(0.0);
        obj1d->Draw();
        cPed->Update();
        cPed->SaveAs(imgName.c_str());
        gPad->SetLogy(0);
	
        delete obj1d;
	
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

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameQual[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameQual[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameAmp[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameAmp[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameAmpoPN[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameAmpoPN[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameAmp[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameAmp[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameAmpoPN[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameAmpoPN[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameTim[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTim[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameShape[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameShape[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameTimav[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTimav[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameTimrms[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTimrms[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameTim[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTim[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameShape[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameShape[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameTimav[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTimav[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameTimrms[4+iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameTimrms[4+iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      htmlFile << "<td colspan=\"2\">Laser " << iCanvas << "</td>" << endl;

    }

    htmlFile << "</tr>" << endl;
    htmlFile << "</table>" << endl;
    htmlFile << "<br>" << endl;

    htmlFile << "<table border=\"0\" cellspacing=\"0\" " << endl;
    htmlFile << "cellpadding=\"10\" align=\"center\"> " << endl;
    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameMEPnQualG01[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPnQualG01[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameMEPnPedG01[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnPedG01[iCanvas-1] << "\"></td>" << endl;
      else
       htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPnRmsPedG01[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnRmsPedG01[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPnG01[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnG01[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      htmlFile << "<td colspan=\"1\"> </td> <td colspan=\"1\">Laser " << iCanvas << " - PN Gain 1</td> <td colspan=\"1\">" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameMEPnQualG16[iCanvas-1].size() != 0 )
        htmlFile << "<td colspan=\"2\"><img src=\"" << imgNameMEPnQualG16[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td colspan=\"2\"><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "<tr>" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      if ( imgNameMEPnPedG16[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnPedG16[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPnRmsPedG16[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnRmsPedG16[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

      if ( imgNameMEPnG16[iCanvas-1].size() != 0 )
        htmlFile << "<td><img src=\"" << imgNameMEPnG16[iCanvas-1] << "\"></td>" << endl;
      else
        htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;

    }

    htmlFile << "<tr align=\"center\">" << endl;

    for ( int iCanvas = 1 ; iCanvas <= 4 ; iCanvas++ ) {

      // skip unused wavelengths
      if ( iCanvas == 2 || iCanvas == 3 ) continue;

      htmlFile << "<td colspan=\"1\"> </td> <td colspan=\"1\">Laser " << iCanvas << " - PN Gain 16</td> <td colspan=\"1\"> </td>" << endl;

    }

    htmlFile << "</tr>" << endl;

    htmlFile << "</tr>" << endl;
    htmlFile << "</table>" << endl;

    htmlFile << "<br>" << endl;

  }

  delete cQual;
  delete cQualPN;
  delete cAmp;
  delete cTim;
  delete cTimav;
  delete cTimrms;
  delete cShape;
  delete cAmpoPN;
  delete cPed;

  // html page footer
  htmlFile << "</body> " << endl;
  htmlFile << "</html> " << endl;

  htmlFile.close();

}

