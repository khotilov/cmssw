// -*- C++ -*-
//
// Package:   EcalCreateTimeCalibrations
// Class:     EcalCreateTimeCalibrations
//
/**\class EcalCreateTimeCalibrations EcalCreateTimeCalibrations.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
 */
//
// Authors:                              Seth Cooper (Minnesota)
//         Created:  Tu Apr 26  10:46:22 CEST 2011
// $Id: $
//
//

#include "CalibCalorimetry/EcalTiming/plugins/EcalCreateTimeCalibrations.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondTools/Ecal/interface/EcalTimeCalibConstantsXMLTranslator.h"
#include "CondTools/Ecal/interface/EcalTimeCalibErrorsXMLTranslator.h"
#include "CondTools/Ecal/interface/EcalTimeOffsetXMLTranslator.h"
#include "CondTools/Ecal/interface/EcalCondHeader.h"

#include "CalibCalorimetry/EcalTiming/interface/EcalCrystalTimingCalibration.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include <fstream>


EcalCreateTimeCalibrations::EcalCreateTimeCalibrations(const edm::ParameterSet& ps) :
  inputFiles_ (ps.getParameter<std::vector<std::string> >("InputFileNames")),
  fileName_ (ps.getParameter<std::string>("FileNameStart")),
  numTotalCrys_ (EBDetId::kSizeForDenseIndexing+EEDetId::kSizeForDenseIndexing),
  subtractDBcalibs_ (ps.getParameter<bool>("SubtractDBcalibs")),
  inBxs_ (ps.getParameter<std::string>("BxIncludeExclude")),
  inOrbits_ (ps.getParameter<std::string>("OrbitIncludeExclude")),
  inTrig_ (ps.getParameter<std::string>("TriggerBitIncludeExclude")),
  inTTrig_ (ps.getParameter<std::string>("TechTriggerBitIncludeExclude")),
  inLumis_ (ps.getParameter<std::string>("LumiIncludeExclude")),
  inRuns_ (ps.getParameter<std::string>("RunIncludeExclude")),
  avgTimeMin_ (ps.getParameter<double>("AvgTimeMin")),
  avgTimeMax_ (ps.getParameter<double>("AvgTimeMax")),
  minHitAmpEB_ (ps.getParameter<double>("MinHitAmpEB")),
  minHitAmpEE_ (ps.getParameter<double>("MinHitAmpEE")),
  maxSwissCrossNoise_ (ps.getParameter<double>("MaxSwissCross")),
  maxHitTimeEB_ (ps.getParameter<double>("MaxHitTimeEB")),
  minHitTimeEB_ (ps.getParameter<double>("MinHitTimeEB")),
  maxHitTimeEE_ (ps.getParameter<double>("MaxHitTimeEE")),
  minHitTimeEE_ (ps.getParameter<double>("MinHitTimeEE"))
{

  edm::Service<TFileService> fileService_;

  // input tree construction
  myInputTree_ = new TChain ("EcalTimeAnalysis") ;
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=inputFiles_.begin(); file_itr!=inputFiles_.end(); file_itr++)
    myInputTree_->Add( (*file_itr).c_str() );

  if(!myInputTree_)
  {
    edm::LogError("EcalCreateTimeCalibrations") << "Couldn't find tree EcalTimeAnalysis";
    produce_ = false;
    return;
  }

  edm::LogInfo("EcalCreateTimeCalibrations") << "Running with options: "
    << "avgTimeMin: " << avgTimeMin_ << " avgTimeMax: " << avgTimeMax_
    << " minAmpEB: " << minHitAmpEB_ << " minHitAmpEE: " << minHitAmpEE_
    << " maxSwissCrossNoise (EB): " << maxSwissCrossNoise_
    << " maxHitTimeEB: " << maxHitTimeEB_ << " minHitTimeEB: " << minHitTimeEB_
    << " maxHitTimeEE: " << maxHitTimeEE_ << " minHitTimeEB: " << minHitTimeEE_
    << " inTrig: " << inTrig_ << " inTTrig: " << inTTrig_ << " inLumi: " << inLumis_
    << " inBxs: " << inBxs_ << " inRuns: " << inRuns_ << " inOrbits: " << inOrbits_;

  setBranchAddresses(myInputTree_,treeVars_);

  //recall: string inBxs, inOrbits, inTrig, inTTrig, inLumi, inRuns;
  genIncludeExcludeVectors(inBxs_,bxIncludeVector,bxExcludeVector);
  genIncludeExcludeVectors(inOrbits_,orbitIncludeVector,orbitExcludeVector);
  genIncludeExcludeVectors(inTrig_,trigIncludeVector,trigExcludeVector);
  genIncludeExcludeVectors(inTTrig_,ttrigIncludeVector,ttrigExcludeVector);
  genIncludeExcludeVectors(inLumis_,lumiIncludeVector,lumiExcludeVector);
  genIncludeExcludeVectors(inRuns_,runIncludeVector,runExcludeVector);

  produce_ = true;
  ttEtaBins_ = {-85, -80, -75, -70, -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86 };
  initEBHists(fileService_);
  initEEHists(fileService_);
}

EcalCreateTimeCalibrations::~EcalCreateTimeCalibrations()
{
}

void
EcalCreateTimeCalibrations::analyze(edm::Event const& evt, edm::EventSetup const& es)
{
  if(!produce_)
    return;

  EcalCrystalTimingCalibration* ebCryCalibs[EBDetId::kSizeForDenseIndexing];
  //XXX: Making calibs with weighted/unweighted mean
  for(int i=0; i < EBDetId::kSizeForDenseIndexing; ++i)
    ebCryCalibs[i] = new EcalCrystalTimingCalibration(); //use weighted mean!
    //ebCryCalibs[i] = new EcalCrystalTimingCalibration(false); //don't use weighted mean!
  EcalCrystalTimingCalibration* eeCryCalibs[EEDetId::kSizeForDenseIndexing];
  //XXX: Making calibs with weighted/unweighted mean
  for(int i=0; i < EEDetId::kSizeForDenseIndexing; ++i)
    eeCryCalibs[i] = new EcalCrystalTimingCalibration(); //use weighted mean!
    //eeCryCalibs[i] = new EcalCrystalTimingCalibration(false); //don't use weighted mean!

  // Loop over the TTree
  int numEventsUsed = 0;
  int nEntries = myInputTree_->GetEntries();
  edm::LogInfo("EcalCreateTimeCalibrations") << "Begin loop over TTree";
  for(int entry = 0; entry < nEntries; ++entry)
  {
    myInputTree_->GetEntry(entry);

    // Loop once to calculate average event time -- use all crys in EB+EE clusters
    float sumTime = 0;
    int numCrys = 0;
    for(int bCluster=0; bCluster < treeVars_.nClusters; bCluster++)
    {
      for(int cryInBC=0; cryInBC < treeVars_.nXtalsInCluster[bCluster]; cryInBC++)
      {
        sumTime += treeVars_.xtalInBCTime[bCluster][cryInBC];
        numCrys++;
      }
    }
    //debug
    //cout << "Number of EB crys in event: " << numEBcrys << endl;

    //XXX: Event cuts
    if(sumTime/numCrys > avgTimeMax_ || sumTime/numCrys < avgTimeMin_)
    {
      //cout << "Average event time: " << sumTime/numCrys  << " so event rejected." << endl;
      continue;
    }
    // check BX, orbit, lumi, run, L1 tech/phys triggers
    bool keepEvent = includeEvent(treeVars_.bx,bxIncludeVector,bxExcludeVector)
      && includeEvent(treeVars_.orbit,orbitIncludeVector,orbitExcludeVector)
      && includeEvent(treeVars_.lumiSection,lumiIncludeVector,lumiExcludeVector)
      && includeEvent(treeVars_.runId,runIncludeVector,runExcludeVector)
      && includeEvent(treeVars_.l1ActiveTriggers,
          treeVars_.l1NActiveTriggers,trigIncludeVector,trigExcludeVector)
      && includeEvent(treeVars_.l1ActiveTechTriggers,
          treeVars_.l1NActiveTechTriggers,ttrigIncludeVector,ttrigExcludeVector);
    if(!keepEvent)
      continue;
      
    numEventsUsed++;

    // Loop over the clusters in the event
    for(int bCluster=0; bCluster < treeVars_.nClusters; bCluster++)
    {
      for(int cryInBC=0; cryInBC < treeVars_.nXtalsInCluster[bCluster]; cryInBC++)
      {
        int hashedIndex = treeVars_.xtalInBCHashedIndex[bCluster][cryInBC];
        float cryTime = treeVars_.xtalInBCTime[bCluster][cryInBC];
        float cryTimeError = treeVars_.xtalInBCTimeErr[bCluster][cryInBC];
        float cryAmp = treeVars_.xtalInBCAmplitudeADC[bCluster][cryInBC];
        //SIC FEB 14,16 2011 - removed E/E9 cut
        //
        float crySwissCrossNoise = treeVars_.xtalInBCSwissCross[bCluster][cryInBC];

        if(treeVars_.xtalInBCIEta[bCluster][0] != -999999)
        {
          EBDetId det = EBDetId::unhashIndex(hashedIndex);
          if(det==EBDetId()) // make sure DetId is valid
            continue;

          int ieta = det.ieta();
          int iphi = det.iphi();

          // RecHit cuts
          bool keepHit = cryAmp >= minHitAmpEB_
            && crySwissCrossNoise < maxSwissCrossNoise_
            && cryTime > minHitTimeEB_
            && cryTime < maxHitTimeEB_;
            if(!keepHit)
              continue;

          //cout << "DEBUG: " << hashedIndex << " cryTime: " << cryTime << " cryTimeError: " << cryTimeError << " cryAmp: " << cryAmp << endl;
          //FIXME
          cryTimeError = 1;
          ebCryCalibs[hashedIndex]->insertEvent(cryAmp,cryTime,cryTimeError,false);
          //SIC Use when we don't have time_error available
          //ebCryCalibs[hashedIndex]->insertEvent(cryAmp,cryTime,35/(cryAmp/1.2),false);
          ampProfileEB_->Fill(hashedIndex,cryAmp);
          ampProfileMapEB_->Fill(iphi,ieta,cryAmp);
          //if(cryTime > 33 || cryTime < -33)
          //  cout << "Crystal: " << det << " event time is over/underflow: " << cryTime << endl;
        }
        else
        {
          EEDetId det = EEDetId::unhashIndex(hashedIndex);
          if(det==EEDetId()) // make sure DetId is valid
            continue;

          int ix = det.ix();
          int iy = det.iy();

          //XXX: RecHit cuts
          bool keepHit = cryAmp >= minHitAmpEE_
            && cryTime > minHitTimeEE_
            && cryTime < maxHitTimeEE_;
          if(!keepHit)
            continue;

          //std::cout << "STUPID DEBUG: EE CRY " << hashedIndex << " cryTime: " << cryTime << " cryTimeError: " << cryTimeError << " cryAmp: " << cryAmp << std::endl;
          //FIXME
          cryTimeError = 1;
          eeCryCalibs[hashedIndex]->insertEvent(cryAmp,cryTime,cryTimeError,false);

          //SIC Use when we don't have time_error available
          //eeCryCalibs[hashedIndex]->insertEvent(cryAmp,cryTime,35/(cryAmp/1.2),false);
          if(det.zside() < 0)
          {
            ampProfileEEM_->Fill(hashedIndex,cryAmp);
            ampProfileMapEEM_->Fill(ix,iy,cryAmp);
          }
          else
          {
            ampProfileEEP_->Fill(hashedIndex,cryAmp);
            ampProfileMapEEP_->Fill(ix,iy,cryAmp);
          }
        }
      }
    }
  }


  //create output text file
  ofstream fileStream;
  std::string fileName1 = fileName_+".calibs.txt";
  fileStream.open(fileName1.c_str());
  if(!fileStream.good() || !fileStream.is_open())
  {
    edm::LogError("EcalCreateTimeCalibrations") << "Couldn't open text file.";
    return;
  }
  //create problem channels text file
  ofstream fileStreamProb;
  std::string fileName2 = fileName_+".problems.txt";
  fileStreamProb.open(fileName2.c_str());
  if(!fileStreamProb.good() || !fileStreamProb.is_open())
  {
    edm::LogError("EcalCreateTimeCalibrations") << "Couldn't open text file.";
    return;
  }

  // Create calibration container objects
  EcalTimeCalibConstants timeCalibConstants;
  EcalTimeCalibErrors timeCalibErrors;
  EcalTimeOffsetConstant timeOffsetConstant;

  edm::LogInfo("EcalCreateTimeCalibrations") << "Using " << numEventsUsed << " out of " << nEntries << " in the tree.";
  edm::LogInfo("EcalCreateTimeCalibrations") << "Creating calibs...";
  float cryCalibAvg = 0;
  int numCrysCalibrated = 0;
  std::vector<int> hashesToCalibrateToZeroEB;
  std::vector<int> hashesToCalibrateToZeroEE;
  std::vector<int> hashesToCalibrateNormallyEB;
  std::vector<int> hashesToCalibrateNormallyEE;
  //Loop over all the crys
  for(int hashedIndex=0; hashedIndex < numTotalCrys_; ++hashedIndex)
  {
    EcalCrystalTimingCalibration cryCalib;
    int ieta = 0;
    int iphi = 0;
    int x = 0;
    int y = 0;
    int zside = 0;
    bool isEB = true;
    if(hashedIndex >= EBDetId::kSizeForDenseIndexing)
      isEB = false;

    if(isEB)
    {
      EBDetId det = EBDetId::unhashIndex(hashedIndex);
      if(det==EBDetId())
        continue;
      ieta = det.ieta();
      iphi = det.iphi();
      cryCalib = *(ebCryCalibs[hashedIndex]);
    }
    else
    {
      EEDetId det = EEDetId::unhashIndex(hashedIndex-EBDetId::kSizeForDenseIndexing);
      if(det==EEDetId())
        continue;
      x = det.ix();
      y = det.iy();
      zside = det.zside();
      cryCalib = *(eeCryCalibs[hashedIndex-EBDetId::kSizeForDenseIndexing]);
    }

    //XXX: Filter events at default 0.5*meanE threshold
    cryCalib.filterOutliers();
    
    //numPointsErasedHist_->Fill(numPointsErased);
    //chiSquaredTotalHist_->Fill(cryCalib.totalChi2);
    double p1 = cryCalib.mean;
    double p1err = cryCalib.meanE;
    
    //Write cryTimingHists
    std::vector<EcalTimingEvent> times = cryCalib.timingEvents;
    for(std::vector<EcalTimingEvent>::const_iterator timeItr = times.begin();
        timeItr != times.end(); ++timeItr)
    {
      float weight = 1/((timeItr->sigmaTime)*(timeItr->sigmaTime));
      if(isEB)
      {
        cryTimingHistsEB_[hashedIndex]->Fill(timeItr->time,weight);
        //expectedStatPresHistEB_->Fill(sqrt(1/expectedPresSumEB));
        //expectedStatPresVsObservedMeanErrHistEB_->Fill(sigmaM,sqrt(1/expectedPresSumEB));
      }
      else
      {
        if(zside < 0)
        {
          float weight = 1/((timeItr->sigmaTime)*(timeItr->sigmaTime));
          cryTimingHistsEEM_[x-1][y-1]->Fill(timeItr->time,weight);
        }
        else
        {
          float weight = 1/((timeItr->sigmaTime)*(timeItr->sigmaTime));
          cryTimingHistsEEP_[x-1][y-1]->Fill(timeItr->time,weight);
        }
      }
    }

    if(isEB)
    {
      //cryDirEB->cd();
      cryTimingHistsEB_[hashedIndex]->Write();
      hitsPerCryHistEB_->SetBinContent(hashedIndex+1,cryCalib.timingEvents.size());
      hitsPerCryMapEB_->Fill(iphi,ieta,cryCalib.timingEvents.size());
    }
    else
    {
      if(zside < 0)
      {
        //cryDirEEM->cd();
        cryTimingHistsEEM_[x-1][y-1]->Write();
        hitsPerCryHistEEM_->SetBinContent(hashedIndex+1,cryCalib.timingEvents.size());
        hitsPerCryMapEEM_->Fill(x,y,cryCalib.timingEvents.size());
      }
      else
      {
        //cryDirEEP->cd();
        cryTimingHistsEEP_[x-1][y-1]->Write();
        hitsPerCryHistEEP_->SetBinContent(hashedIndex+1,cryCalib.timingEvents.size());
        hitsPerCryMapEEP_->Fill(x,y,cryCalib.timingEvents.size());
      }

    }

    if(p1err < 0.5 && p1err > 0)
    {
      cryCalibAvg+=p1;
      if(isEB)
      {
        hashesToCalibrateNormallyEB.push_back(hashedIndex);

        calibHistEB_->Fill(p1);
        //calibMapEEMFlip_->Fill(y-85,x+1,p1);
        calibMapEB_->Fill(iphi,ieta,p1);
        calibTTMapEB_->Fill(iphi,ieta,p1);
        //calibMapEEMPhase_->Fill(x+1,y-85,p1/25-floor(p1/25));
        //errorOnMeanVsNumEvtsHist_->Fill(times.size(),p1err);
      }
      else
      {
        hashesToCalibrateNormallyEE.push_back(hashedIndex-EBDetId::kSizeForDenseIndexing);
        if(zside < 0)
        {
          calibHistEEM_->Fill(p1);
          //calibMapEEMFlip_->Fill(y-85,x+1,p1);
          calibMapEEM_->Fill(x,y,p1);
          //calibMapEEMPhase_->Fill(x+1,y-85,p1/25-floor(p1/25));
          //errorOnMeanVsNumEvtsHistEE_->Fill(times.size(),p1err);
        }
        else
        {
          calibHistEEP_->Fill(p1);
          //calibMapEEPFlip_->Fill(y-85,x+1,p1);
          calibMapEEP_->Fill(x,y,p1);
          //calibMapEEPPhase_->Fill(x+1,y-85,p1/25-floor(p1/25));
          //errorOnMeanVsNumEvtsHist_->Fill(times.size(),p1err);
        }
      }
    }
    else
    {
      //std::cout << "Cry: " << ieta <<", " << iphi << ", hash: " << itr->first
      //  << "\t Calib: " << p1 << "\t Error: " << p1err << std::endl;
      if(isEB)
      {
        hashesToCalibrateToZeroEB.push_back(hashedIndex);
        fileStreamProb << "EB Cry ( " << cryCalib.timingEvents.size() << " events) was calibrated to avg: " << ieta <<", " << iphi << ", hash: " << hashedIndex
          << "\t Calib: " << p1 << "\t Error: " << p1err << std::endl;
      }
      else
      {
        hashesToCalibrateToZeroEE.push_back(hashedIndex-61200);
        fileStreamProb << "EE Cry ( " << cryCalib.timingEvents.size() << " events) was calibrated to avg: " << ieta <<", " << iphi << ", hash: " << hashedIndex-61200
          << "\t Calib: " << p1 << "\t Error: " << p1err << std::endl;
      }
    }

    // Fill other hists regardless
    if(isEB)
    {
      //calibsVsErrorsEB_->Fill(p1err, p1 > 0 ? p1 : -1*p1);
      calibErrorHistEB_->Fill(p1err);
      calibErrorMapEB_->Fill(iphi,ieta,p1err);
      sigmaHistEB_->Fill(cryCalib.stdDev);
      sigmaMapEB_->Fill(iphi,ieta,cryCalib.stdDev);
    }
    else
    {
      sigmaHistEE_->Fill(cryCalib.stdDev);
      if(zside < 0)
      {
        //calibsVsErrorsEEM->Fill(p1err, p1 > 0 ? p1 : -1*p1);
        calibErrorHistEEM_->Fill(p1err);
        calibErrorMapEEM_->Fill(x,y,p1err);
        sigmaMapEEM_->Fill(x,y,cryCalib.stdDev);
      }
      else
      {
        //calibsVsErrorsEEP->Fill(p1err, p1 > 0 ? p1 : -1*p1);
        calibErrorHistEEP_->Fill(p1err);
        calibErrorMapEEP_->Fill(x,y,p1err);
        sigmaMapEEP_->Fill(x,y,cryCalib.stdDev);
      }
    }
  }
  fileStreamProb.close();

  edm::LogInfo("EcalCreateTimeCalibrations") << "Calibrating " << hashesToCalibrateNormallyEB.size()
    << " EB crys normally and " << hashesToCalibrateToZeroEB.size() << " EB crys to zero.";
  edm::LogInfo("EcalCreateTimeCalibrations") << "Calibrating " << hashesToCalibrateNormallyEE.size()
    << " EE crys normally and " << hashesToCalibrateToZeroEE.size() << " EE crys to zero.";

  // Calc average
  numCrysCalibrated = hashesToCalibrateNormallyEB.size()+hashesToCalibrateNormallyEE.size();
  if(numCrysCalibrated > 0)
    cryCalibAvg/=numCrysCalibrated;
  //cryCalibAvg-= 2.0833; // Global phase shift

  //Loop over all the crys to calibrate normally -- EB
  for(std::vector<int>::const_iterator hashItr = hashesToCalibrateNormallyEB.begin();
      hashItr != hashesToCalibrateNormallyEB.end(); ++hashItr)
  {
    EBDetId det = EBDetId::unhashIndex(*hashItr);
    if(det==EBDetId())
      continue;
    EcalCrystalTimingCalibration cryCalib = *(ebCryCalibs[*hashItr]);
    // Make timing calibs
    double p1 = cryCalib.mean-cryCalibAvg;
    double p1err = cryCalib.meanE;
    // Make it so we can add calib to reco time
    p1*=-1;
    fileStream << "EB\t" << *hashItr << "\t" << p1 << "\t\t" << p1err << std::endl;
    //Store in timeCalibration container
    EcalTimeCalibConstant tcConstant = p1;
    EcalTimeCalibError tcError = p1err;
    uint32_t rawId = EBDetId::unhashIndex(*hashItr);
    timeCalibConstants[rawId] = tcConstant;
    timeCalibErrors[rawId] = tcError;
  }
  //Loop over all the crys to calibrate normally -- EE
  for(std::vector<int>::const_iterator hashItr = hashesToCalibrateNormallyEE.begin();
      hashItr != hashesToCalibrateNormallyEE.end(); ++hashItr)
  {
    EEDetId det = EEDetId::unhashIndex(*hashItr);
    if(det==EEDetId())
      continue;
    EcalCrystalTimingCalibration cryCalib = *(eeCryCalibs[*hashItr]);
    // Make timing calibs
    double p1 = cryCalib.mean-cryCalibAvg;
    double p1err = cryCalib.meanE;
    // Make it so we can add calib to reco time
    p1*=-1;
    fileStream << "EE\t" << *hashItr << "\t" << p1 << "\t\t" << p1err << std::endl;
    //Store in timeCalibration container
    EcalTimeCalibConstant tcConstant = p1;
    EcalTimeCalibError tcError = p1err;
    uint32_t rawId = EEDetId::unhashIndex(*hashItr);
    timeCalibConstants[rawId] = tcConstant;
    timeCalibErrors[rawId] = tcError;
  }
  
  fileStream.close();
  // calibrate uncalibratable crys -- EB
  for(std::vector<int>::const_iterator hashItr = hashesToCalibrateToZeroEB.begin();
      hashItr != hashesToCalibrateToZeroEB.end(); ++hashItr)
  {
    EBDetId det = EBDetId::unhashIndex(*hashItr);
    if(det==EBDetId())
      continue;
    //Store in timeCalibration container
    double p1 = 0;
    double p1err = 999;
    EcalTimeCalibConstant tcConstant = p1;
    EcalTimeCalibError tcError = p1err;
    uint32_t rawId = det.rawId();
    timeCalibConstants[rawId] = tcConstant;
    //timeCalibConstants.insert(std::pair<uint32_t,EcalTimeCalibConstant>(rawId,tcConstant));
    timeCalibErrors[rawId] = tcError;
  }
  // calibrate uncalibratable crys -- EE
  for(std::vector<int>::const_iterator hashItr = hashesToCalibrateToZeroEE.begin();
      hashItr != hashesToCalibrateToZeroEE.end(); ++hashItr)
  {
    EEDetId det = EEDetId::unhashIndex(*hashItr);
    if(det==EEDetId())
      continue;
    //Store in timeCalibration container
    EcalTimeCalibConstant tcConstant = 0;
    EcalTimeCalibError tcError = 999;
    uint32_t rawId = EEDetId::unhashIndex(*hashItr);
    timeCalibConstants[rawId] = tcConstant;
    timeCalibErrors[rawId] = tcError;
  }

  // global time offset
  timeOffsetConstant.setEBValue(-1*cryCalibAvg);
  timeOffsetConstant.setEEValue(-1*cryCalibAvg);

  //TODO at the end
  if(subtractDBcalibs_)
  {
    set(es);
    // Loop over all crys
    for(int hashedIndex=0; hashedIndex<numTotalCrys_; ++hashedIndex)
    {
      DetId det = 0;
      bool isEB = true;
      if(hashedIndex >= EBDetId::kSizeForDenseIndexing)
        isEB = false;

      if(isEB)
      {
        det = EBDetId::unhashIndex(hashedIndex);
        if(det==EBDetId())
          continue;
      }
      else
      {
        det = EEDetId::unhashIndex(hashedIndex-EBDetId::kSizeForDenseIndexing);
        if(det==EEDetId())
          continue;
      }

      // get orig time calibration coefficient
      const EcalTimeCalibConstantMap & origTimeMap = origTimeCalibConstHandle->getMap();
      EcalTimeCalibConstant itimeconstOrig = 0;
      EcalTimeCalibConstantMap::const_iterator itimeItrOrig = origTimeMap.find(det);
      if( itimeItrOrig!=origTimeMap.end() )
      {
        itimeconstOrig = (*itimeItrOrig);
      }
      // if no time constant in DB, use zero
      timeCalibConstants[det.rawId()]-=itimeconstOrig;
    }

  }


  //Write XML files
  edm::LogInfo("EcalCreateTimeCalibrations") << "Writing XML files.";
  EcalCondHeader header;
  header.method_="testmethod";
  header.version_="testversion";
  header.datasource_="testdata";
  header.since_=123;
  header.tag_="testtag";
  header.date_="Mar 24 1973";
  std::string timeCalibFile = "EcalTimeCalibs.xml";
  std::string timeCalibErrFile = "EcalTimeCalibErrors.xml";
  std::string timeOffsetFile = "EcalTimeOffset.xml";
  //NOTE: Hack should now be unnecessary...
  // Hack to prevent seg fault
  //EcalTimeCalibConstant tcConstant = 0;
  //EcalTimeCalibError tcError = 0;
  //uint32_t rawId = EEDetId::unhashIndex(0);
  //timeCalibConstants[rawId] = tcConstant;
  //timeCalibErrors[rawId] = tcError;
  // End hack
  EcalTimeCalibConstantsXMLTranslator::writeXML(timeCalibFile,header,timeCalibConstants);
  EcalTimeCalibErrorsXMLTranslator::writeXML(timeCalibErrFile,header,timeCalibErrors);
  EcalTimeOffsetXMLTranslator::writeXML(timeOffsetFile,header,timeOffsetConstant);

  //Move empty bins out of the way
  int nxbins = calibMapEEM_->GetNbinsX();
  int nybins = calibMapEEM_->GetNbinsY();
  for(int i=0;i<=(nxbins+2)*(nybins+2); ++i)
  {
    double binentsM = calibMapEEM_->GetBinContent(i);
    if(binentsM==0)
    {
      calibMapEEM_->SetBinContent(i,-1000);
    }
    double binentsP = calibMapEEP_->GetBinContent(i);
    if(binentsP==0)
    {
      calibMapEEP_->SetBinContent(i,-1000);
    }
  }

}

void EcalCreateTimeCalibrations::set(edm::EventSetup const& eventSetup)
{
  eventSetup.get<EcalTimeCalibConstantsRcd>().get(origTimeCalibConstHandle);
}

void EcalCreateTimeCalibrations::beginRun(edm::EventSetup const& eventSetup)
{
}

void EcalCreateTimeCalibrations::beginJob()
{
}

void EcalCreateTimeCalibrations::endJob()
{
}

void EcalCreateTimeCalibrations::initEBHists(edm::Service<TFileService>& fileService_)
{
  calibHistEB_ = fileService_->make<TH1F>("timingCalibDiffEB","timingCalib diff EB [ns]",2000,-100,100);
  calibErrorHistEB_ = fileService_->make<TH1F>("calibDiffErrorEB","timingCalibDiffError EB [ns]",500,0,5);
  calibHistEB_->Sumw2();
  calibErrorHistEB_->Sumw2();
  calibsVsErrors_ = fileService_->make<TH2F>("timingCalibDiffsAndErrors","TimingCalibDiffs vs. errors [ns]",500,0,5,100,0,10);
  calibsVsErrors_->Sumw2();
  expectedStatPresHistEB_ = fileService_->make<TH1F>("expectedStatPresEB","Avg. expected statistical precision EB [ns], all crys",200,0,2);
  expectedStatPresVsObservedMeanErrHistEB_ = fileService_->make<TH2F>("expectedStatPresVsObsEB","Expected stat. pres. vs. obs. error on mean each event EB [ns]",200,0,2,200,0,2);
  expectedStatPresEachEventHistEB_ = fileService_->make<TH1F>("expectedStatPresSingleEventEB","Expected stat. pres. each event EB [ns]",200,0,2);
  errorOnMeanVsNumEvtsHist_ = fileService_->make<TH2F>("errorOnMeanVsNumEvts","Error_on_mean vs. number of events",50,0,50,200,0,2);
  errorOnMeanVsNumEvtsHist_->Sumw2();
  hitsPerCryHistEB_ = fileService_->make<TH1F>("hitsPerCryEB","Hits used in each crystal;hashedIndex",61200,0,61200);
  hitsPerCryMapEB_ = fileService_->make<TH2F>("hitsPerCryMapEB","Hits used in each crystal;i#phi;i#eta",360,1.,361.,172,-86,86);
  ampProfileMapEB_ = fileService_->make<TProfile2D>("ampProfileMapEB","amp profile map [ADC];i#phi;i#eta",360,1.,361.,172,-86,86);
  ampProfileEB_ = fileService_->make<TProfile>("ampProfileEB","Average amplitude in cry [ADC];hashedIndex",61200,0,61200);
  sigmaHistEB_ = fileService_->make<TH1F>("sigmaCalibsEB"," Sigma of calib distributions EB [ns]",100,0,1);
  for (int i = 0; i < 79; ++i)
  {
    timingBins_[i]=-7.+double(i)*14./78.;
    if (i<73)
    {
      ttPhiBins_[i]=1+5*i;
      if ( i < 19) 
      {
        modPhiBins_[i]=1+20*i;
        if (i < 11)
        {
          highEBins_[i]=10.+double(i)*20.;
        }
      }
    }
  } 
  calibMapEB_ = fileService_->make<TH2F>("calibDiffMapEB","time calib diff map EB [ns];i#phi;i#eta",360,1.,361.,172,-86,86);
  calibMapEB_->Sumw2();
  sigmaMapEB_ = fileService_->make<TH2F>("sigmaDiffMapEB","Sigma of time calib diff map EB [ns];i#phi;i#eta",360,1.,361.,172,-86,86);
  calibErrorMapEB_ = fileService_->make<TH2F>("calibDiffErrorMapEB","Error of time calib diff map EB [ns];i#phi;i#eta",360,1.,361.,172,-86,86);
  calibTTMapEB_ = fileService_->make<TProfile2D>("calibDiffTTMapEB","time calib diff map EB (TT) [ns];i#phi;i#eta",360/5,ttPhiBins_,35, ttEtaBins_);
  TFileDirectory cryDirEB = fileService_->mkdir("crystalTimingHistsEB");
  EBDetId det;
  for(int hi=0; hi < 61200; ++hi)
  {
    det = EBDetId::unhashIndex(hi);
    if(det==EBDetId())
      continue;
    std::string histname = "EB_cryTiming_ieta";
    histname+=intToString(det.ieta());
    histname+="_iphi";
    histname+=intToString(det.iphi());
    cryTimingHistsEB_[hi] = cryDirEB.make<TH1F>(histname.c_str(),histname.c_str(),200,-10,10);
    cryTimingHistsEB_[hi]->Sumw2();
  }
}

void EcalCreateTimeCalibrations::initEEHists(edm::Service<TFileService>& fileService_)
{
  calibHistEE_ = fileService_->make<TH1F>("timingCalibDiffsEE","timingCalibDiffs EE [ns]",2000,-100,100);
  calibErrorHistEE_ = fileService_->make<TH1F>("timingCalibDiffErrorEE","timingCalibDiffError EE [ns]",500,0,5);
  calibHistEE_->Sumw2();
  calibErrorHistEE_->Sumw2();
  calibsVsErrorsEE_ = fileService_->make<TH2F>("timingCalibDiffsAndErrors","TimingCalibDiffs vs. errors [ns]",500,0,5,100,0,10);
  calibsVsErrorsEE_->Sumw2();
  //calibMapEE_ = fileService_->make<TH2F>("calibMapEE","time calib diff map EE",360,1,361,170,-86,86);
  //calibMapEEFlip_ = fileService_->make<TH2F>("calibMapEEFlip","time calib diff map EE",170,-86,86,360,1,361);
  //calibMapEEPhase_ = fileService_->make<TH2F>("calibMapEEPhase","time calib diff map EE (phase of Tmax)",360,1,361,170,-86,86);
  //calibMapEtaAvgEE_ = fileService_->make<TH2F>("calibMapEtaAvgEE","time calib diffs raw eta avg map EE",360,1,361,170,-86,86);
  //calibHistEtaAvgEE_ = fileService_->make<TH1F>("timingCalibsEtaAvgEE","EtaAvgTimingCalibDiffs EE [ns]",2000,-100,100);
  hitsPerCryMapEEM_ = fileService_->make<TH2F>("hitsPerCryMapEEM","Hits per cry EEM;ix;iy",100,1,101,100,1,101);
  hitsPerCryMapEEP_ = fileService_->make<TH2F>("hitsPerCryMapEEP","Hits per cry EEP;ix;iy",100,1,101,100,1,101);
  hitsPerCryHistEEM_ = fileService_->make<TH1F>("hitsPerCryHistEEM","Hits per cry EEM;hashedIndex",14648,0,14648);
  hitsPerCryHistEEP_ = fileService_->make<TH1F>("hitsPerCryHistEEP","Hits per cry EEP;hashedIndex",14648,0,14648);
  //eventsEEMHist_ = new TH1C("numEventsEEM","Number of events, EEM",100,0,100);
  //eventsEEPHist_ = new TH1C("numEventsEEP","Number of events, EEP",100,0,100);
  ampProfileEEM_ = new TProfile("ampProfileEEM","Amp profile EEM;hashedIndex",14648,0,14648);
  ampProfileEEP_ = new TProfile("ampProfileEEP","Amp profile EEP;hashedIndex",14648,0,14648);
  ampProfileMapEEP_ = new TProfile2D("ampProfileMapEEP","Amp profile EEP;ix;iy",100,1,101,100,1,101);
  ampProfileMapEEM_ = new TProfile2D("ampProfileMapEEM","Amp profile EEM;ix;iy",100,1,101,100,1,101);
  //eventsEEHist_ = fileService_->make<TH1F>("numEventsEE","Number of events, EE",100,0,100);
  //calibSigmaHist_ = fileService_->make<TH1F>("timingSpreadEE","Crystal timing spread [ns]",1000,-5,5);
  sigmaHistEE_ = fileService_->make<TH1F>("sigmaCalibDiffsEE"," Sigma of calib diff distributions EE [ns]",100,0,1);
  //chiSquaredEachEventHist_ = fileService_->make<TH1F>("chi2eachEvent","Chi2 of each event",500,0,500);
  //chiSquaredVsAmpEachEventHist_ = fileService_->make<TH2F>("chi2VsAmpEachEvent","Chi2 vs. amplitude of each event",500,0,500,750,0,750);
  //chiSquaredHighMap_ = fileService_->make<TH2F>("chi2HighMap","Channels with event #Chi^{2} > 100",360,1,361,170,-86,86);
  //chiSquaredTotalHist_ = fileService_->make<TH1F>("chi2Total","Total chi2 of all events in each crystal",500,0,500);
  //chiSquaredSingleOverTotalHist_ = fileService_->make<TH1F>("chi2SingleOverTotal","Chi2 of each event over total chi2",100,0,1);
  //ampEachEventHist_ = fileService_->make<TH1F>("energyEachEvent","Energy of all events [GeV]",1000,0,10);
  //numPointsErasedHist_ = fileService_->make<TH1F>("numPointsErased","Number of points erased per crystal",25,0,25);
  //myAmpProfile_ = (TProfile2D*)EBampProfile->Clone();
  //myAmpProfile_->Write();
  expectedStatPresHistEEM_ = fileService_->make<TH1F>("expectedStatPresEEM","Avg. expected statistical precision EEM [ns], all crys",200,0,2);
  expectedStatPresVsObservedMeanErrHistEEM_ = fileService_->make<TH2F>("expectedStatPresVsObsEEM","Expected stat. pres. vs. obs. error on mean each event EEM [ns]",200,0,2,200,0,2);
  expectedStatPresEachEventHistEEM_ = fileService_->make<TH1F>("expectedStatPresSingleEventEEM","Expected stat. pres. each event EEM [ns]",200,0,2);
  expectedStatPresHistEEP_ = fileService_->make<TH1F>("expectedStatPresEEP","Avg. expected statistical precision EEP [ns], all crys",200,0,2);
  expectedStatPresVsObservedMeanErrHistEEP_ = fileService_->make<TH2F>("expectedStatPresVsObsEEP","Expected stat. pres. vs. obs. error on mean each event [ns] EEP",200,0,2,200,0,2);
  expectedStatPresEachEventHistEEP_ = fileService_->make<TH1F>("expectedStatPresSingleEventEEP","Expected stat. pres. each event EEP [ns]",200,0,2);
  errorOnMeanVsNumEvtsHist_ = fileService_->make<TH2F>("errorOnMeanVsNumEvts","Error_on_mean vs. number of events",50,0,50,200,0,2);
  errorOnMeanVsNumEvtsHist_->Sumw2();
  calibHistEEM_ = fileService_->make<TH1F>("timingCalibDiffsEEM","timingCalibDiffs EEM [ns]",500,-25,25);
  calibHistEEP_ = fileService_->make<TH1F>("timingCalibDiffsEEP","timingCalibDiffs EEP [ns]",500,-25,25);
  calibErrorHistEEM_ = fileService_->make<TH1F>("calibDiffErrorEEM","timingCalibDiffError EEM [ns]",250,0,5);
  calibErrorHistEEP_ = fileService_->make<TH1F>("calibDiffErrorEEP","timingCalibDiffError EEP [ns]",250,0,5);
  calibHistEEM_->Sumw2();
  calibHistEEP_->Sumw2();
  calibErrorHistEEM_->Sumw2();
  calibErrorHistEEP_->Sumw2();
  calibMapEEM_ = fileService_->make<TH2F>("calibDiffMapEEM","time calib diff map EEM",100,1,101,100,1,101);
  calibMapEEP_ = fileService_->make<TH2F>("calibDiffMapEEP","time calib diff map EEP",100,1,101,100,1,101);
  calibMapEEM_->Sumw2();
  calibMapEEP_->Sumw2();
  sigmaMapEEM_ = fileService_->make<TH2F>("sigmaMapEEM","Sigma of time calib diff map EEM [ns];ix;iy",100,1.,101.,100,1,101);
  sigmaMapEEP_ = fileService_->make<TH2F>("sigmaMapEEP","Sigma of time calib diff map EEP [ns];ix;iy",100,1.,101.,100,1,101);
  calibErrorMapEEM_ = fileService_->make<TH2F>("calibErrorMapEEM","Error of time calib diff map EEM [ns];ix;iy",100,1.,101.,100,1,101);
  calibErrorMapEEP_ = fileService_->make<TH2F>("calibErrorMapEEP","Error of time calib diff map EEP [ns];ix;iy",100,1.,101.,100,1,101);

  TFileDirectory cryDirEEP = fileService_->mkdir("crystalTimingHistsEEP");
  for(int x=0; x < 100; ++x)
  {
    for(int y=0; y < 100; ++y)
    {
      if(!EEDetId::validDetId(x+1,y+1,1))
        continue;
      std::string histname = "EEP_cryTiming_ix";
      histname+=intToString(x+1);
      histname+="_iy";
      histname+=intToString(y+1);
      cryTimingHistsEEP_[x][y] = cryDirEEP.make<TH1F>(histname.c_str(),histname.c_str(),200,-10,10);
      cryTimingHistsEEP_[x][y]->Sumw2();
    }
  }
  TFileDirectory cryDirEEM = fileService_->mkdir("crystalTimingHistsEEM");
  for(int x=0; x < 100; ++x)
  {
    for(int y=0; y < 100; ++y)
    {
      if(!EEDetId::validDetId(x+1,y+1,-1))
        continue;
      std::string histname = "EEM_cryTiming_ix";
      histname+=intToString(x+1);
      histname+="_iy";
      histname+=intToString(y+1);
      cryTimingHistsEEM_[x][y] = cryDirEEM.make<TH1F>(histname.c_str(),histname.c_str(),200,-10,10);
      cryTimingHistsEEM_[x][y]->Sumw2();
    }
  }
}

// ****************************************************************
// Helper functions
std::string EcalCreateTimeCalibrations::intToString(int num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}

//
std::vector<std::string> EcalCreateTimeCalibrations::split(std::string msg, std::string separator)
{
  boost::char_separator<char> sep(separator.c_str());
  boost::tokenizer<boost::char_separator<char> > tok(msg, sep );
  std::vector<std::string> token ;
  for ( boost::tokenizer<boost::char_separator<char> >::const_iterator i = tok.begin(); i != tok.end(); ++i ) {
    token.push_back(std::string(*i)) ;
  }
  return token ;
}

//
void EcalCreateTimeCalibrations::genIncludeExcludeVectors(std::string optionString,
    std::vector<std::vector<double> >& includeVector,
    std::vector<std::vector<double> >& excludeVector)
{
  std::vector<std::string> rangeStringVector;
  std::vector<double> rangeIntVector;

  if(optionString != "-1"){
    std::vector<std::string> stringVector = split(optionString,",") ;

    for (uint i=0 ; i<stringVector.size() ; i++) {
      bool exclude = false;

      if(stringVector[i].at(0)=='x'){
        exclude = true;
        stringVector[i].erase(0,1);
      }
      rangeStringVector = split(stringVector[i],"-") ;

      rangeIntVector.clear();
      for(uint j=0; j<rangeStringVector.size();j++) {
        rangeIntVector.push_back(atof(rangeStringVector[j].c_str()));
      }
      if(exclude) excludeVector.push_back(rangeIntVector);
      else includeVector.push_back(rangeIntVector);

    }
  }
}


//
bool EcalCreateTimeCalibrations::includeEvent(double eventParameter,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
  bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for(uint i=0; i!=includeVector.size();++i){
    if(includeVector[i].size()==1 && eventParameter==includeVector[i][0])
      keepEvent=true;
    else if(includeVector[i].size()==2 && (eventParameter>=includeVector[i][0] && eventParameter<=includeVector[i][1]))
      keepEvent=true;
  }
  if(!keepEvent) // if it's not in our include list, skip it
    return false;

  keepEvent = true;
  for(uint i=0; i!=excludeVector.size();++i){
    if(excludeVector[i].size()==1 && eventParameter==excludeVector[i][0])
      keepEvent=false;
    else if(excludeVector[i].size()==2 && (eventParameter>=excludeVector[i][0] && eventParameter<=excludeVector[i][1]))
      keepEvent=false;
  }

  return keepEvent; // if someone includes and excludes, exclusion will overrule

}


//
bool EcalCreateTimeCalibrations::includeEvent(int* triggers,
    int numTriggers,
    std::vector<std::vector<double> > includeVector,
    std::vector<std::vector<double> > excludeVector)
{
  bool keepEvent = false;
  if(includeVector.size()==0) keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=includeVector.size();++i){
      if(includeVector[i].size()==1 && triggers[ti]==includeVector[i][0]) keepEvent=true;
      else if(includeVector[i].size()==2 && (triggers[ti]>=includeVector[i][0] && triggers[ti]<=includeVector[i][1])) keepEvent=true;
    }
  }
  if(!keepEvent)
    return false;

  keepEvent = true;
  for (int ti = 0; ti < numTriggers; ++ti) {
    for(uint i=0; i!=excludeVector.size();++i){
      if(excludeVector[i].size()==1 && triggers[ti]==excludeVector[i][0]) keepEvent=false;
      else if(excludeVector[i].size()==2 && (triggers[ti]>=excludeVector[i][0] && triggers[ti]<=excludeVector[i][1])) keepEvent=false;
    }
  }

  return keepEvent;
}
