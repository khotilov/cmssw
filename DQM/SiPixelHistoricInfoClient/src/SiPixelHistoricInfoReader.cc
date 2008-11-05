#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelPerformanceSummary.h"
#include "CondFormats/DataRecord/interface/SiPixelPerformanceSummaryRcd.h"

#include "DQM/SiPixelHistoricInfoClient/interface/SiPixelHistoricInfoReader.h"

#include <math.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <sys/time.h>

#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"


using namespace cms;
using namespace std;


SiPixelHistoricInfoReader::SiPixelHistoricInfoReader(const edm::ParameterSet& pSet) {
  parameterSet_ = pSet;  

  variables_ = parameterSet_.getUntrackedParameter<vstring>("variables");
  for (int i=0; i<15; i++) variable_[i] = false; 
  for (vector<string>::const_iterator variable = variables_.begin(); variable!=variables_.end(); ++variable) {
    if (variable->compare("errorType")==0)   variable_[ 0] = true; 
    if (variable->compare("ndigis")==0)      variable_[ 1] = true; 
    if (variable->compare("adc")==0)         variable_[ 2] = true; 
    if (variable->compare("nclusters")==0)   variable_[ 3] = true; 
    if (variable->compare("charge")==0)      variable_[ 4] = true; 
    if (variable->compare("size")==0)        variable_[ 5] = true; 
    if (variable->compare("sizeX")==0)       variable_[ 6] = true; 
    if (variable->compare("sizeY")==0)       variable_[ 7] = true; 
    if (variable->compare("nRecHits")==0)    variable_[ 8] = true; 
    if (variable->compare("residualX")==0)   variable_[ 9] = true; 
    if (variable->compare("residualY")==0)   variable_[10] = true; 
    if (variable->compare("nPixHitsTrk")==0) variable_[11] = true; 
    if (variable->compare("nNoisPixels")==0) variable_[12] = true; 
    if (variable->compare("nDeadPixels")==0) variable_[13] = true; 
  }
  normEvents_ = parameterSet_.getUntrackedParameter<bool>("normEvents",false);  
  printDebug_ = parameterSet_.getUntrackedParameter<bool>("printDebug",false); 
   makePlots_ = parameterSet_.getUntrackedParameter<bool>("makePlots", true); 
   typePlots_ = parameterSet_.getUntrackedParameter<string>("typePlots", "gif"); 
  outputDir_  = parameterSet_.getUntrackedParameter<string>("outputDir", ".");
  outputFile_ = parameterSet_.getUntrackedParameter<string>("outputFile","pixelhistory.root");
  
  firstBeginRun_ = true; 
} 


SiPixelHistoricInfoReader::~SiPixelHistoricInfoReader() {}


void SiPixelHistoricInfoReader::beginJob(const edm::EventSetup& iSetup) {
  string outputDirFile = outputDir_ + "/" + outputFile_; 
  outputDirFile_ = new TFile(outputDirFile.data(), "RECREATE");
} 


string SiPixelHistoricInfoReader::getMEregionString(uint32_t detID) const {
  uint32_t localMEdetID = detID; 
  TString regionStr; 
       if (localMEdetID>100000000) { regionStr = "det"; regionStr += localMEdetID; }
  else if (localMEdetID<40)        { regionStr = "FED"; regionStr += localMEdetID; }
  else if (localMEdetID>99 && localMEdetID<120) { 
    localMEdetID -= 100; 
    if (localMEdetID<12) { 
      regionStr = "Barrel/"; 
    	   if (localMEdetID<3) { regionStr += "Shell_mI/";	              }
      else if (localMEdetID<6) { regionStr += "Shell_mO/"; localMEdetID -= 3; }
      else if (localMEdetID<9) { regionStr += "Shell_pI/"; localMEdetID -= 6; } 
      else		       { regionStr += "Shell_pO/"; localMEdetID -= 9; }
      regionStr += "Layer_"; regionStr += (localMEdetID+1); 
    }
    else { 
      regionStr = "Endcap/"; localMEdetID -= 12; 
    	   if (localMEdetID<2) { regionStr += "HalfCylinder_mI/";		     }
      else if (localMEdetID<4) { regionStr += "HalfCylinder_mO/"; localMEdetID -= 2; }
      else if (localMEdetID<6) { regionStr += "HalfCylinder_pI/"; localMEdetID -= 4; } 
      else		       { regionStr += "HalfCylinder_pO/"; localMEdetID -= 6; }
      regionStr += "Disk_"; regionStr += (localMEdetID+1); 
    }
  } 
  else {
    localMEdetID -= 1000; 
    if (localMEdetID<192) { 
      regionStr = "Barrel/"; 
    	   if (localMEdetID<48)  { regionStr += "Shell_mI/";			  }
      else if (localMEdetID<96)  { regionStr += "Shell_mO/"; localMEdetID -= 48;  }
      else if (localMEdetID<144) { regionStr += "Shell_pI/"; localMEdetID -= 96;  } 
      else			 { regionStr += "Shell_pO/"; localMEdetID -= 144; }
           if (localMEdetID<10) { regionStr += "Layer_1/";		       }
      else if (localMEdetID<26) { regionStr += "Layer_2/"; localMEdetID -= 10; }
      else			{ regionStr += "Layer_3/"; localMEdetID -= 26; }
      regionStr += "Ladder_"; regionStr += (localMEdetID+1); 
    }
    else { 
      regionStr = "Endcap/"; localMEdetID -= 192; 
    	   if (localMEdetID<24) { regionStr += "HalfCylinder_mI/";		       }
      else if (localMEdetID<48) { regionStr += "HalfCylinder_mO/"; localMEdetID -= 24; }
      else if (localMEdetID<72) { regionStr += "HalfCylinder_pI/"; localMEdetID -= 48; } 
      else			{ regionStr += "HalfCylinder_pO/"; localMEdetID -= 72; }
      if (localMEdetID<12) { regionStr += "Disk_1/"; }
      else		   { regionStr += "Disk_2/"; localMEdetID -= 12; } 
      regionStr += "Blade_"; regionStr += (localMEdetID+1); 
    }
  }
  return regionStr.Data(); 
}


void SiPixelHistoricInfoReader::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {
  edm::ESHandle<SiPixelPerformanceSummary> pSummary;
  iSetup.get<SiPixelPerformanceSummaryRcd>().get(pSummary); 

  if (firstBeginRun_) {
    firstBeginRun_ = false;

    allDetIds.clear(); 
    allDetIds = pSummary->getAllDetIds(); // allDetIds.push_back(369345800);
      
    AllDetHistograms = new TObjArray();

    for (vector<uint32_t>::const_iterator iDet = allDetIds.begin(); iDet!=allDetIds.end(); ++iDet) {  
      string detRegion = getMEregionString(*iDet); 
      
      if (variable_[0] && *iDet<40) { 
        for (int pBin=0; pBin<15; pBin++) {
  	  hisID = "errorType"; hisID += (pBin+25); hisID += "_"; hisID += *iDet; 	     
          title = "errorType"; title += (pBin+25); title += " "; title += detRegion; 	        
  	  AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1));
  	  ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin); 
        } 
      }
      if (variable_[1] && *iDet>99) {
        hisID = "nDigis_"; hisID += *iDet; 	        
        title = "nDigis "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "nDigis_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in nDigis "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[2] && *iDet>99) {
        hisID = "adc_"; hisID += *iDet; 		        
        title = "adc "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "adc_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in adc "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[3] && *iDet>99) {
        hisID = "nClusters_"; hisID += *iDet;         
        title = "nClusters "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "nClusters_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in nClusters "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[4] && *iDet>99) {
        hisID = "charge_"; hisID += *iDet; 	        
        title = "charge "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "charge_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in charge "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[5] && *iDet>99) {
        hisID = "clusterSize_"; hisID += *iDet; 	       
        title = "clusterSize "; title += detRegion;	       
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "clusterSize_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in clusterSize "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[6] && *iDet>99) {
        hisID = "clusterSizeX_"; hisID += *iDet; 	        
        title = "clusterSizeX "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "clusterSizeX_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in clusterSizeX "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[7] && *iDet>99) {
        hisID = "clusterSizeY_"; hisID += *iDet; 	        
        title = "clusterSizeY "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "clusterSizeY_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in clusterSizeY "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[8] && *iDet>99) {
        hisID = "nRecHits_"; hisID += *iDet;          
        title = "nRecHits "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "nRecHits_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in nRecHits "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[9] && *iDet>99) {
        hisID = "residualX_"; hisID += *iDet; 	        
        title = "residualX "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "residualX_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in residualX "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      }
      if (variable_[10] && *iDet>99) {
        hisID = "residualY_"; hisID += *iDet; 	        
        title = "residualY "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  

        hisID = "residualY_mFr_"; hisID += *iDet; 	        
        title = "Fraction of Empty Modules in residualY "; title += detRegion; 	        
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      } 
      if (variable_[11] && *iDet>99) {
        hisID = "nPixHitsTrk_"; hisID += *iDet;	       
        title = "nPixHitsTrk "; title += detRegion;	       
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      } 
      if (variable_[12] && *iDet>99) {
        hisID = "nNoisPixels_"; hisID += *iDet;	       
        title = "nNoisPixels "; title += detRegion;	       
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      } 
      if (variable_[13] && *iDet>99) {
        hisID = "nDeadPixels_"; hisID += *iDet;	       
        title = "nDeadPixels "; title += detRegion;	       
        AllDetHistograms->Add(new TH1F(hisID, title, 1, 0, 1)); 	        
        ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBit(TH1::kCanRebin);  
      } 
    }
  }
  if (pSummary->getRunNumber()==run.run()) { // pSummary's run changes only when the table is newly retrieved 
    TString sRun; sRun += pSummary->getRunNumber();     
    float nEvents = pSummary->getNumberOfEvents(); 
    float SF = 1.0; if (normEvents_) SF = 100000.0/nEvents; 
    
    if (printDebug_) cout << "run "<< sRun.Data() <<" with "<< nEvents <<" events" << endl; 
    
    for (vector<uint32_t>::const_iterator iDet = allDetIds.begin(); iDet!=allDetIds.end(); ++iDet) {
      vector<float> performances = pSummary->getDetSummary(*iDet); 
      
      if (*iDet<40) {
      	if (variable_[0]) {
      	  for (int pBin=0; pBin<15; pBin++) {
      	    hisID = "errorType"; hisID += (pBin+25); hisID += "_"; hisID += *iDet;	
      	    ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[pBin]/nEvents);
	    int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
            ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, sqrt(performances[pBin])/nEvents);        
      	  } 
      	} 
      } 
      if (*iDet>99) {
      	if (variable_[1]) {
      	  hisID = "nDigis_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[15]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[16]*SF);

      	  hisID = "nDigis_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[17]);
      	} 
      	if (variable_[2]) {      
      	  hisID = "adc_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[18]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[19]*SF);

      	  hisID = "adc_mFr_"; hisID += *iDet; 	 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[20]);
      	} 
      	if (variable_[3]) { 	  
      	  hisID = "nClusters_"; hisID += *iDet;        
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[21]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[22]*SF);    

      	  hisID = "nClusters_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[23]);
      	} 
      	if (variable_[4]) {      
      	  hisID = "charge_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[24]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[25]*SF);    

      	  hisID = "charge_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[26]);
      	} 
      	if (variable_[5]) {      
      	  hisID = "clusterSize_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[27]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[28]*SF);

      	  hisID = "clusterSize_mFr_"; hisID += *iDet; 	 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[29]);
      	} 
      	if (variable_[6]) {      
      	  hisID = "clusterSizeX_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[30]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[31]*SF);

      	  hisID = "clusterSizeX_mFr_"; hisID += *iDet;	 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[32]);
      	} 
      	if (variable_[7]) {      
      	  hisID = "clusterSizeY_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[33]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[34]*SF);    

      	  hisID = "clusterSizeY_mFr_"; hisID += *iDet;	 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[35]);
      	} 
      	if (variable_[8]) {      
      	  hisID = "nRecHits_"; hisID += *iDet;         
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[36]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[37]*SF);    

      	  hisID = "nRecHits_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[38]);
      	} 
      	if (variable_[9]) {      
      	  hisID = "residualX_"; hisID += *iDet;       
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[39]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[40]*SF);    

      	  hisID = "residualX_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[41]);
      	} 
      	if (variable_[10]) {      
      	  hisID = "residualY_"; hisID += *iDet;       
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[42]*SF);
      	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, performances[43]*SF);    

      	  hisID = "residualY_mFr_"; hisID += *iDet;		 
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[44]);
      	} 
      	if (variable_[11]) {      
      	  hisID = "nPixHitsTrk_"; hisID += *iDet;	
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[45]*SF);
      	} 
      	if (variable_[12]) {      
      	  hisID = "nNoisPixels_"; hisID += *iDet;	
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[46]*SF);
	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, sqrt(performances[46])/sqrt(nEvents));	    
      	} 
      	if (variable_[13]) {      
      	  hisID = "nDeadPixels_"; hisID += *iDet;	
      	  ((TH1F*)AllDetHistograms->FindObject(hisID))->Fill(sRun, performances[47]*SF);
	  int iBin = ((TH1F*)AllDetHistograms->FindObject(hisID))->GetXaxis()->FindBin(sRun);
          ((TH1F*)AllDetHistograms->FindObject(hisID))->SetBinError(iBin, sqrt(performances[47])/sqrt(nEvents));	    
      	} 
      }
    }    
  }
}


void SiPixelHistoricInfoReader::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {}


void SiPixelHistoricInfoReader::endRun(const edm::Run& run, const edm::EventSetup& iSetup) {}


void SiPixelHistoricInfoReader::endJob() {
  for (vector<uint32_t>::const_iterator iDet=allDetIds.begin(); iDet!=allDetIds.end(); ++iDet) {
    if (*iDet<40) {
      if (variable_[0]) {
        for (int pBin=0; pBin<15; pBin++) {				  	   		       
          hisID = "errorType"; hisID += (pBin+25); hisID += "_"; hisID += *iDet;     
          ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");	   		       
        }								  	   		       
      } 
    } 
    if (*iDet>99) {
      if (variable_[1]) {
        hisID = "nDigis_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "nDigis_mFr_"; hisID += *iDet;			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[2]) {
        hisID = "adc_"; hisID += *iDet; 				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X"); 

        hisID = "adc_mFr_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[3]) {
        hisID = "nClusters_"; hisID += *iDet;				 
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "nClusters_mFr_"; hisID += *iDet;			    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[4]) {
        hisID = "charge_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "charge_mFr_"; hisID += *iDet;			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[5]) {
        hisID = "clusterSize_"; hisID += *iDet; 			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "clusterSize_mFr_"; hisID += *iDet;			    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[6]) {
        hisID = "clusterSizeX_"; hisID += *iDet;			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "clusterSizeX_mFr_"; hisID += *iDet;  		    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[7]) {
        hisID = "clusterSizeY_"; hisID += *iDet;			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "clusterSizeY_mFr_"; hisID += *iDet;  		    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[8]) {
        hisID = "nRecHits_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   														      

        hisID = "nRecHits_mFr_"; hisID += *iDet;			    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[9]) {
        hisID = "residualX_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   														      

        hisID = "residualX_mFr_"; hisID += *iDet;			    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      } 
      if (variable_[10]) {
        hisID = "residualY_"; hisID += *iDet;				       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       

        hisID = "residualY_mFr_"; hisID += *iDet;			    	       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      }
      if (variable_[11]) {
        hisID = "nPixHitsTrk_"; hisID += *iDet; 			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      }
      if (variable_[12]) {
        hisID = "nNoisPixels_"; hisID += *iDet; 			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      }
      if (variable_[13]) {
        hisID = "nDeadPixels_"; hisID += *iDet; 			       
        ((TH1F*)AllDetHistograms->FindObject(hisID))->LabelsDeflate("X");   		       
      }
    } 
  } 
  if (makePlots_) plot(); 
  
  outputDirFile_->Write();
  outputDirFile_->Close();
}


void SiPixelHistoricInfoReader::plot() { 
  ostringstream spa; for (unsigned int i=0; i<typePlots_.size(); i++) spa << " "; 
  cout << "************************************ "<< spa.str() <<" *************" << endl; 
  cout << "SiPixelHistoricInfoReader::plot() in "<< typePlots_<<" with patience" << endl;
  cout << "************************************ "<< spa.str() <<" *************" << endl;

  TStyle* sty1 = new TStyle("sty1","My Style.1 for Pixel History Plots"); 
	  sty1->SetCanvasDefW(1080); 
	  sty1->SetCanvasDefH(250); 
	  sty1->SetCanvasBorderMode(0); 
	  sty1->SetPadBorderMode(0);
	  sty1->SetPadLeftMargin(0.05);
  	  sty1->SetPadRightMargin(0.03);
  	  sty1->SetPadTopMargin(0.16);
  	  sty1->SetPadBottomMargin(0.15);
	  sty1->SetTitleBorderSize(0); 
  	  sty1->SetTitleFontSize(0.07); 
	  sty1->SetOptStat(0); 
          sty1->cd(); 

  TCanvas* c1 = new TCanvas("c1","c1",1080,250); 
           c1->UseCurrentStyle(); 
  
  for (vector<uint32_t>::const_iterator iDet=allDetIds.begin(); iDet!=allDetIds.end(); ++iDet) {
    if (*iDet<40) {
      if (variable_[0]) {
        for (int pBin=0; pBin<15; pBin++) {				  	   		       
          hisID = "errorType"; hisID += (pBin+25); hisID += "_"; hisID += *iDet;     
          TH1F* errFED = (TH1F*)AllDetHistograms->FindObject(hisID);	   		       
          if (errFED->Integral()>0.0) {
	    errFED->SetMinimum(0.0); 
      	    errFED->GetXaxis()->SetLabelSize(0.08); 
      	    errFED->GetYaxis()->SetLabelSize(0.06); 
      	    errFED->SetMarkerStyle(8); 
      	    errFED->SetMarkerSize(0.2); 
	    errFED->Draw(); 
	    title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	    c1->SaveAs(title); 
          }
	}								  	   		       
      } 
    } 
    if (*iDet>99) {
      if (variable_[1]) {
        hisID = "nDigis_"; hisID += *iDet;				       
        TH1F* nDigis = (TH1F*)AllDetHistograms->FindObject(hisID);		      
	      nDigis->SetMinimum(0.0); 
	      nDigis->SetMarkerStyle(8); 
	      nDigis->SetMarkerSize(0.2); 
	      nDigis->GetXaxis()->SetLabelSize(0.08); 
	      nDigis->GetYaxis()->SetLabelSize(0.06); 
	      nDigis->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[2]) {
        hisID = "adc_"; hisID += *iDet; 				       
	TH1F* adc = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      adc->SetMinimum(0.0); 
	      adc->SetMarkerStyle(8); 
	      adc->SetMarkerSize(0.2); 
	      adc->GetXaxis()->SetLabelSize(0.08); 
	      adc->GetYaxis()->SetLabelSize(0.06); 
	      adc->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[3]) {
        hisID = "nClusters_"; hisID += *iDet;				 
	TH1F* nClusters = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      nClusters->SetMinimum(0.0); 
	      nClusters->SetMarkerStyle(8); 
	      nClusters->SetMarkerSize(0.2); 
	      nClusters->GetXaxis()->SetLabelSize(0.08); 
	      nClusters->GetYaxis()->SetLabelSize(0.06); 
	      nClusters->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[4]) {
        hisID = "charge_"; hisID += *iDet;				       
	TH1F* charge = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      charge->SetMinimum(0.0); 
	      charge->SetMarkerStyle(8); 
	      charge->SetMarkerSize(0.2); 
	      charge->GetXaxis()->SetLabelSize(0.08); 
	      charge->GetYaxis()->SetLabelSize(0.06); 
	      charge->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[5]) {
        hisID = "clusterSize_"; hisID += *iDet; 			       
	TH1F* clusterSize = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      clusterSize->SetMinimum(0.0); 
	      clusterSize->SetMarkerStyle(8); 
	      clusterSize->SetMarkerSize(0.2); 
	      clusterSize->GetXaxis()->SetLabelSize(0.08); 
	      clusterSize->GetYaxis()->SetLabelSize(0.06); 
	      clusterSize->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[6]) {
        hisID = "clusterSizeX_"; hisID += *iDet;			       
	TH1F* sizeX = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      sizeX->SetMinimum(0.0); 
	      sizeX->SetMarkerStyle(8); 
	      sizeX->SetMarkerSize(0.2); 
	      sizeX->GetXaxis()->SetLabelSize(0.08); 
	      sizeX->GetYaxis()->SetLabelSize(0.06); 
	      sizeX->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[7]) {
        hisID = "clusterSizeY_"; hisID += *iDet;			       
	TH1F* sizeY = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      sizeY->SetMinimum(0.0); 
	      sizeY->SetMarkerStyle(8); 
	      sizeY->SetMarkerSize(0.2); 
	      sizeY->GetXaxis()->SetLabelSize(0.08); 
	      sizeY->GetYaxis()->SetLabelSize(0.06); 
	      sizeY->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      }
      if (variable_[8]) {
        hisID = "nRecHits_"; hisID += *iDet;				       
	TH1F* nRecHits = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      nRecHits->SetMinimum(0.0); 
	      nRecHits->SetMarkerStyle(8); 
	      nRecHits->SetMarkerSize(0.2); 
	      nRecHits->GetXaxis()->SetLabelSize(0.08); 
	      nRecHits->GetYaxis()->SetLabelSize(0.06); 
	      nRecHits->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[9]) {
        hisID = "residualX_"; hisID += *iDet;				       
	TH1F* residualX = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      residualX->SetMinimum(0.0); 
	      residualX->SetMarkerStyle(8); 
	      residualX->SetMarkerSize(0.2); 
	      residualX->GetXaxis()->SetLabelSize(0.08); 
	      residualX->GetYaxis()->SetLabelSize(0.06); 
	      residualX->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      } 
      if (variable_[10]) {
        hisID = "residualY_"; hisID += *iDet;				       
	TH1F* residualY = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      residualY->SetMinimum(0.0); 
	      residualY->SetMarkerStyle(8); 
	      residualY->SetMarkerSize(0.2); 
	      residualY->GetXaxis()->SetLabelSize(0.08); 
	      residualY->GetYaxis()->SetLabelSize(0.06); 
	      residualY->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      }
      if (variable_[11]) {
        hisID = "nPixHitsTrk_"; hisID += *iDet; 			       
	TH1F* nPixHitsTrk = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      nPixHitsTrk->SetMinimum(0.0); 
	      nPixHitsTrk->SetMarkerStyle(8); 
	      nPixHitsTrk->SetMarkerSize(0.2); 
	      nPixHitsTrk->GetXaxis()->SetLabelSize(0.08); 
	      nPixHitsTrk->GetYaxis()->SetLabelSize(0.06); 
	      nPixHitsTrk->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      }
      if (variable_[12]) {
        hisID = "nNoisPixels_"; hisID += *iDet; 			       
	TH1F* nNoisPixels = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      nNoisPixels->SetMinimum(0.0); 
	      nNoisPixels->SetMarkerStyle(8); 
	      nNoisPixels->SetMarkerSize(0.2); 
	      nNoisPixels->GetXaxis()->SetLabelSize(0.08); 
	      nNoisPixels->GetYaxis()->SetLabelSize(0.06); 
	      nNoisPixels->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      }
      if (variable_[13]) {
        hisID = "nDeadPixels_"; hisID += *iDet; 			       
	TH1F* nDeadPixels = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      nDeadPixels->SetMinimum(0.0); 
	      nDeadPixels->SetMarkerStyle(8); 
	      nDeadPixels->SetMarkerSize(0.2); 
	      nDeadPixels->GetXaxis()->SetLabelSize(0.08); 
	      nDeadPixels->GetYaxis()->SetLabelSize(0.06); 
	      nDeadPixels->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c1->SaveAs(title); 
      }
    } 
  } 
  TStyle* sty2 = new TStyle("mSty2","My Style.2 for Pixel History Plots"); 
	  sty2->SetCanvasDefW(1080); 
	  sty2->SetCanvasDefH(135); 
	  sty2->SetCanvasBorderMode(0); 
	  sty2->SetPadBorderMode(0);
	  sty2->SetPadLeftMargin(0.05);
  	  sty2->SetPadRightMargin(0.03);
  	  sty2->SetPadTopMargin(0.18);
  	  sty2->SetPadBottomMargin(0.31);
	  sty2->SetTitleBorderSize(0); 
  	  sty2->SetTitleFontSize(0.144); 
  	  sty2->SetTitleX(0.409); 
	  sty2->SetOptStat(0); 
          sty2->cd(); 
 
  TCanvas* c2 = new TCanvas("c2","c2",1080,135); 
           c2->UseCurrentStyle(); 

  for (vector<uint32_t>::const_iterator iDet=allDetIds.begin(); iDet!=allDetIds.end(); ++iDet) {
    if (*iDet>99) {
      if (variable_[1]) {
        hisID = "nDigis_mFr_"; hisID += *iDet;			       
        (TH1F*)AllDetHistograms->FindObject(hisID);   		       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[2]) {
        hisID = "adc_mFr_"; hisID += *iDet;				       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[3]) {
        hisID = "nClusters_mFr_"; hisID += *iDet;			    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[4]) {
        hisID = "charge_mFr_"; hisID += *iDet;			       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[5]) {
        hisID = "clusterSize_mFr_"; hisID += *iDet;			    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[6]) {
        hisID = "clusterSizeX_mFr_"; hisID += *iDet;  		    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[7]) {
        hisID = "clusterSizeY_mFr_"; hisID += *iDet;  		    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[8]) {
        hisID = "nRecHits_mFr_"; hisID += *iDet;			    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[9]) {
        hisID = "residualX_mFr_"; hisID += *iDet;			    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      } 
      if (variable_[10]) {
        hisID = "residualY_mFr_"; hisID += *iDet;			    	       
	TH1F* mFr = (TH1F*)AllDetHistograms->FindObject(hisID); 
	      mFr->SetMinimum(0.0); 
	      mFr->GetXaxis()->SetLabelOffset(0.009); 
	      mFr->GetXaxis()->SetLabelSize(0.165); 
	      mFr->GetYaxis()->SetLabelSize(0.120); 
	      mFr->Draw(); 
	title = outputDir_; title += "/"; title += hisID; title += "."; title += typePlots_; 
	c2->SaveAs(title); 
      }
    } 
  } 
}
