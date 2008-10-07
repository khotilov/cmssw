

/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/10/07 09:56:58 $
 *  $Revision: 1.17 $
 *  \author G. Mila - INFN Torino
 */


#include <DQM/DTMonitorClient/src/DTSegmentAnalysisTest.h>

// Framework
#include <FWCore/Framework/interface/Event.h>
#include "DataFormats/Common/interface/Handle.h" 
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>


// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>


using namespace edm;
using namespace std;


DTSegmentAnalysisTest::DTSegmentAnalysisTest(const edm::ParameterSet& ps){

  edm::LogVerbatim ("segment") << "[DTSegmentAnalysisTest]: Constructor";
  parameters = ps;

  dbe = edm::Service<DQMStore>().operator->();

  // get the cfi parameters
  detailedAnalysis = parameters.getUntrackedParameter<bool>("detailedAnalysis","false");
}


DTSegmentAnalysisTest::~DTSegmentAnalysisTest(){

  edm::LogVerbatim ("segment") << "DTSegmentAnalysisTest: analyzed " << nevents << " events";
}


void DTSegmentAnalysisTest::beginJob(const edm::EventSetup& context){

  edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: BeginJob"; 

  nevents = 0;
  // Get the geometry
  context.get<MuonGeometryRecord>().get(muonGeom);

  // book the histos
  bookHistos();  

}


void DTSegmentAnalysisTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: Begin of LS transition";

}


void DTSegmentAnalysisTest::analyze(const edm::Event& e, const edm::EventSetup& context){
 
  nevents++;
  edm::LogVerbatim ("segment") << "[DTSegmentAnalysisTest]: "<<nevents<<" events";

}


void DTSegmentAnalysisTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {
 
  edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: End of LS transition, performing the DQM client operation";

  // counts number of lumiSegs 
  nLumiSegs = lumiSeg.id().luminosityBlock();
 
  edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: "<<nLumiSegs<<" updates";

  summaryHistos[3]->Reset();
  vector<DTChamber*>::const_iterator ch_it = muonGeom->chambers().begin();
  vector<DTChamber*>::const_iterator ch_end = muonGeom->chambers().end();
 
  for (; ch_it != ch_end; ++ch_it) {
    DTChamberId chID = (*ch_it)->id();
    
    MonitorElement * segm_histo = dbe->get(getMEName(chID, "h4DSegmNHits"));
    MonitorElement * summary_histo = dbe->get(getMEName(chID, "numberOfSegments"));
   
    if (segm_histo && summary_histo){
      edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: I've got the recHits histo and the summary!!";
      
      TH1F * segmHit_histo_root = segm_histo->getTH1F();
      TH2F * segm_histo_root = summary_histo->getTH2F();
      TH2F * summary_histo_root = summaryHistos[3]->getTH2F();
      
      int sector = chID.sector();
      if(sector == 13) sector=4;
      if(sector == 14) sector=10;
      
      
      if((chID.station()!=4 && segmHit_histo_root->GetMaximumBin() != 12)||
	 (chID.station()==4 &&  segmHit_histo_root->GetMaximumBin() != 8)){
	summaryHistos[chID.wheel()]->setBinContent(sector, chID.station(),1);
	if(summary_histo_root->GetBinContent(sector, chID.wheel()+3)<1)
	  summaryHistos[3]->setBinContent(sector, chID.wheel()+3,1);  
      }
      else
	summaryHistos[chID.wheel()]->setBinContent(sector, chID.station(),0);
    
      if(detailedAnalysis){
	if(chID.station()!=4)
	  segmRecHitHistos[make_pair(chID.wheel(),chID.sector())]->Fill(chID.station(),abs(12-segmHit_histo_root->GetMaximumBin()));
	else
	   segmRecHitHistos[make_pair(chID.wheel(),chID.sector())]->Fill(chID.station(),abs(8-segmHit_histo_root->GetMaximumBin()));
      }

      TH2F * summary2_histo_root = summaryHistos[3]->getTH2F();
      
      if(segm_histo_root->GetBinContent(sector,chID.station())==0){
	summaryHistos[chID.wheel()]->setBinContent(sector, chID.station(),2);
	if(summary2_histo_root->GetBinContent(sector, chID.wheel()+3)<2)
	  summaryHistos[3]->setBinContent(sector, chID.wheel()+3,2);
      }
      
    }

    if(detailedAnalysis){ // switch on detailed analysis
   
      //test on chi2 segment quality
      MonitorElement * chi2_histo = dbe->get(getMEName(chID, "h4DChi2"));
      if (chi2_histo) {
	edm::LogVerbatim ("segment") <<"[DTSegmentAnalysisTest]: I've got the histo of the segment chi2!!";
	TH1F * chi2_histo_root = chi2_histo->getTH1F();
	double threshold = parameters.getUntrackedParameter<double>("chi2Threshold", 5);
	double maximum = chi2_histo_root->GetXaxis()->GetXmax();
	double minimum = chi2_histo_root->GetXaxis()->GetXmin();
	int nbins = chi2_histo_root->GetXaxis()->GetNbins();
	int thresholdBin = int(threshold/((maximum-minimum)/nbins));
	
	double badSegments=0;
	for(int bin=thresholdBin; bin<=nbins; bin++){
	  badSegments+=chi2_histo_root->GetBinContent(bin);
	}
      
	if(chi2_histo_root->GetEntries()!=0){
	  double badSegmentsPercentual= badSegments/double(chi2_histo_root->GetEntries());
	  chi2Histos[make_pair(chID.wheel(),chID.sector())]->Fill(chID.station(),badSegmentsPercentual);
	}
      }      
    } // end of switch for detailed analysis
    
  } //loop over all the chambers
  

  if(detailedAnalysis){
    
    string chi2CriterionName = parameters.getUntrackedParameter<string>("chi2TestName","chi2InRange");
    for(map<pair<int, int>, MonitorElement*> ::const_iterator histo = chi2Histos.begin();
	histo != chi2Histos.end();
	histo++) {

      const QReport * theChi2QReport = (*histo).second->getQReport(chi2CriterionName);
      if(theChi2QReport) {
	vector<dqm::me_util::Channel> badChannels = theChi2QReport->getBadChannels();
	for (vector<dqm::me_util::Channel>::iterator channel = badChannels.begin(); 
	     channel != badChannels.end(); channel++) {
	  edm::LogError ("dtSegmentsTest") << "Wheel: "<<(*histo).first.first<< " Sector: "<<(*histo).first.second<< " Bad stations: "<<(*channel).getBin()<<"  Contents : "<<(*channel).getContents();
	}
      }
    }
    
    string segmRecHitCriterionName = parameters.getUntrackedParameter<string>("segmRecHitTestName","segmRecHitInRange");
    for(map<pair<int, int>, MonitorElement*> ::const_iterator histo = segmRecHitHistos.begin();
	histo != segmRecHitHistos.end();
	histo++) {

      const QReport * theSegmRecHitQReport = (*histo).second->getQReport(segmRecHitCriterionName);
      if(theSegmRecHitQReport) {
	vector<dqm::me_util::Channel> badChannels = theSegmRecHitQReport->getBadChannels();
	for (vector<dqm::me_util::Channel>::iterator channel = badChannels.begin(); 
	     channel != badChannels.end(); channel++) {
	  edm::LogError ("dtSegmentsTest") << "Wheel: "<<(*histo).first.first<< " Sector: "<<(*histo).first.second<< " Bad stations on recHit number: "<<(*channel).getBin()<<"  Contents : "<<(*channel).getContents();
	}
      }
    }

  } // end of detailedAnalysis

}


string DTSegmentAnalysisTest::getMEName(const DTChamberId & chID, string histoTag) {
  
  stringstream wheel; wheel << chID.wheel();	
  stringstream station; station << chID.station();	
  stringstream sector; sector << chID.sector();	
  
  string folderRoot = parameters.getUntrackedParameter<string>("folderRoot", "Collector/FU0/");
  string folderName = 
    folderRoot + "DT/02-Segments/Wheel" +  wheel.str() +
    "/Station" + station.str() +
    "/Sector" + sector.str() + "/";

  string histoname = folderName + histoTag  
    + "_W" + wheel.str() 
    + "_St" + station.str() 
    + "_Sec" + sector.str(); 
  
  if(histoTag == "numberOfSegments")
    histoname = 
      folderRoot + "DT/02-Segments/Wheel" +  wheel.str() + "/" +
      histoTag  + + "_W" + wheel.str();

  return histoname;
  
}


void DTSegmentAnalysisTest::bookHistos() {

  for(int wh=-2; wh<=2; wh++){
      stringstream wheel; wheel << wh;
      string histoName =  "segmentSummary_W" + wheel.str();
      dbe->setCurrentFolder("DT/02-Segments");
      summaryHistos[wh] = dbe->book2D(histoName.c_str(),histoName.c_str(),12,1,13,4,1,5);
      summaryHistos[wh]->setAxisTitle("Sector",1);
      summaryHistos[wh]->setBinLabel(1,"MB1",2);
      summaryHistos[wh]->setBinLabel(2,"MB2",2);
      summaryHistos[wh]->setBinLabel(3,"MB3",2);
      summaryHistos[wh]->setBinLabel(4,"MB4",2);

      if(detailedAnalysis){
	for(int sect=1; sect<=14; sect++){
	  stringstream sector; sector << sect;
	  string chi2HistoName =  "chi2BadSegmPercentual_W" + wheel.str() + "_Sec" + sector.str();
	  dbe->setCurrentFolder("DT/Tests/Segments/Wheel" + wheel.str() +
				"/Sector" + sector.str());
	  chi2Histos[make_pair(wh,sect)] = dbe->book1D(chi2HistoName.c_str(),chi2HistoName.c_str(),4,1,5);
	  chi2Histos[make_pair(wh,sect)]->setBinLabel(1,"MB1");
	  chi2Histos[make_pair(wh,sect)]->setBinLabel(2,"MB2");
	  chi2Histos[make_pair(wh,sect)]->setBinLabel(3,"MB3");
	  chi2Histos[make_pair(wh,sect)]->setBinLabel(4,"MB4");
	  
	  string segmHistoName =  "residualsOnSegmRecHitNumber_W" + wheel.str() + "_Sec" + sector.str();
	  segmRecHitHistos[make_pair(wh,sect)] = dbe->book1D(segmHistoName.c_str(),segmHistoName.c_str(),4,1,5);
	  segmRecHitHistos[make_pair(wh,sect)]->setBinLabel(1,"MB1");
	  segmRecHitHistos[make_pair(wh,sect)]->setBinLabel(2,"MB2");
	  segmRecHitHistos[make_pair(wh,sect)]->setBinLabel(3,"MB3");
	  segmRecHitHistos[make_pair(wh,sect)]->setBinLabel(4,"MB4");
	  
	}
      }
  }
  
  string histoName =  "segmentSummary";
  dbe->setCurrentFolder("DT/02-Segments");
  summaryHistos[3] = dbe->book2D(histoName.c_str(),histoName.c_str(),12,1,13,5,-2,3);
  summaryHistos[3]->setAxisTitle("Sector",1);
  summaryHistos[3]->setAxisTitle("Wheel",2); 

}
  

  

