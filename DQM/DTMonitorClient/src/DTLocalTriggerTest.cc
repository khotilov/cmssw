
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2007/09/20 07:46:22 $
 *  $Revision: 1.10 $
 *  \author C. Battilana S. Marcellini - INFN Bologna
 */


// This class header
#include "DQM/DTMonitorClient/src/DTLocalTriggerTest.h"

// Framework headers
#include "FWCore/Framework/interface/EventSetup.h"
#include "DQMServices/Core/interface/MonitorElementBaseT.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

//C++ headers
#include <iostream>
#include <sstream>

using namespace edm;
using namespace std;


DTLocalTriggerTest::DTLocalTriggerTest(const edm::ParameterSet& ps){

  edm::LogVerbatim ("localTrigger") << "[DTLocalTriggerTest]: Constructor";

  sourceFolder = ps.getUntrackedParameter<string>("folderRoot", ""); 
  hwSource = ps.getUntrackedParameter<bool>("dataFromDDU", false) ? "DDU" : "DCC" ; 
  parameters = ps;
  dbe = edm::Service<DaqMonitorBEInterface>().operator->();
  dbe->setVerbose(1);

  prescaleFactor = parameters.getUntrackedParameter<int>("diagnosticPrescale", 1);

}


DTLocalTriggerTest::~DTLocalTriggerTest(){

  edm::LogVerbatim ("localTrigger") << "[DTLocalTriggerTest]: analyzed " << nevents << " events";

}


void DTLocalTriggerTest::beginJob(const edm::EventSetup& context){

  edm::LogVerbatim ("localTrigger") << "[DTLocalTriggerTest]: BeginJob";
  nevents = 0;
  context.get<MuonGeometryRecord>().get(muonGeom);
  
}


void DTLocalTriggerTest::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  edm::LogVerbatim ("localTrigger") <<"[DTLocalTriggerTest]: Begin of LS transition";

  // Get the run number
  run = lumiSeg.run();

}


void DTLocalTriggerTest::analyze(const edm::Event& e, const edm::EventSetup& context){

  nevents++;
  edm::LogVerbatim ("localTrigger") << "[DTLocalTriggerTest]: "<<nevents<<" events";

}


void DTLocalTriggerTest::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {
  
  // counts number of updats (online mode) or number of events (standalone mode)
  //nevents++;
  // if running in standalone perform diagnostic only after a reasonalbe amount of events
  //if ( parameters.getUntrackedParameter<bool>("runningStandalone", false) && 
  //   nevents%parameters.getUntrackedParameter<int>("diagnosticPrescale", 1000) != 0 ) return;


  edm::LogVerbatim ("localTrigger") <<"[DTLocalTriggerTest]: End of LS transition, performing the DQM client operation";

  // counts number of lumiSegs 
  nLumiSegs = lumiSeg.id().luminosityBlock();

  // prescale factor
  if ( nLumiSegs%prescaleFactor != 0 ) return;

  edm::LogVerbatim ("localTrigger") <<"[DTLocalTriggerTest]: "<<nLumiSegs<<" updates";
  

  // Loop over the TriggerUnits
  for (int stat=1; stat<=4; ++stat){
    for (int wh=-2; wh<=2; ++wh){
      for (int sect=1; sect<=12; ++sect){
	DTChamberId chId(wh,stat,sect);
	int sector_id = (wh+3)+(sect-1)*5;
	uint32_t indexCh = chId.rawId();

	// Perform DCC/DDU common plot analysis (Phi ones)
	TH2F * BXvsQual    = getHisto<TH2F>(dbe->get(getMEName("BXvsQual","LocalTriggerPhi", chId)));
	TH1F * BestQual    = getHisto<TH1F>(dbe->get(getMEName("BestQual","LocalTriggerPhi", chId)));
	TH2F * Flag1stvsBX = getHisto<TH2F>(dbe->get(getMEName("Flag1stvsBX","LocalTriggerPhi", chId)));
	if (BXvsQual && Flag1stvsBX && BestQual) {
	      
	  TH1D* BXHH    = BXvsQual->ProjectionY("",7,7,"");
	  TH1D* Flag1st = Flag1stvsBX->ProjectionY();
	  int BXOK_bin = BXHH->GetMaximumBin();
	  double BX_OK =  BXvsQual->GetYaxis()->GetBinCenter(BXOK_bin);
	  double trigsFlag2nd = Flag1st->GetBinContent(2);
	  double trigs = Flag1st->GetEntries();
	  double besttrigs = BestQual->GetEntries();
	  double besttrigsCorr = 0;
	  for (int i=5;i<=7;++i)
	    besttrigsCorr+=BestQual->GetBinContent(i);
	      
	  if( secME[sector_id].find("CorrectBX_Phi") == secME[sector_id].end() ){
	    bookSectorHistos(wh,sect,"LocalTriggerPhi","CorrectBX_Phi");
	    bookSectorHistos(wh,sect,"LocalTriggerPhi","CorrFraction_Phi");
	    bookSectorHistos(wh,sect,"LocalTriggerPhi","2ndFraction_Phi");
	  }
	  std::map<std::string,MonitorElement*> innerME = secME[sector_id];
	  innerME.find("CorrectBX_Phi")->second->setBinContent(stat,BX_OK);
	  innerME.find("CorrFraction_Phi")->second->setBinContent(stat,besttrigsCorr/besttrigs);
	  innerME.find("2ndFraction_Phi")->second->setBinContent(stat,trigsFlag2nd/trigs);
	    
	}

	// Perform analysis on DCC exclusive plots (Phi)	  
	TH2F * QualvsPhirad  = getHisto<TH2F>(dbe->get(getMEName("QualvsPhirad","LocalTriggerPhi", chId)));
	TH2F * QualvsPhibend = getHisto<TH2F>(dbe->get(getMEName("QualvsPhibend","LocalTriggerPhi", chId)));
	if (QualvsPhirad && QualvsPhibend) {
	      
	  TH1D* phiR = QualvsPhirad->ProjectionX();
	  TH1D* phiB = QualvsPhibend->ProjectionX();

	  if( chambME[indexCh].find("TrigDirection_Phi") == chambME[indexCh].end() ){
	    bookChambHistos(chId,"TrigDirection_Phi");
	    bookChambHistos(chId,"TrigPosition_Phi");
	  }
	  std::map<std::string,MonitorElement*> innerME = chambME[indexCh];
	  for (int i=-1;i<(phiB->GetNbinsX()+1);i++)
	    innerME.find("TrigDirection_Phi")->second->setBinContent(i,phiB->GetBinContent(i));
	  for (int i=-1;i<(phiR->GetNbinsX()+1);i++)
	    innerME.find("TrigPosition_Phi")->second->setBinContent(i,phiR->GetBinContent(i));
	     
	}

	// Perform DCC/DDU common plot analysis (Theta ones)	    
	TH2F * ThetaBXvsQual = getHisto<TH2F>(dbe->get(getMEName("ThetaBXvsQual","LocalTriggerTheta", chId)));
	TH1F * ThetaBestQual = getHisto<TH1F>(dbe->get(getMEName("ThetaBestQual","LocalTriggerTheta", chId)));
	
	if (ThetaBXvsQual && ThetaBestQual) {
	      
	  TH1D* projBXH   = ThetaBXvsQual->ProjectionY("",4,4,"");
	  int    BXOK_bin = projBXH->GetMaximumBin();
	  double BX_OK    = ThetaBXvsQual->GetYaxis()->GetBinCenter(BXOK_bin);
	  double trigs    = ThetaBestQual->GetEntries(); 
	  double trigsH   = ThetaBestQual->GetBinContent(4);
	      
	  if( secME[sector_id].find("HFraction_Theta") == secME[sector_id].end() ){
	    bookSectorHistos(wh,sect,"LocalTriggerTheta","CorrectBX_Theta");
	    bookSectorHistos(wh,sect,"LocalTriggerTheta","HFraction_Theta");
	  }
	  std::map<std::string,MonitorElement*> innerME = secME.find(sector_id)->second;
	  innerME.find("CorrectBX_Theta")->second->setBinContent(stat,BX_OK);
	  innerME.find("HFraction_Theta")->second->setBinContent(stat,trigsH/trigs);
	    
	}

	// Perform Efficiency analysis (Phi+Segments 2D)
	TH2F * TrackPosvsAngle            = getHisto<TH2F>(dbe->get(getMEName("TrackPosvsAngle","Segment", chId)));
	TH2F * TrackPosvsAngleandTrig     = getHisto<TH2F>(dbe->get(getMEName("TrackPosvsAngleandTrig","Segment", chId)));
	TH2F * TrackPosvsAngleandTrigHHHL = getHisto<TH2F>(dbe->get(getMEName("TrackPosvsAngleandTrigHHHL","Segment", chId)));
	    
	if (TrackPosvsAngle && TrackPosvsAngleandTrig && TrackPosvsAngleandTrigHHHL) {
	      
	  // Fill client histos
	  if( chambME[indexCh].find("TrigEffAngle_Phi") == chambME[indexCh].end()){
	    bookChambHistos(chId,"TrigEffPosvsAngle_Phi");
	    bookChambHistos(chId,"TrigEffPosvsAngleHHHL_Phi");
	    bookChambHistos(chId,"TrigEffPos_Phi");
	    bookChambHistos(chId,"TrigEffPosHHHL_Phi");
	    bookChambHistos(chId,"TrigEffAngle_Phi");
	    bookChambHistos(chId,"TrigEffAngleHHHL_Phi");
	  }
	  if( secME[sector_id].find("TrigEff_Phi") == secME[sector_id].end() ){
	    bookSectorHistos(wh,sect,"TriggerAndSeg","TrigEff_Phi");  
	  }

	  std::map<std::string,MonitorElement*> innerME = secME[sector_id];
	  TH1D* TrackPos               = TrackPosvsAngle->ProjectionY();
	  TH1D* TrackAngle             = TrackPosvsAngle->ProjectionX();
	  TH1D* TrackPosandTrig        = TrackPosvsAngleandTrig->ProjectionY();
	  TH1D* TrackAngleandTrig      = TrackPosvsAngleandTrig->ProjectionX();
	  TH1D* TrackPosandTrigHHHL    = TrackPosvsAngleandTrigHHHL->ProjectionY();
	  TH1D* TrackAngleandTrigHHHL  = TrackPosvsAngleandTrigHHHL->ProjectionX();
	  
 	  MonitorElement* globalEff = innerME.find("TrigEff_Phi")->second;
 	  float binEff = float(TrackPosandTrig->GetEntries())/TrackPos->GetEntries();
 	  float binErr = sqrt(binEff*(1-binEff)/TrackPos->GetEntries());
 	  globalEff->setBinContent(stat,binEff);
 	  globalEff->setBinError(stat,binErr);

	  innerME = chambME[indexCh];
	  makeEfficiencyME(TrackPosandTrig,TrackPos,innerME.find("TrigEffPos_Phi")->second);
	  makeEfficiencyME(TrackPosandTrigHHHL,TrackPos,innerME.find("TrigEffPosHHHL_Phi")->second);
	  makeEfficiencyME(TrackAngleandTrig,TrackAngle,innerME.find("TrigEffAngle_Phi")->second);
	  makeEfficiencyME(TrackAngleandTrigHHHL,TrackAngle,innerME.find("TrigEffAngleHHHL_Phi")->second);
	  makeEfficiencyME2D(TrackPosvsAngleandTrig,TrackPosvsAngle,innerME.find("TrigEffPosvsAngle_Phi")->second);
 	  makeEfficiencyME2D(TrackPosvsAngleandTrigHHHL,TrackPosvsAngle,innerME.find("TrigEffPosvsAngleHHHL_Phi")->second);
	     
	}
	
 	// Perform Efficiency analysis (Theta+Segments 2D)
	TH2F * TrackThetaPosvsAngle            = getHisto<TH2F>(dbe->get(getMEName("TrackThetaPosvsAngle","Segment", chId)));
	TH2F * TrackThetaPosvsAngleandTrig     = getHisto<TH2F>(dbe->get(getMEName("TrackThetaPosvsAngleandTrig","Segment", chId)));
	TH2F * TrackThetaPosvsAngleandTrigH    = getHisto<TH2F>(dbe->get(getMEName("TrackThetaPosvsAngleandTrigH","Segment", chId)));
	    
	if (TrackThetaPosvsAngle && 
	    TrackThetaPosvsAngleandTrig && 
	    TrackThetaPosvsAngleandTrigH) {
	      
	  // Fill client histos
 	  if( chambME[indexCh].find("TrigEffAngle_Theta") == chambME[indexCh].end()){
 	    bookChambHistos(chId,"TrigEffPosvsAngle_Theta");
 	    bookChambHistos(chId,"TrigEffPosvsAngleH_Theta");
	    bookChambHistos(chId,"TrigEffPos_Theta");
	    bookChambHistos(chId,"TrigEffPosH_Theta");
	    bookChambHistos(chId,"TrigEffAngle_Theta");
	    bookChambHistos(chId,"TrigEffAngleH_Theta");
	  }
	  if( secME[sector_id].find("TrigEff_Theta") == secME[sector_id].end() ){
	    bookSectorHistos(wh,sect,"TriggerAndSeg","TrigEff_Theta");  
	  }

	  std::map<std::string,MonitorElement*> innerME = secME[sector_id];
	  TH1D* TrackThetaPos               = TrackThetaPosvsAngle->ProjectionY();
	  TH1D* TrackThetaAngle             = TrackThetaPosvsAngle->ProjectionX();
	  TH1D* TrackThetaPosandTrig        = TrackThetaPosvsAngleandTrig->ProjectionY();
	  TH1D* TrackThetaAngleandTrig      = TrackThetaPosvsAngleandTrig->ProjectionX();
	  TH1D* TrackThetaPosandTrigH       = TrackThetaPosvsAngleandTrigH->ProjectionY();
	  TH1D* TrackThetaAngleandTrigH     = TrackThetaPosvsAngleandTrigH->ProjectionX();
	  
 	  MonitorElement* globalEff = innerME.find("TrigEff_Theta")->second;
 	  float binEff = float(TrackThetaPosandTrig->GetEntries())/TrackThetaPos->GetEntries();
 	  float binErr = sqrt(binEff*(1-binEff)/TrackThetaPos->GetEntries());
 	  globalEff->setBinContent(stat,binEff);
 	  globalEff->setBinError(stat,binErr);

	  innerME = chambME[indexCh];
	  makeEfficiencyME(TrackThetaPosandTrig,TrackThetaPos,innerME.find("TrigEffPos_Theta")->second);
	  makeEfficiencyME(TrackThetaPosandTrigH,TrackThetaPos,innerME.find("TrigEffPosH_Theta")->second);
	  makeEfficiencyME(TrackThetaAngleandTrig,TrackThetaAngle,innerME.find("TrigEffAngle_Theta")->second);
	  makeEfficiencyME(TrackThetaAngleandTrigH,TrackThetaAngle,innerME.find("TrigEffAngleH_Theta")->second);
       	  makeEfficiencyME2D(TrackThetaPosvsAngleandTrig,TrackThetaPosvsAngle,innerME.find("TrigEffPosvsAngle_Theta")->second);
  	  makeEfficiencyME2D(TrackThetaPosvsAngleandTrigH,TrackThetaPosvsAngle,innerME.find("TrigEffPosvsAngleH_Theta")->second);	     
	}

      }
    }
  }	
  
  // Efficiency test (performed on chamber plots)
  for(map<uint32_t,map<string,MonitorElement*> >::const_iterator imapIt = chambME.begin();
      imapIt != chambME.end();
      ++imapIt){
    for (map<string,MonitorElement*>::const_iterator effME = (*imapIt).second.begin();
	 effME!=(*imapIt).second.end();
	 ++effME){
      if ((*effME).second->getName().find("TrigEffPos_Phi") == 0) {
	const QReport *effQReport = (*effME).second->getQReport("ChambTrigEffInRangePhi");
	if (effQReport) {
	  if (effQReport->getBadChannels().size())
	    edm::LogError ("localTrigger") << (*effME).second->getName() <<" has " << effQReport->getBadChannels().size() << " channels out of expected efficiency range";
	  edm::LogWarning ("localTrigger") << "-------" << effQReport->getMessage() << " ------- " << effQReport->getStatus();
	}
      }
      if ((*effME).second->getName().find("TrigEffPos_Theta") == 0) {
	const QReport *effQReport = (*effME).second->getQReport("ChambTrigEffInRangeTheta");
	if (effQReport) {
	  if (effQReport->getBadChannels().size())
	    edm::LogError ("localTrigger") << (*effME).second->getName() <<" has " << effQReport->getBadChannels().size() << " channels out of expected efficiency range";
	  edm::LogWarning ("localTrigger") << "-------" << effQReport->getMessage() << " ------- " << effQReport->getStatus();
	}
      }
    }
  }

  // Efficiency test (performed on wheel plots)
  for(map<int,map<string,MonitorElement*> >::const_iterator imapIt = secME.begin();
      imapIt != secME.end();
      ++imapIt){
    for (map<string,MonitorElement*>::const_iterator effME = (*imapIt).second.begin();
	 effME!=(*imapIt).second.end();
	 ++effME){
      if ((*effME).second->getName().find("TrigEff_Phi") == 0) {
	const QReport *effQReport = (*effME).second->getQReport("SectorTrigEffInRangePhi");
	if (effQReport) {
	  edm::LogWarning ("localTrigger") << "-------" << effQReport->getMessage() << " ------- " << effQReport->getStatus();
	}
      }
      if ((*effME).second->getName().find("TrigEff_Theta") == 0) {
	const QReport *effQReport = (*effME).second->getQReport("SectorTrigEffInRangeTheta");
	if (effQReport) {
	  edm::LogWarning ("localTrigger") << "-------" << effQReport->getMessage() << " ------- " << effQReport->getStatus();
	}
      }
    }
  }

}


void DTLocalTriggerTest::endJob(){

  edm::LogVerbatim ("localTrigger") << "[DTLocalTriggerTest] endjob called!";
  dbe->rmdir("DT/Tests/DTLocalTrigger");

}


void DTLocalTriggerTest::makeEfficiencyME(TH1D* numerator, TH1D* denominator, MonitorElement* result){
  
  MonitorElementT<TNamed>* efficiencyME = dynamic_cast<MonitorElementT<TNamed>*>(result);
  TH1F* efficiency = dynamic_cast<TH1F*> (efficiencyME->operator->());
  efficiency->Divide(numerator,denominator,1,1,"");
  
  int nbins = efficiency->GetNbinsX();
  for (int bin=1; bin<=nbins; ++bin){
    float error = 0;
    float bineff = efficiency->GetBinContent(bin);

    if (denominator->GetBinContent(bin)){
      error = sqrt(bineff*(1-bineff)/denominator->GetBinContent(bin));
    }
    else {
      error = 1;
      efficiency->SetBinContent(bin,1.);
    }
 
    efficiency->SetBinError(bin,error);
  }

}


void DTLocalTriggerTest::makeEfficiencyME2D(TH2F* numerator, TH2F* denominator, MonitorElement* result){
  
  MonitorElementT<TNamed>* efficiencyME = dynamic_cast<MonitorElementT<TNamed>*>(result);
  TH2F* efficiency = dynamic_cast<TH2F*> (efficiencyME->operator->());
  efficiency->Divide(numerator,denominator,1,1,"");
  
  int nbinsx = efficiency->GetNbinsX();
  int nbinsy = efficiency->GetNbinsY();
  for (int binx=1; binx<=nbinsx; ++binx){
    for (int biny=1; biny<=nbinsy; ++biny){
      float error = 0;
      float bineff = efficiency->GetBinContent(binx,biny);

      if (denominator->GetBinContent(binx,biny)){
	error = sqrt(bineff*(1-bineff)/denominator->GetBinContent(binx,biny));
      }
      else {
	error = 1;
	efficiency->SetBinContent(binx,biny,0.);
      }
 
      efficiency->SetBinError(binx,biny,error);
    }
  }

}    


string DTLocalTriggerTest::getMEName(string histoTag, string subfolder, const DTChamberId & chambid) {

  stringstream wheel; wheel << chambid.wheel();
  stringstream station; station << chambid.station();
  stringstream sector; sector << chambid.sector();

  if (subfolder == "Segment" && histoTag.find("Trig") == string::npos) 
    histoTag = "SEG_" + histoTag;
  else
    histoTag = hwSource + "_" + histoTag;

  string folderName = 
    "DT/DTLocalTriggerTask/Wheel" +  wheel.str() +
    "/Sector" + sector.str() +
    "/Station" + station.str() + "/" +  subfolder + "/";  

  string histoname = sourceFolder + folderName + histoTag  
    + "_W" + wheel.str()  
    + "_Sec" + sector.str()
    + "_St" + station.str();
  
  return histoname;
  
}

template <class T>
T* DTLocalTriggerTest::getHisto(MonitorElement* me) {

  if (!me)
    return 0;
  MonitorElementT<TNamed>* meT = dynamic_cast<MonitorElementT<TNamed>*>(me);
  if (!meT)
    return 0;
  return dynamic_cast<T*> (meT->operator->());

}


void DTLocalTriggerTest::bookChambHistos(DTChamberId chambId, string htype) {
  
  stringstream wheel; wheel << chambId.wheel();
  stringstream station; station << chambId.station();	
  stringstream sector; sector << chambId.sector();

  string HistoName = htype + "_W" + wheel.str() + "_Sec" + sector.str() + "_St" + station.str();

  dbe->setCurrentFolder("DT/Tests/DTLocalTrigger/Wheel" + wheel.str() +
			"/Sector" + sector.str() +
			"/Station" + station.str());
  
  uint32_t indexChId = chambId.rawId();
  if (htype.find("TrigEffAngle_Phi") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.);
    return;
  }
  if (htype.find("TrigEffAngleHHHL_Phi") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.);
    return;
  }
  if (htype.find("TrigEffAngle_Theta") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.);
    return;
  }
  if (htype.find("TrigEffAngleH_Theta") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.);
    return;
  }
  if (htype.find("TrigPosition_Phi") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),100,-500.,500.);
    return;
  }
  if (htype.find("TrigDirection_Phi") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),200,-2.,2.);
    return;
  }
  if (htype.find("TrigEffPos_Phi") == 0 ){
    pair<float,float> range = phiRange(chambId);
    int nbins = int((range.second - range.first)/15);
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),nbins,range.first,range.second);
    return;
  }
  if (htype.find("TrigEffPosvsAngle_Phi") == 0 ){
    pair<float,float> range = phiRange(chambId);
    int nbins = int((range.second - range.first)/15);
    chambME[indexChId][htype] = dbe->book2D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.,nbins,range.first,range.second);
    return;
  }
  if (htype.find("TrigEffPosvsAngleHHHL_Phi") == 0 ){
    pair<float,float> range = phiRange(chambId);
    int nbins = int((range.second - range.first)/15);
    chambME[indexChId][htype] = dbe->book2D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.,nbins,range.first,range.second);
    return;
  }
  if (htype.find("TrigEffPosHHHL_Phi") == 0 ){
    pair<float,float> range = phiRange(chambId);
    int nbins = int((range.second - range.first)/15);
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),nbins,range.first,range.second);
    return;
  }
  if (htype.find("TrigEffPos_Theta") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),20,-117.5,117.5);
    return;
  }
  if (htype.find("TrigEffPosH_Theta") == 0){
    chambME[indexChId][htype] = dbe->book1D(HistoName.c_str(),HistoName.c_str(),20,-117.5,117.5);
    return;
  }
  if (htype.find("TrigEffPosvsAngle_Theta") == 0 ){
    chambME[indexChId][htype] = dbe->book2D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.,20,-117.5,117.5);
    return;
  }
  if (htype.find("TrigEffPosvsAngleH_Theta") == 0 ){
    chambME[indexChId][htype] = dbe->book2D(HistoName.c_str(),HistoName.c_str(),16,-40.,40.,20,-117.5,117.5);
    return;
  }

}


void DTLocalTriggerTest::bookSectorHistos(int wheel,int sector,string folder, string htype) {
  
  stringstream wh; wh << wheel;
  stringstream sc; sc << sector;
  int sectorid = (wheel+3) + (sector-1)*5;
  dbe->setCurrentFolder("DT/Tests/DTLocalTrigger/Wheel"+ wh.str()+"/Sector"+ sc.str()+"/"+folder);

  if (htype.find("Phi") != string::npos){    
    string hname = htype + "_W" + wh.str()+"_Sec" +sc.str();
    MonitorElement* me = dbe->book1D(hname.c_str(),hname.c_str(),4,1,5);
    //setLabelPh(me);
    secME[sectorid][htype] = me;
    return;
  }
  
  if (htype.find("Theta") != string::npos){
    string hname = htype + "_W" + wh.str()+ "_Sec"+sc.str();
    MonitorElement* me =dbe->book1D(hname.c_str(),hname.c_str(),3,1,4);
    //setLabelTh(me);
    secME[sectorid][htype] = me;
    return;
  }
  
}


pair<float,float> DTLocalTriggerTest::phiRange(const DTChamberId& id){

  float min,max;
  int station = id.station();
  int sector  = id.sector(); 
  int wheel   = id.wheel();
  
  const DTLayer  *layer = muonGeom->layer(DTLayerId(id,1,1));
  DTTopology topo = layer->specificTopology();
  min = topo.wirePosition(topo.firstChannel());
  max = topo.wirePosition(topo.lastChannel());

  if (station == 4){
    
    const DTLayer *layer2;
    float lposx;
    
    if (sector == 4){
      layer2  = muonGeom->layer(DTLayerId(wheel,station,13,1,1));
      lposx = layer->toLocal(layer2->position()).x();
    }
    else if (sector == 10){
      layer2 = muonGeom->layer(DTLayerId(wheel,station,14,1,1));
      lposx = layer->toLocal(layer2->position()).x();
    }
    else
      return make_pair(min,max);
    
    DTTopology topo2 = layer2->specificTopology();

    if (lposx>0){
      max = lposx*.5+topo2.wirePosition(topo2.lastChannel());
      min -= lposx*.5;
    }
    else{
      min = lposx*.5+topo2.wirePosition(topo2.firstChannel());
      max -= lposx*.5;
    }
      
  }
  
  return make_pair(min,max);

}
