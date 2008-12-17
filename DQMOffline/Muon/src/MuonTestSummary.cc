
/*
 *  See header file for a description of this class.
 *
 *  $Date: 2008/12/17 13:57:16 $
 *  $Revision: 1.9 $
 *  \author G. Mila - INFN Torino
 */


#include <DQMOffline/Muon/src/MuonTestSummary.h>

// Framework
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "DQMOffline/Muon/test/langauFit.C"
#include <string>

using namespace edm;
using namespace std;


MuonTestSummary::MuonTestSummary(const edm::ParameterSet& ps){

  dbe = Service<DQMStore>().operator->();

  // parameter initialization for kinematics test
  etaExpected = ps.getParameter<double>("etaExpected");
  phiExpected = ps.getParameter<double>("phiExpected");
  chi2Fraction = ps.getParameter<double>("chi2Fraction");
  chi2Spread = ps.getParameter<double>("chi2Spread");
  resEtaSpread_tkGlb = ps.getParameter<double>("resEtaSpread_tkGlb");
  resEtaSpread_glbSta = ps.getParameter<double>("resEtaSpread_glbSta");
  resPhiSpread_tkGlb = ps.getParameter<double>("resPhiSpread_tkGlb");
  resPhiSpread_glbSta = ps.getParameter<double>("resPhiSpread_glbSta");
  resOneOvPSpread_tkGlb = ps.getParameter<double>("resOneOvPSpread_tkGlb");
  resOneOvPSpread_glbSta = ps.getParameter<double>("resOneOvPSpread_glbSta");
  resChargeLimit_tkGlb = ps.getParameter<double>("resChargeLimit_tkGlb");
  resChargeLimit_glbSta = ps.getParameter<double>("resChargeLimit_glbSta");
  resChargeLimit_tkSta = ps.getParameter<double>("resChargeLimit_tkSta");
  numMatchedExpected = ps.getParameter<double>("numMatchedExpected");
  expMolteplicityGlb = ps.getParameter<double>("expMolteplicityGlb");
  expMolteplicityTk = ps.getParameter<double>("expMolteplicityTk");
  expMolteplicitySta = ps.getParameter<double>("expMolteplicitySta");
}

MuonTestSummary::~MuonTestSummary(){}

void MuonTestSummary::beginJob(const edm::EventSetup& context){

  metname = "muonTestSummary";
  LogTrace(metname)<<"[MuonTestSummary] Histo booking";

  // book the summary histos
  dbe->setCurrentFolder("Muons/TestSummary"); 

  // kinematics test report
  kinematicsSummaryMap = dbe->book2D("kinematicsSummaryMap","Kinematics test summary",5,1,6,3,1,4);
  kinematicsSummaryMap->setAxisTitle("track monitored",1);
  kinematicsSummaryMap->setBinLabel(1,"GLB",1);
  kinematicsSummaryMap->setBinLabel(2,"TKfromGLB",1);
  kinematicsSummaryMap->setBinLabel(3,"STAfromGLB",1);
  kinematicsSummaryMap->setBinLabel(4,"TK",1);
  kinematicsSummaryMap->setBinLabel(5,"STA",1);
  kinematicsSummaryMap->setAxisTitle("parameter tested",2);
  kinematicsSummaryMap->setBinLabel(1,"#chi^{2}",2);
  kinematicsSummaryMap->setBinLabel(2,"#eta",2);
  kinematicsSummaryMap->setBinLabel(3,"#phi",2);
   
  // residuals test report
  residualsSummaryMap = dbe->book2D("residualsSummaryMap","Residuals test summary",3,1,4,4,1,5);
  residualsSummaryMap->setAxisTitle("residuals",1);
  residualsSummaryMap->setBinLabel(1,"TK-GLB",1);
  residualsSummaryMap->setBinLabel(2,"GLB-STA",1);
  residualsSummaryMap->setBinLabel(3,"TK-STA",1);
  residualsSummaryMap->setAxisTitle("parameter tested",2);
  residualsSummaryMap->setBinLabel(1,"#eta",2);
  residualsSummaryMap->setBinLabel(2,"#phi",2);
  residualsSummaryMap->setBinLabel(3,"1/p",2);
  residualsSummaryMap->setBinLabel(4,"q",2);

  // muonId test report
  muonIdSummaryMap = dbe->book1D("muonIdSummaryMap","muonId test summary",3,1,4);
  muonIdSummaryMap->setAxisTitle("test");
  muonIdSummaryMap->setBinLabel(1,"#matchCh");
  muonIdSummaryMap->setBinLabel(2,"#assSegm");
  muonIdSummaryMap->setBinLabel(3,"resTrackSegm");

  // energy test report
  energySummaryMap = dbe->book2D("energySummaryMap","Energy deposits test summary",3,1,4,3,1,4);
  energySummaryMap->setAxisTitle("muons",1);
  energySummaryMap->setBinLabel(1,"GLB",1);
  energySummaryMap->setBinLabel(2,"TK",1);
  energySummaryMap->setBinLabel(3,"STA",1);
  energySummaryMap->setAxisTitle("calorimeter tested",2);
  energySummaryMap->setBinLabel(1,"ECAL",2);
  energySummaryMap->setBinLabel(2,"HAD",2);
  energySummaryMap->setBinLabel(3,"H0",2);

  // molteplicity tests report
  molteplicitySummaryMap = dbe->book1D("molteplicitySummaryMap","muon molteplicity test summary",3,1,4);
  molteplicitySummaryMap->setAxisTitle("muon");
  molteplicitySummaryMap->setBinLabel(1,"GLB");
  molteplicitySummaryMap->setBinLabel(2,"TK");
  molteplicitySummaryMap->setBinLabel(3,"STA");

  
  // summary test report
  dbe->setCurrentFolder("Muons/EventInfo"); 
  summaryReport = dbe->bookFloat("reportSummary");

  summaryReportMap = dbe->book2D("reportSummaryMap","Muon Report Summary Map",3,1,4,7,1,8);
  summaryReportMap->setAxisTitle("muons",1);
  summaryReportMap->setBinLabel(1,"GLB",1);
  summaryReportMap->setBinLabel(2,"TK",1);
  summaryReportMap->setBinLabel(3,"STA",1);
  summaryReportMap->setAxisTitle("test",2);
  summaryReportMap->setBinLabel(1,"#chi^{2}/Df",2);
  summaryReportMap->setBinLabel(2,"#eta",2);
  summaryReportMap->setBinLabel(3,"#phi",2);
  summaryReportMap->setBinLabel(4,"residuals",2);
  summaryReportMap->setBinLabel(5,"muonId",2);
  summaryReportMap->setBinLabel(6,"energyDeposits",2);
  summaryReportMap->setBinLabel(7,"molteplicity",2);

  dbe->setCurrentFolder("Muons/EventInfo/reportSummaryContents");
  theSummaryContents.push_back(dbe->bookFloat("kinematics_GLB"));
  theSummaryContents.push_back(dbe->bookFloat("residuals_GLB"));
  theSummaryContents.push_back(dbe->bookFloat("kinematics_TK"));
  theSummaryContents.push_back(dbe->bookFloat("muonId_TK"));
  theSummaryContents.push_back(dbe->bookFloat("GLB"));
  theSummaryContents.push_back(dbe->bookFloat("TK"));
  theSummaryContents.push_back(dbe->bookFloat("STA"));
  theSummaryContents.push_back(dbe->bookFloat("energyDeposits"));
  theSummaryContents.push_back(dbe->bookFloat("molteciplicity"));

}


void MuonTestSummary::beginLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  // initialisation of histo bins
  for(int xBin=1; xBin<=5; xBin++){
    for(int yBin=1; yBin<=3; yBin++){
      kinematicsSummaryMap->Fill(xBin,yBin,0);
    }
  }
  for(int xBin=1; xBin<=3; xBin++){
    for(int yBin=1; yBin<=2; yBin++){
      residualsSummaryMap->Fill(xBin,yBin,0);
    }
  }
  for(int xBin=1; xBin<=3; xBin++){
    muonIdSummaryMap->Fill(xBin,0);
  }
  for(int xBin=1; xBin<=3; xBin++){
    for(int yBin=1; yBin<=3; yBin++){
      energySummaryMap->Fill(xBin,yBin,0);
    }
  }
  for(int xBin=1; xBin<=3; xBin++){
    molteplicitySummaryMap->Fill(xBin,0);
  }
  for(int xBin=1; xBin<=3; xBin++){
    for(int yBin=1; yBin<=7; yBin++){
      summaryReportMap->Fill(xBin,yBin,0);
    }
  }
}



void MuonTestSummary::endLuminosityBlock(LuminosityBlock const& lumiSeg, EventSetup const& context) {

  // fill the kinematics report summary
  doKinematicsTests("GlbMuon_Glb_", 1);
  doKinematicsTests("GlbMuon_Tk_", 2);
  doKinematicsTests("GlbMuon_Sta_",3);
  doKinematicsTests("TkMuon_", 4);
  doKinematicsTests("StaMuon_", 5);

  // fill the residuals report summary
  doResidualsTests("TkGlb", "eta", 1);
  doResidualsTests("GlbSta", "eta", 2);
  doResidualsTests("TkSta", "eta", 3);
  doResidualsTests("TkGlb", "phi", 1);
  doResidualsTests("GlbSta", "phi", 2);
  doResidualsTests("TkSta", "phi", 3);
  doResidualsTests("TkGlb", "oneOverp", 1);
  doResidualsTests("GlbSta", "oneOverp", 2);
  doResidualsTests("TkSta", "oneOverp", 3);
  doResidualsTests("GlbMuon", "qComparison", -1);

  // fill the muonID report summary
  doMuonIDTests();

  // fill the energy report summary
  doEnergyTests("ecalS9PointingMuDepositedEnergy_","Glb_muons", 1);
  doEnergyTests("hadS9PointingMuDepositedEnergy_", "Glb_muons", 1);
  doEnergyTests("hoS9PointingMuDepositedEnergy_", "Glb_muons", 1);
  doEnergyTests("ecalS9PointingMuDepositedEnergy_", "Tk_muons", 2);
  doEnergyTests("hadS9PointingMuDepositedEnergy_", "Tk_muons", 2);
  doEnergyTests("hoS9PointingMuDepositedEnergy_", "Tk_muons", 2);
  doEnergyTests("ecalS9PointingMuDepositedEnergy_", "Sta_muons", 3);
  doEnergyTests("hadS9PointingMuDepositedEnergy_", "Sta_muons", 3);
  doEnergyTests("hoS9PointingMuDepositedEnergy_", "Sta_muons", 3);

  // fill the molteplicity test summary
  doMolteplicityTests();
  
  // fill the final report summary
  summaryReportMap->setBinContent(1,1,double(kinematicsSummaryMap->getBinContent(1,1)+kinematicsSummaryMap->getBinContent(2,1)+kinematicsSummaryMap->getBinContent(3,1))/3.0);
  summaryReportMap->setBinContent(2,1,kinematicsSummaryMap->getBinContent(4,1));
  summaryReportMap->setBinContent(3,1,kinematicsSummaryMap->getBinContent(5,1));
  summaryReportMap->setBinContent(1,2,double(kinematicsSummaryMap->getBinContent(1,2)+kinematicsSummaryMap->getBinContent(2,2)+kinematicsSummaryMap->getBinContent(3,2))/3.0);
  summaryReportMap->setBinContent(2,2,kinematicsSummaryMap->getBinContent(4,2));
  summaryReportMap->setBinContent(3,2,kinematicsSummaryMap->getBinContent(5,2));
  summaryReportMap->setBinContent(1,3,double(kinematicsSummaryMap->getBinContent(1,3)+kinematicsSummaryMap->getBinContent(2,3)+kinematicsSummaryMap->getBinContent(3,3))/3.0);
  summaryReportMap->setBinContent(2,3,kinematicsSummaryMap->getBinContent(4,3));
  summaryReportMap->setBinContent(3,3,kinematicsSummaryMap->getBinContent(5,3));
  summaryReportMap->setBinContent(1,4,double(residualsSummaryMap->getBinContent(1,1)+residualsSummaryMap->getBinContent(1,2)+residualsSummaryMap->getBinContent(1,3)+residualsSummaryMap->getBinContent(1,4)+residualsSummaryMap->getBinContent(2,1)+residualsSummaryMap->getBinContent(2,2)+residualsSummaryMap->getBinContent(2,3)+residualsSummaryMap->getBinContent(2,4)+residualsSummaryMap->getBinContent(3,1)+residualsSummaryMap->getBinContent(3,2)+residualsSummaryMap->getBinContent(3,3)+residualsSummaryMap->getBinContent(3,4))/12);
  summaryReportMap->setBinContent(2,4,-1.0/6.0);
  summaryReportMap->setBinContent(3,4,-1.0/6.0);
  summaryReportMap->setBinContent(1,5,-1.0/6.0);
  summaryReportMap->setBinContent(2,5,double(muonIdSummaryMap->getBinContent(1)+muonIdSummaryMap->getBinContent(2)+muonIdSummaryMap->getBinContent(3))/3.0);
  summaryReportMap->setBinContent(3,5,-1.0/6.0);
  summaryReportMap->setBinContent(1,6,double(energySummaryMap->getBinContent(1,1)+energySummaryMap->getBinContent(1,2)+energySummaryMap->getBinContent(1,3))/3.0);
  summaryReportMap->setBinContent(2,6,double(energySummaryMap->getBinContent(2,1)+energySummaryMap->getBinContent(2,2)+energySummaryMap->getBinContent(2,3))/3.0);
  summaryReportMap->setBinContent(3,6,double(energySummaryMap->getBinContent(3,1)+energySummaryMap->getBinContent(3,2)+energySummaryMap->getBinContent(3,3))/3.0);
  summaryReportMap->setBinContent(1,7,molteplicitySummaryMap->getBinContent(1));
  summaryReportMap->setBinContent(2,7,molteplicitySummaryMap->getBinContent(2));
  summaryReportMap->setBinContent(3,7,molteplicitySummaryMap->getBinContent(3));

  double kinematics_GLB = double(summaryReportMap->getBinContent(1,1)+summaryReportMap->getBinContent(1,2)+summaryReportMap->getBinContent(1,3))/3.0;
  theSummaryContents[0]->Fill(kinematics_GLB);
  double residuals_GLB = summaryReportMap->getBinContent(1,4);
  theSummaryContents[1]->Fill(residuals_GLB);
  double kinematics_TK = double(summaryReportMap->getBinContent(2,1)+summaryReportMap->getBinContent(2,2)+summaryReportMap->getBinContent(2,3))/3.0;
  theSummaryContents[2]->Fill(kinematics_TK);
  double muonId_TK = summaryReportMap->getBinContent(2,5);
  theSummaryContents[3]->Fill(muonId_TK);
  double GLB = (kinematics_GLB+residuals_GLB)/2.0;
  theSummaryContents[4]->Fill(GLB);
  double TK = double(kinematics_TK+muonId_TK)/2.0;
  theSummaryContents[5]->Fill(TK);
  double STA = double(summaryReportMap->getBinContent(3,1)+summaryReportMap->getBinContent(3,2)+summaryReportMap->getBinContent(3,3))/3.0;
  theSummaryContents[6]->Fill(STA);
  double energyDeposits = double(summaryReportMap->getBinContent(1,6)+summaryReportMap->getBinContent(2,6)+summaryReportMap->getBinContent(3,6))/3.0;
  theSummaryContents[7]->Fill(energyDeposits);
  double molteplicity = double(summaryReportMap->getBinContent(1,7)+summaryReportMap->getBinContent(2,7)+summaryReportMap->getBinContent(3,7))/3.0;
  theSummaryContents[8]->Fill(molteplicity);

  summaryReport->Fill((GLB+TK*STA+energyDeposits+molteplicity)/5.0);

 }


void MuonTestSummary::doKinematicsTests(string muonType, int bin){
  

  // chi2 test
  string path = "Muons/MuonRecoAnalyzer/" + muonType + "chi2OverDf";
  MonitorElement * chi2Histo = dbe->get(path);

  if(chi2Histo){

    TH1F * chi2Histo_root = chi2Histo->getTH1F();
    int maxBin = chi2Histo_root->GetMaximumBin();
    if(chi2Histo_root->Integral(maxBin+1,chi2Histo_root->GetNbinsX())!=0){
      double fraction = double(chi2Histo_root->Integral(1,maxBin))/double(chi2Histo_root->Integral(maxBin+1,chi2Histo_root->GetNbinsX()));
      LogTrace(metname)<<"chi2 fraction for "<<muonType<<" : "<<fraction<<endl;
      if(fraction>(chi2Fraction-chi2Spread) && fraction<(chi2Fraction+chi2Spread))
	kinematicsSummaryMap->setBinContent(bin,1,1);
      else
	kinematicsSummaryMap->setBinContent(bin,1,0);
    }
   }



  // pseudorapidity test
  path = "Muons/MuonRecoAnalyzer/" + muonType + "eta";
  MonitorElement * etaHisto = dbe->get(path);
  
  if(etaHisto){

    TH1F * etaHisto_root = etaHisto->getTH1F();
    double binSize = (etaHisto_root->GetXaxis()->GetXmax()-etaHisto_root->GetXaxis()->GetXmin())/etaHisto_root->GetNbinsX();
    int binZero = int((0-etaHisto_root->GetXaxis()->GetXmin())/binSize);
    if(etaHisto_root->Integral(1,binZero-1)!=0 && etaHisto_root->Integral(binZero,etaHisto_root->GetNbinsX())!=0){
      double symmetryFactor = 
	double(etaHisto_root->Integral(1,binZero-1)) / double(etaHisto_root->Integral(binZero,etaHisto_root->GetNbinsX()));
      double errSymmetryFactor =
	symmetryFactor*sqrt(1.0/double(etaHisto_root->Integral(1,binZero-1)) + 1.0/double(etaHisto_root->Integral(binZero,etaHisto_root->GetNbinsX())));
      LogTrace(metname)<<"eta symmetryFactor for "<<muonType<<" : "<<symmetryFactor<<endl;
      LogTrace(metname)<<"eta errSymmetryFactor for "<<muonType<<" : "<<errSymmetryFactor<<endl;
      double tParameter;
      if((symmetryFactor-etaExpected)>0) tParameter=double(symmetryFactor-etaExpected)/errSymmetryFactor;
      else tParameter=double(-symmetryFactor+etaExpected)/errSymmetryFactor;
      if (tParameter<1.95) //2sigma rejection
	kinematicsSummaryMap->setBinContent(bin,2,1);
      else
	kinematicsSummaryMap->setBinContent(bin,2,0);
    }
  }


  // phi test
  path = "Muons/MuonRecoAnalyzer/" + muonType + "phi";
  MonitorElement * phiHisto = dbe->get(path);

  if(phiHisto){

    TH1F * phiHisto_root = phiHisto->getTH1F();
    double binSize = (phiHisto_root->GetXaxis()->GetXmax()-phiHisto_root->GetXaxis()->GetXmin())/phiHisto_root->GetNbinsX();
    int binZero = int((0-phiHisto_root->GetXaxis()->GetXmin())/binSize);
    if(phiHisto_root->Integral(binZero+1,phiHisto_root->GetNbinsX())!=0 && phiHisto_root->Integral(1,binZero)!=0){
      double symmetryFactor = 
	double(phiHisto_root->Integral(binZero+1,phiHisto_root->GetNbinsX())) / double(phiHisto_root->Integral(1,binZero));
      double errSymmetryFactor = 
	symmetryFactor*sqrt(1.0/double(phiHisto_root->Integral(binZero+1,phiHisto_root->GetNbinsX())) + 1.0/double(phiHisto_root->Integral(1,binZero)));
      LogTrace(metname)<<"phi symmetryFactor for "<<muonType<<" : "<<symmetryFactor<<endl;
      LogTrace(metname)<<"phi errSymmetryFactor for "<<muonType<<" : "<<errSymmetryFactor<<endl;
      double tParameter;
      if((symmetryFactor-phiExpected)>0) tParameter=double(symmetryFactor-phiExpected)/errSymmetryFactor;
      else tParameter=double(-symmetryFactor+phiExpected)/errSymmetryFactor;
      if (tParameter<1.95) //2sigma rejection
	kinematicsSummaryMap->setBinContent(bin,3,1);
      else
	kinematicsSummaryMap->setBinContent(bin,3,0);
    }
  }

}




void MuonTestSummary::doResidualsTests(string type, string parameter, int bin){

  // residuals test
  if(type!="GlbMuon"){
    string path = "Muons/MuonRecoAnalyzer/Res_" + type + "_" + parameter;
    MonitorElement * residualsHisto = dbe->get(path);
    
    if(residualsHisto){
      
      // Gaussian Fit
      float statMean = residualsHisto->getMean(1);
      float statSigma = residualsHisto->getRMS(1);
      Double_t mean = -1;
      Double_t sigma = -1;
      TH1F * histo_root = residualsHisto->getTH1F();
      if(histo_root->GetEntries()>20){
	TF1 *gfit = new TF1("Gaussian","gaus",(statMean-(2*statSigma)),(statMean+(2*statSigma)));
	try {
	  histo_root->Fit(gfit, "Q");
	} catch (...) {
	  edm::LogError (metname)<< "[MuonTestSummary]: Exception when fitting Res_"<<type<<"_"<<parameter;
	}
	if(gfit){
	  mean = gfit->GetParameter(1); 
	  sigma = gfit->GetParameter(2);
	  LogTrace(metname)<<"mean: "<<mean<<"for Res_"<<type<<"_"<<parameter<<endl;
	  LogTrace(metname)<<"sigma: "<<sigma<<"for Res_"<<type<<"_"<<parameter<<endl;
	}
      }
      else{
	LogTrace(metname) << "[MuonTestSummary]: Test of  Res_"<<type<<"_"<<parameter<< " not performed because # entries < 20 ";
      }
      
      if(sigma!=-1 && parameter=="eta" && type=="TkGlb"){
	if(sigma<resEtaSpread_tkGlb) residualsSummaryMap->setBinContent(bin, 1, 1);
	else residualsSummaryMap->setBinContent(bin, 1, 0);
      }
      if(sigma!=-1 && parameter=="eta" && (type=="GlbSta" || type=="TkSta")) {
	if(sigma<resEtaSpread_glbSta) residualsSummaryMap->setBinContent(bin, 1, 1);
	else residualsSummaryMap->setBinContent(bin, 1, 0);
      }
      if(sigma!=-1 && parameter=="phi" && type=="TkGlb"){
	if(sigma<resPhiSpread_tkGlb) residualsSummaryMap->setBinContent(bin, 2, 1);     
	else residualsSummaryMap->setBinContent(bin, 2, 0);
      }
      if(sigma!=-1 && parameter=="phi" && (type=="GlbSta" || type=="TkSta")){ 
	if(sigma<resPhiSpread_glbSta) residualsSummaryMap->setBinContent(bin, 2, 1); 
	else residualsSummaryMap->setBinContent(bin, 2, 0); 
      }
      if(sigma!=-1 && parameter=="oneOverp" && type=="TkGlb"){
	if(sigma<resOneOvPSpread_tkGlb) residualsSummaryMap->setBinContent(bin, 3, 1);
	else residualsSummaryMap->setBinContent(bin, 3, 0);
      }
      if(sigma!=-1 && parameter=="oneOverp" && (type=="GlbSta" || type=="TkSta")) {
	if(sigma<resOneOvPSpread_glbSta) residualsSummaryMap->setBinContent(bin, 3, 1);
	else residualsSummaryMap->setBinContent(bin, 3, 0);
      }
    }
  }
  else{
    string path = "Muons/MuonRecoAnalyzer/" + type + "_" + parameter;
    MonitorElement * residualsHisto = dbe->get(path);
    
    if(residualsHisto){
      if((residualsHisto->getBinContent(3)+residualsHisto->getBinContent(4))!=0){
	LogTrace(metname)<<"charge comparison TkGlb: "<<residualsHisto->getBinContent(4)/double(residualsHisto->getBinContent(3)+residualsHisto->getBinContent(4))<<endl;
	if(residualsHisto->getBinContent(4)/double(residualsHisto->getBinContent(3)+residualsHisto->getBinContent(4)) < resChargeLimit_tkGlb)
	  residualsSummaryMap->setBinContent(1, 4, 1);
	else
	  residualsSummaryMap->setBinContent(1, 4, 0);
      }
      if((residualsHisto->getBinContent(1)+residualsHisto->getBinContent(2))!=0){
	LogTrace(metname)<<"charge comparison GlbSta: "<<residualsHisto->getBinContent(2)/double(residualsHisto->getBinContent(1)+residualsHisto->getBinContent(2))<<endl;
	if(residualsHisto->getBinContent(2)/double(residualsHisto->getBinContent(1)+residualsHisto->getBinContent(2)) < resChargeLimit_glbSta)
	  residualsSummaryMap->setBinContent(2, 4, 1);
	else
	  residualsSummaryMap->setBinContent(2, 4, 0);
      }
      if(residualsHisto->getBinContent(5)+residualsHisto->getBinContent(6)!=0){
	LogTrace(metname)<<"charge comparison TkSta: "<<residualsHisto->getBinContent(6)/double(residualsHisto->getBinContent(5)+residualsHisto->getBinContent(6))<<endl;
	if(residualsHisto->getBinContent(6)/double(residualsHisto->getBinContent(5)+residualsHisto->getBinContent(6)) < resChargeLimit_tkSta)
	  residualsSummaryMap->setBinContent(3, 4, 1);
	else
	  residualsSummaryMap->setBinContent(3, 4, 0);
      }
    }
  }
  
}

void MuonTestSummary::doMuonIDTests(){

  // num matches test
  string path = "Muons/MuonIdDQM/TrackerMuons/hNumMatches";
  MonitorElement * matchesHisto = dbe->get(path);

  if(matchesHisto){
    TH1F * matchesHisto_root = matchesHisto->getTH1F();
    if(matchesHisto_root->GetMaximumBin() == numMatchedExpected || matchesHisto_root->GetMaximumBin() == numMatchedExpected+1)
      muonIdSummaryMap->setBinContent(1,1);
    else
      muonIdSummaryMap->setBinContent(1,0);
  }

  
  // num of associated segments (limits computed from simulated data)
  double numOneSegm_dt = 0;
  MonitorElement * DT1Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hDT1NumSegments");
  if(DT1Histo) numOneSegm_dt+=DT1Histo->getBinContent(2);
  MonitorElement * DT2Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hDT2NumSegments");
  if(DT2Histo) numOneSegm_dt+=DT2Histo->getBinContent(2);
  MonitorElement * DT3Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hDT3NumSegments");
  if(DT3Histo) numOneSegm_dt+=DT3Histo->getBinContent(2);
  MonitorElement * DT4Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hDT4NumSegments"); 
  if(DT4Histo) numOneSegm_dt+=DT4Histo->getBinContent(2);
  LogTrace(metname)<<"numOneSegm_dt: "<<numOneSegm_dt<<endl;
  double fraction_dt=0;
  if(numOneSegm_dt!=0){
    fraction_dt = double(DT1Histo->getEntries())/numOneSegm_dt;
    LogTrace(metname)<<"fraction_dt: "<<fraction_dt<<endl;
  }

  double numOneSegm_csc = 0;
  MonitorElement * CSC1Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hCSC1NumSegments");
  if(CSC1Histo) numOneSegm_csc+=CSC1Histo->getBinContent(2);
  MonitorElement * CSC2Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hCSC2NumSegments");
  if(CSC2Histo) numOneSegm_csc+=CSC2Histo->getBinContent(2);
  MonitorElement * CSC3Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hCSC3NumSegments");
  if(CSC3Histo) numOneSegm_csc+=CSC3Histo->getBinContent(2);
  MonitorElement * CSC4Histo = dbe->get("Muons/MuonIdDQM/TrackerMuons/hCSC4NumSegments");
  if(CSC4Histo) numOneSegm_csc+=CSC4Histo->getBinContent(2);
  LogTrace(metname)<<"numOneSegm_csc: "<<numOneSegm_csc<<endl;
  double fraction_csc=0;
  if(numOneSegm_csc!=0){
    fraction_csc = double(CSC1Histo->getEntries())/numOneSegm_csc;
    LogTrace(metname)<<"fraction_csc: "<<fraction_csc<<endl;
  }

  if((fraction_dt>0.7 && fraction_dt<0.8) && (fraction_csc>0.57 && fraction_csc<0.67))
    muonIdSummaryMap->setBinContent(2,1);
  else{
    if((fraction_dt>0.7 && fraction_dt<0.8) || (fraction_csc>0.57 && fraction_csc<0.67))
      muonIdSummaryMap->setBinContent(2,1.0/2.0);
    else
      muonIdSummaryMap->setBinContent(2,0);
  }

  // residuals test
  vector<string> resHistos;
  resHistos.push_back("hDT1Pullx");
  resHistos.push_back("hDT2Pullx");
  resHistos.push_back("hDT3Pullx");
  resHistos.push_back("hDT4Pullx");
  resHistos.push_back("hDT1Pully");
  resHistos.push_back("hDT2Pully");
  resHistos.push_back("hDT3Pully");
  resHistos.push_back("hCSC1Pullx");
  resHistos.push_back("hCSC2Pullx");
  resHistos.push_back("hCSC3Pullx");
  resHistos.push_back("hCSC4Pullx");
  resHistos.push_back("hCSC1Pully");
  resHistos.push_back("hCSC2Pully");
  resHistos.push_back("hCSC3Pully");
  resHistos.push_back("hCSC4Pully");

  double dtSigma=0; 
  double cscSigma=0;
  int numPlot_dt=0;
  int numPlot_csc=0;
  for(int name=0; name<=14; name++){   
    MonitorElement * resHisto = dbe->get("Muons/MuonIdDQM/TrackerMuons/"+resHistos[name]);
    if(resHisto){
 
      TH1F * resHisto_root = resHisto->getTH1F();
      if(resHisto_root->GetEntries()>20){
	TF1 *gfit = new TF1("Gaussian","gaus",-2,2);
	try {
	  resHisto_root->Fit(gfit, "Q");
	} catch (...) {
	  edm::LogError (metname)<< "[MuonTestSummary]: Exception when fitting "<<resHistos[name];
	}
	if(gfit){
	  double mean = gfit->GetParameter(1); 
	  double sigma = gfit->GetParameter(2);
	  LogTrace(metname)<<"meanRes: "<<mean<<" for "<<resHistos[name]<<endl;
	  LogTrace(metname)<<"sigmaRes: "<<sigma<<" for "<<resHistos[name]<<endl;
	  if(name<=6) {
	    dtSigma+=sigma;
	    numPlot_dt++;
	  }
	  else {
	    cscSigma+=sigma;
	    numPlot_csc++;
	  }
	}
      }
      else{
      LogTrace(metname) << "[MuonTestSummary]: Test of "<<resHistos[name]<< " not performed because # entries < 20 ";
    }

    }
  } // loop over residuals histos
  
  LogTrace(metname)<<"meanDtSigma: "<<dtSigma/double(numPlot_dt)<<endl;
  LogTrace(metname)<<"meanCscSigma: "<<cscSigma/double(numPlot_csc)<<endl;
  if(numPlot_dt!=0 && numPlot_csc!=0){
    if(((dtSigma/double(numPlot_dt))>0.95 && (dtSigma/double(numPlot_dt))<1.05) 
       && ((cscSigma/double(numPlot_csc))>0.95 && (cscSigma/double(numPlot_csc))<1.05))
      muonIdSummaryMap->setBinContent(3,1);
    else{
      if(((dtSigma/double(numPlot_dt))>0.95 && (dtSigma/double(numPlot_dt))<1.05) 
	 || ((cscSigma/double(numPlot_csc))>0.95 && (cscSigma/double(numPlot_csc))<1.05))
	muonIdSummaryMap->setBinContent(3,1.0/2.0);
      else
	muonIdSummaryMap->setBinContent(3,0);
    }
  }
  if(numPlot_dt!=0 && numPlot_csc!=0){
    if((dtSigma/double(numPlot_dt))>0.95 && (dtSigma/double(numPlot_dt))<1.05) 
      muonIdSummaryMap->setBinContent(3,1.0/2.0);
    else
      muonIdSummaryMap->setBinContent(3,0);
  }
  if(numPlot_dt==0 && numPlot_csc!=0){
    if((cscSigma/double(numPlot_csc))>0.95 && (cscSigma/double(numPlot_csc))<1.05)
      muonIdSummaryMap->setBinContent(3,1.0/2.0);
    else
      muonIdSummaryMap->setBinContent(3,0);
  }
  
}


void MuonTestSummary::doEnergyTests(string histname, string muonType, int binNumber){

  // num matches test
  string path = "Muons/MuonEnergyDepositAnalyzer/"+histname+muonType;
  MonitorElement * energyHisto = dbe->get(path);
  Double_t hPeak=-1, hFWHM=-1;
  if(energyHisto){
    TH1F * energyHisto_root = energyHisto->getTH1F();
    
    // Setting fit range and start values
    Double_t fitRange[2];
    Double_t startValues[4], parlimitslo[4], parlimitshi[4], fitPar[4], fitParErr[4];
    fitRange[0]=0.0;
    fitRange[1]=3.0;

    startValues[0]=0.036; startValues[1]=0.193; startValues[2]=110.0; startValues[3]=0.06;
    parlimitslo[0]=0.0; parlimitslo[1]=0.; parlimitslo[2]=1.0; parlimitslo[3]=0.;
    parlimitshi[0]=0.05; parlimitshi[1]=0.5; parlimitshi[2]=80000.0; parlimitshi[3]=0.1;
    
    Double_t chisqr;
    Int_t    ndf;
    TF1 *fit = langaufit(energyHisto_root,fitRange,startValues,parlimitslo,parlimitshi,fitPar,fitParErr,&chisqr,&ndf);
    if(fit){
      langaupro(fitPar,hPeak,hFWHM);
      LogTrace(metname)<<"hPeak from langau fit: "<<hPeak<<" for: "<<histname+muonType<<endl;
      LogTrace(metname)<<"hFWHM from langau fit: "<<hFWHM<<" for: "<<histname+muonType<<endl;
    }
   }

  if(histname=="ecalS9PointingMuDepositedEnergy_" && hPeak>0.2 && hPeak<0.3)
    energySummaryMap->setBinContent(binNumber,1, 1);
  if(histname=="ecalS9PointingMuDepositedEnergy_" && !(hPeak>0.2 && hPeak<0.3))
    energySummaryMap->setBinContent(binNumber,1, 0);
    
  if(histname=="hadS9PointingMuDepositedEnergy_" && hPeak>2 && hPeak<3.5)
    energySummaryMap->setBinContent(binNumber,2, 1);
  if(histname=="hadS9PointingMuDepositedEnergy_" && !(hPeak>2 && hPeak<3.5))
    energySummaryMap->setBinContent(binNumber,2, 0);

  if(histname=="hoS9PointingMuDepositedEnergy_" && hPeak>2 && hPeak<3.5)
    energySummaryMap->setBinContent(binNumber,3, 1);
  if(histname=="hoS9PointingMuDepositedEnergy_" && !(hPeak>2 && hPeak<3.5))
    energySummaryMap->setBinContent(binNumber,3, 0);

}


void MuonTestSummary::doMolteplicityTests(){

  MonitorElement* molteplicityHisto = dbe->get("Muons/MuonRecoAnalyzer/muReco");
  
  double molteplicity_GLB = double(molteplicityHisto->getBinContent(1)+molteplicityHisto->getBinContent(2))/double(molteplicityHisto->getEntries());
  LogTrace(metname)<<"molteplicity_GLB: "<<molteplicity_GLB<<endl;
  double molteplicity_TK = double(molteplicityHisto->getBinContent(3)+molteplicityHisto->getBinContent(4))/double(molteplicityHisto->getEntries());
  LogTrace(metname)<<"molteplicity_TK: "<<molteplicity_TK<<endl;
  double molteplicity_STA = double(molteplicityHisto->getBinContent(3)+molteplicityHisto->getBinContent(5))/double(molteplicityHisto->getEntries());
  LogTrace(metname)<<"molteplicity_STA: "<<molteplicity_STA<<endl;

  if(molteplicity_GLB<expMolteplicityGlb+0.06 && molteplicity_GLB>expMolteplicityGlb-0.06)
    molteplicitySummaryMap->setBinContent(1,1);
  else 
    molteplicitySummaryMap->setBinContent(1,0);

  if(molteplicity_TK<expMolteplicityTk+0.02 && molteplicity_TK>expMolteplicityTk-0.02)
    molteplicitySummaryMap->setBinContent(2,1);
  else 
    molteplicitySummaryMap->setBinContent(2,0);

  if(molteplicity_STA<expMolteplicitySta+0.1 && molteplicity_STA>expMolteplicitySta-0.1)
    molteplicitySummaryMap->setBinContent(3,1);
  else 
    molteplicitySummaryMap->setBinContent(3,0);

}
