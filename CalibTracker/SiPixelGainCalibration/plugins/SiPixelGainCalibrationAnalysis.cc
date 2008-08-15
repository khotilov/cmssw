// -*- C++ -*-
//
// Package:    SiPixelGainCalibrationAnalysis
// Class:      SiPixelGainCalibrationAnalysis
// 
/**\class SiPixelGainCalibrationAnalysis SiPixelGainCalibrationAnalysis.cc CalibTracker/SiPixelGainCalibrationAnalysis/src/SiPixelGainCalibrationAnalysis.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Freya Blekman
//         Created:  Wed Nov 14 15:02:06 CET 2007
// $Id: SiPixelGainCalibrationAnalysis.cc,v 1.25 2008/08/14 12:01:36 fblekman Exp $
//
//

// user include files
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "SiPixelGainCalibrationAnalysis.h"
#include <sstream>
#include <vector>
#include <math.h>
#include "TGraph.h"
#include "TMath.h"
//
// constructors and destructor
//
SiPixelGainCalibrationAnalysis::SiPixelGainCalibrationAnalysis(const edm::ParameterSet& iConfig):
  SiPixelOfflineCalibAnalysisBase(iConfig),
  conf_(iConfig),
  bookkeeper_(),
  bookkeeper_pixels_(),
  nfitparameters_(iConfig.getUntrackedParameter<int>("numberOfFitParameters",2)),
  fitfunction_(iConfig.getUntrackedParameter<std::string>("fitFunctionRootFormula","pol1")),
  reject_plateaupoints_(iConfig.getUntrackedParameter<bool>("suppressPlateauInFit",true)),
  reject_single_entries_(iConfig.getUntrackedParameter<bool>("suppressPointsWithOneEntryOrLess",true)),
  plateau_max_slope_(iConfig.getUntrackedParameter<double>("plateauSlopeMax",1.0)),
  reject_first_point_(iConfig.getUntrackedParameter<bool>("rejectVCalZero",true)),
  reject_badpoints_frac_(iConfig.getUntrackedParameter<double>("suppressZeroAndPlateausInFitFrac",0)),
  bookBIGCalibPayload_(iConfig.getUntrackedParameter<bool>("saveFullPayloads",false)),
  savePixelHists_(iConfig.getUntrackedParameter<bool>("savePixelLevelHists",false)),
  chi2Threshold_(iConfig.getUntrackedParameter<double>("minChi2NDFforHistSave",10)),
  chi2ProbThreshold_(iConfig.getUntrackedParameter<double>("minChi2ProbforHistSave",0.05)),
  maxGainInHist_(iConfig.getUntrackedParameter<double>("maxGainInHist",10)),
  maxChi2InHist_(iConfig.getUntrackedParameter<double>("maxChi2InHist",25)),
  saveALLHistograms_(iConfig.getUntrackedParameter<bool>("saveAllHistograms",false)),
  filldb_(iConfig.getUntrackedParameter<bool>("writeDatabase",false)),
  recordName_(conf_.getParameter<std::string>("record")),
  appendMode_(conf_.getUntrackedParameter<bool>("appendMode",true)),
  listofdetids_(conf_.getUntrackedParameter<std::vector<uint32_t> >("listOfDetIDs")),
  theGainCalibrationDbInput_(0),
  theGainCalibrationDbInputOffline_(0),
  theGainCalibrationDbInputHLT_(0),
  theGainCalibrationDbInputService_(iConfig),
  gainlow_(10.),gainhi_(0.),pedlow_(255.),pedhi_(0.)
{
  if(reject_single_entries_)
    min_nentries_=1;
  else
    min_nentries_=0;
  ::putenv("CORAL_AUTH_USER=me");
  ::putenv("CORAL_AUTH_PASSWORD=test");   
  edm::LogInfo("SiPixelGainCalibrationAnalysis") << "now using fit function " << fitfunction_ << ", which has " << nfitparameters_ << " free parameters. " << std::endl;
  func_= new TF1("func",fitfunction_.c_str());
}

SiPixelGainCalibrationAnalysis::~SiPixelGainCalibrationAnalysis()
{
}
// member functions
//
// ------------ method called once each job just before starting event loop  ------------

std::vector<float> SiPixelGainCalibrationAnalysis::CalculateAveragePerColumn(uint32_t detid, std::string label){
  std::vector<float> result;
  int ncols= bookkeeper_[detid][label]->getNbinsX();
  int nrows= bookkeeper_[detid][label]->getNbinsY();
  for(int icol=1; icol<=ncols; ++icol){
    float val=0;
    float ntimes =0;
    for(int irow=1; irow<=nrows; ++irow){
      val+= bookkeeper_[detid][label]->getBinContent(icol,irow);
      ntimes++;
    }
    val/= ntimes;
    result.push_back(val);
  }
  return result;
}

bool
SiPixelGainCalibrationAnalysis::checkCorrectCalibrationType()
{
  if(calibrationMode_=="GainCalibration")
    return true;
  else if(calibrationMode_=="unknown"){
    edm::LogInfo("SiPixelGainCalibrationAnalysis") <<  "calibration mode is: " << calibrationMode_ << ", continuing anyway..." ;
    return true;
  }
  else{
    //    edm::LogError("SiPixelGainCalibrationAnalysis") << "unknown calibration mode for Gain calibration, should be \"Gain\" and is \"" << calibrationMode_ << "\"";
  }
  return false;
}

void SiPixelGainCalibrationAnalysis::calibrationSetup(const edm::EventSetup&)
{
}
//------- summary printing method. Very verbose.
void
SiPixelGainCalibrationAnalysis::printSummary(){

  uint32_t detid=0;
  for(std::map<uint32_t,std::map<std::string,MonitorElement *> >::const_iterator idet = bookkeeper_.begin(); idet != bookkeeper_.end(); ++idet){
    if(detid==idet->first)
      continue;// only do things once per detid
    detid=idet->first;
    std::vector<float> gainvec=CalculateAveragePerColumn(detid,"gain_2d");
    std::vector<float> pedvec =CalculateAveragePerColumn(detid,"ped_2d");
    std::vector<float> chi2vec = CalculateAveragePerColumn(detid,"chi2_2d");
    std::ostringstream summarytext;

    summarytext << "Summary for det ID " << detid << "(" << translateDetIdToString(detid) << ")\n";
    summarytext << "\t Following: values per column: column #, gain, pedestal, chi2\n";
    for(uint32_t i=0; i<gainvec.size(); i++)
      summarytext << "\t " << i << " \t" << gainvec[i] << " \t" << pedvec[i] << " \t" << chi2vec[i] << "\n";
    summarytext << "\t list of pixels with high chi2 (chi2> " << chi2Threshold_ << "): \n";

    
    for(std::map<std::string, MonitorElement *>::const_iterator ipix = bookkeeper_pixels_[detid].begin(); ipix!=bookkeeper_pixels_[detid].end(); ++ipix)
      summarytext << "\t " << ipix->first << "\n";
    edm::LogInfo("SiPixelGainCalibrationAnalysis") << summarytext.str() << std::endl;

  }

}

// ------------ method called once each job just after ending the event loop  ------------

void 
SiPixelGainCalibrationAnalysis::calibrationEnd() {

  //  printSummary();
  
  // this is where we loop over all histograms and save the database objects
  if(filldb_)
    fillDatabase();
}
//-----------method to fill the database
void SiPixelGainCalibrationAnalysis::fillDatabase(){
 // only create when necessary.
  // process the minimum and maximum gain & ped values...
  edm::LogError("SiPixelGainCalibration::fillDatabase()") << "PLEASE do not fill the database directly from the gain calibration analyzer. This function is currently disabled and no DB payloads will be produced!" << std::endl;

}
// ------------ method called to do fits to all objects available  ------------
bool
SiPixelGainCalibrationAnalysis::doFits(uint32_t detid, std::vector<SiPixelCalibDigi>::const_iterator ipix)
{
  bool makehistopersistent = saveALLHistograms_;
  std::vector<uint32_t>::const_iterator detidfinder=find(listofdetids_.begin(),listofdetids_.end(),detid);
  if(detidfinder!=listofdetids_.end())
    makehistopersistent=true;
  // first, fill the input arrays to the TLinearFitter.
  double xvals[201];
  double yvals[200];
  double yerrvals[200];
  double xvalsall[201];
  float  xvalsasfloatsforDQM[201];
  double yvalsall[200];
  double yerrvalsall[200];
  int npoints=0;
  int nallpoints=0;
  bool use_point=true;
  for(uint32_t ii=0; ii< ipix->getnpoints() && ii<200; ii++){
    nallpoints++;
    use_point=true;
    xvalsasfloatsforDQM[ii]=xvalsall[ii]=vCalValues_[ii];
    yerrvalsall[ii]=yvalsall[ii]=0;
    if(ipix->getnentries(ii)>min_nentries_){
      yvalsall[ii]=ipix->getsum(ii)/(float)ipix->getnentries(ii);
      yerrvalsall[ii]=ipix->getsumsquares(ii)/(float)(ipix->getnentries(ii));
      yerrvalsall[ii]-=pow(yvalsall[ii],2);
      yerrvalsall[ii]=sqrt(yerrvalsall[ii])*sqrt(ipix->getnentries(ii));
    }
  }
  
  // calculate plateau value from last 3 full entries
  double plateauval=0;
  for(int ii=nallpoints-1; ii>=0 && npoints<10; --ii){
    if(yvalsall[ii]>0 && yerrvalsall[ii]>0){
      plateauval+=yvalsall[ii];
      npoints++;
    }
  }
  plateauval/=npoints;
  double maxgoodvalinfit=plateauval*(1.-reject_badpoints_frac_);
  if(maxgoodvalinfit<1)
    maxgoodvalinfit=255*(1.-reject_badpoints_frac_);
  npoints=0;
  for(int ii=0; ii<nallpoints; ++ii){
    // now selecting the appropriate points for the fit.
    use_point=true;
    if(reject_first_point_ && xvalsall[ii]<0.1)
      use_point=false;
    if(ipix->getnentries(ii)<=min_nentries_ && reject_single_entries_)
      use_point=false;
    if(ipix->getnentries(ii)==0 && reject_badpoints_)
      use_point=false;
    if(yvalsall[ii]>maxgoodvalinfit)
      use_point=false;
    if(ii>1 && fabs(yvalsall[ii]-yvalsall[ii-1])<1. && yvalsall[ii]>0.8*maxgoodvalinfit && reject_plateaupoints_)
      use_point=false;
    
    if(use_point){
      xvals[npoints]=xvalsall[ii];
      yvals[npoints]=yvalsall[ii];
      yerrvals[npoints]=yerrvalsall[ii];
      npoints++;
    }
  }
  int result=-1;
  float chi2,slope,intercept,prob;
  prob=chi2=-1;
  slope=intercept=0;
  TLinearFitter fitter(nfitparameters_,fitfunction_.c_str());
  //  fitter.SetFitOption("V");
  if(npoints<2)
    result=-2;
  else{
    fitter.AssignData(npoints,1,xvals,yvals,yerrvals);
    
    // and do the fit:
    result = fitter.Eval();
    
    if(result!=0)
      result=-3;
    else if(result==0)
      result=1;
    slope = fitter.GetParameter(1);
    intercept = fitter.GetParameter(0);

    chi2 = fitter.GetChisquare()/fitter.GetNumberFreeParameters();
    prob = TMath::Prob(fitter.GetChisquare(),fitter.GetNumberFreeParameters());
    
    if(isnan(slope) || isnan(intercept) ){
      //      std::cout << "fit result : "<< result << " "  << fitter.GetChisquare() << slope << " " << intercept << std::endl;
      makehistopersistent=true;
      //      for(int ibla=0; ibla<npoints; ibla++){
      //	std::cout << ibla << " " << xvals[ibla] << " " << yvals[ibla] << std::endl;
      //      }
      // and do the fit another way:
      TGraph gr(npoints,xvals,yvals);
      Int_t tempresult = gr.Fit("pol1","Q0");
      slope= gr.GetFunction("pol1")->GetParameter(1);
      intercept = gr.GetFunction("pol1")->GetParameter(0);
      if(tempresult==0)
	result=2;
      else
	result = -4;
      chi2=gr.GetFunction("pol1")->GetChisquare()/gr.GetFunction("pol1")->GetNumberFreeParameters();
      prob= TMath::Prob(gr.GetFunction("pol1")->GetChisquare(),gr.GetFunction("pol1")->GetNumberFreeParameters());
      //      std::cout << "modified fit result : "<< result << " "  << gr.GetFunction("pol1")->GetChisquare() << slope << " " << intercept << std::endl;
     
    }
    for(int i=0; i< func_->GetNpar();i++)
      func_->SetParameter(i,fitter.GetParameter(i));
    
    // convert the gain and pedestal parameters to functional form y= x/gain+ ped

    if(slope<0)
      makehistopersistent=true;
    if(chi2>chi2Threshold_ && chi2Threshold_>=0)
      makehistopersistent=true;
    if(prob<chi2ProbThreshold_)
      makehistopersistent=true;
    if(result<0)
      makehistopersistent=true;

    if(result>0){
      if(slope<gainlow_)
	gainlow_=slope;
      if(slope>gainhi_)
	gainhi_=slope;
      if(intercept>pedhi_)
	pedhi_=intercept;
      if(intercept<pedlow_)
	pedlow_=intercept;
      bookkeeper_[detid]["gain_1d"]->Fill(slope);
      bookkeeper_[detid]["gain_2d"]->setBinContent(ipix->col()+1,ipix->row()+1,slope);
      bookkeeper_[detid]["ped_1d"]->Fill(intercept);
      bookkeeper_[detid]["ped_2d"]->setBinContent(ipix->col()+1,ipix->row()+1,intercept);
      bookkeeper_[detid]["chi2_1d"]->Fill(chi2);
      bookkeeper_[detid]["chi2_2d"]->setBinContent(ipix->col()+1,ipix->row()+1,chi2);
      bookkeeper_[detid]["prob_1d"]->Fill(prob);
      bookkeeper_[detid]["prob_2d"]->setBinContent(ipix->col()+1,ipix->row()+1,prob);
    }
  }
  bookkeeper_[detid]["status_2d"]->setBinContent(ipix->col()+1,ipix->row()+1,result);
  
  if(!savePixelHists_)
    return true;
  if(makehistopersistent){
    setDQMDirectory(detid);
    std::ostringstream pixelinfo;
    pixelinfo << "GainCurve_row_" << ipix->row() << "_col_" << ipix->col();
    std::string tempname=translateDetIdToString(detid);
    tempname+="_";
    tempname+=pixelinfo.str();
    // and book the histo
    // fill the last value of the vcal array...
    xvalsasfloatsforDQM[nallpoints]=256;
    if(nallpoints>2)
      xvalsasfloatsforDQM[nallpoints]=2*xvalsasfloatsforDQM[nallpoints-1]-xvalsasfloatsforDQM[nallpoints-2];
    bookkeeper_pixels_[detid][pixelinfo.str()] =  bookDQMHistogram1D(detid,pixelinfo.str(),tempname,nallpoints,xvalsasfloatsforDQM);
    edm::LogInfo("SiPixelGainCalibrationAnalysis") << "now saving histogram for pixel " << tempname << ", gain = " << slope << ", pedestal = " << intercept << ", chi2/NDF=" << chi2 << "(prob:" << prob << "), fit status " << result;
    for(int ii=0; ii<nallpoints; ++ii){
      bookkeeper_pixels_[detid][pixelinfo.str()]->setBinContent(ii+1,yvalsall[ii]);
      bookkeeper_pixels_[detid][pixelinfo.str()]->setBinError(ii+1,yerrvalsall[ii]);
    }
    
    addTF1ToDQMMonitoringElement(bookkeeper_pixels_[detid][pixelinfo.str()],func_);
  } 
  return true;
}
// ------------ method called to do fill new detids  ------------
void 
SiPixelGainCalibrationAnalysis::newDetID(uint32_t detid)
{
  setDQMDirectory(detid);
  std::string tempname=translateDetIdToString(detid);
  bookkeeper_[detid]["gain_1d"] = bookDQMHistogram1D(detid,"Gain1d","gain for "+tempname,100,0.,maxGainInHist_);
  bookkeeper_[detid]["gain_2d"] = bookDQMHistoPlaquetteSummary2D(detid, "Gain2d","gain for "+tempname);
  bookkeeper_[detid]["ped_1d"] = bookDQMHistogram1D(detid,"Pedestal1d","pedestal for "+tempname,256,0.,256.);
  bookkeeper_[detid]["ped_2d"] = bookDQMHistoPlaquetteSummary2D(detid,"Pedestal2d","pedestal for "+tempname);
  bookkeeper_[detid]["chi2_1d"] = bookDQMHistogram1D(detid,"GainChi2NDF1d","#chi^{2}/NDOF for "+tempname,100,0.,maxChi2InHist_);
  bookkeeper_[detid]["chi2_2d"] = bookDQMHistoPlaquetteSummary2D(detid,"GainChi2NDF2d","#chi^{2}/NDOF for "+tempname);
  bookkeeper_[detid]["prob_1d"] = bookDQMHistogram1D(detid,"GainChi2Prob1d","P(#chi^{2},NDOF) for "+tempname,100,0.,1.0);
  bookkeeper_[detid]["prob_2d"] = bookDQMHistoPlaquetteSummary2D(detid,"GainChi2Prob2d","P(#chi^{2},NDOF) for "+tempname);
  bookkeeper_[detid]["status_2d"] = bookDQMHistoPlaquetteSummary2D(detid,"GainFitResult2d","Fit result for "+tempname);

}
//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelGainCalibrationAnalysis);
