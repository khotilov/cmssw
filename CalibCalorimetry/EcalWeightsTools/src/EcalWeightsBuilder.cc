/**
 * \file EcalWeightsBuilder.cc 
 * This module is used to produce the weights needed by
 * the reconstruction of amplitude (RecHits) in the ECAL.
 * 
 * $Date: 2006/07/27 17:42:33 $
 * $Revision: 1.3 $
 * Author Alexandre Zabi
 *
 * Modif-Alex Zabi 25/07/2006
 * This module has been modified to generate an extra sets of weights
 * to be used for the reconstruction at high energy
 * At high energy 2 configurations can happen:
 * 1st: samples 5, 6, 7, 8, 9 are in gain 6
 * 2nd: samples 4, 5, 6, 7, 8 are in gain 6
 * An extra set of 25 weights is produced. The number of TdcBin 
 * is increased to 50 but the parameter nTdcBin must remain 25.
 *
*/

#include "TFile.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "CalibCalorimetry/EcalWeightsTools/interface/ComputeWeights.h"
#include "CalibCalorimetry/EcalWeightsTools/interface/EcalWeightsBuilder.h"

//include reference signal representation
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h"
#include "SimCalorimetry/EcalSimAlgos/interface/EcalShape.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

EcalWeightsBuilder::EcalWeightsBuilder(edm::ParameterSet const& pSet)
{
  //verbosity
  verbosity_        = pSet.getUntrackedParameter("verbosity", 1U);

  //Prints out the individual weights files: ampWeights.txt, pedWeights.txt etc..
  debug_            = pSet.getUntrackedParameter<bool>("debug", false);

  //GENERATION OF GROUPID WEIGHTS:
  // This generates weights for Test Beam analysis. The Module EcalGroupIdBuilder
  // must have been called priored to EcalWeightsBuilder in order to assign each
  // crystal to a certain group Id number.
  // When generating group_id weights, the supermodule number must be specified
  // and the corresponding group id file produced by EcalGroupIdBuilder must be
  // accessible: GroupId_SMXX.out 
  gen_groupID_      = pSet.getUntrackedParameter<bool>("gen_groupID", false);
  SuperModule_      = pSet.getUntrackedParameter<unsigned int>("SuperModule", 0);

  //REFERENCE TIMING: This is used to generated group id weights. This measures
  // the difference between the reference and mean timing of SuperModule.
  // tMaxRef_ = 5.5 if default setting for the Supermodule
  tMaxRef_          = pSet.getUntrackedParameter<double>("tMaxRef", 5.5); 
  BinOfMax_         = int(tMaxRef_);

  //GENERATING DEFAULT WEIGHTS:
  // This generates the default sets of weights for CMS and Test Beam (simulation)
  // NOTE: if this option is chosen all, other options refering to groupIds are
  // discraded
  default_weights_  = pSet.getUntrackedParameter<bool>("default_weights", true);
  if(default_weights_) 
    {tMaxRef_ = 5.5;
      if(gen_groupID_) {
	gen_groupID_ = false; 
	if(verbosity_ > 1) std::cout << "EcalWeightsBuilder: WARNING: " 
				     << " the group id weights will not be generated"
				     << " default_weights option should be false" << std::endl;
      }//group id check
    }//default weight`s

  //SAMPLES: number of samples in the electronics readout
  nSamples_          = pSet.getUntrackedParameter<unsigned int>("nSamples", 10);

  //TDCBINS: number of Tdc bins
  nTdcBins_          = pSet.getUntrackedParameter<unsigned int>("nTdcBins", 25);

  //WEIGHTS OPTIONS:
  // WEIGHTS OPTIONS BEFORE GAIN SWITCHING
  // number of signal sample and pre-sample used for amplitude reconstruction
  nPulseSamples_     = pSet.getUntrackedParameter<unsigned int>("nPulseSamples", 5);
  nPrePulseSamples_  = pSet.getUntrackedParameter<unsigned int>("nPrePulseSamples", 3);
  //pedestal subtracting weights
  doFitBaseline_     = pSet.getUntrackedParameter<bool>("doFitBaseline", true);
  //jitter compensating weights
  doFitTime_         = pSet.getUntrackedParameter<bool>("doFitTime", true);

  // WEIGHTS OPTIONS AFTER GAIN SWITCHING
  // number of signal sample and pre-sample used for amplitude reconstruction after gain switching
  nPulseSamples_gain_    = pSet.getUntrackedParameter<unsigned int>("nPulseSamples_gain", 5);
  nPrePulseSamples_gain_ = pSet.getUntrackedParameter<unsigned int>("nPrePulseSamples_gain", 3);
  doFitBaseline_gain_    = pSet.getUntrackedParameter<bool>("doFitBaseline_gain", false);
  doFitTime_gain_        = pSet.getUntrackedParameter<bool>("doFitTime_gain", true);

  if(default_weights_)
    {
      if(verbosity_ > 1) std::cout << "EcalWeightsBuilder: WARNING: " 
				   << " PRODUCING DEFAULT WEIGHTS FOR CMS and TEST BEAM SIMULATION" 
				   << std::endl;
      
      //If default weights are produced (for data taking)
      //Jitter compensating weights are not necessary
      nPulseSamples_         = 5;
      nPrePulseSamples_      = 3;
      doFitBaseline_         = true;
      doFitTime_             = false;
      nPulseSamples_gain_    = 5;
      nPrePulseSamples_gain_ = 3;
      doFitBaseline_gain_    = false;
      doFitTime_gain_        = false;

      //default number of Tdcbins
      nTdcBins_              = 25;
    }//default weights

}//CONSTRUCTOR

EcalWeightsBuilder::~EcalWeightsBuilder()
{

}//DESTRUCTOR

void EcalWeightsBuilder::analyze(const edm::Event& evt, const edm::EventSetup& evtSetup) 
{
  std::cout << "STARTING ECAL WEIGHTS PRODUCTION"       << std::endl; 
  std::cout << "REFERENCE TIME OF MAXIMUM=" << tMaxRef_ << std::endl;

  //WEIGHTS FILES
  std::ofstream WeightsFileCMS("WeightsFileCMS.out"); //weights file for CMS
  std::ofstream WeightsFileTB("WeightsFileTB.out");   //weights file for Test Beam

  //LOADING REFERENCE SIGNAL REPRESENTATION
  EcalSimParameterMap parameterMap;
  
  EBDetId   barrel(1,1);
  double    thisPhase = parameterMap.simParameters(barrel).timePhase();
  EcalShape theShape(thisPhase);

  std::cout << "Parameters for the ECAL MGPA shape \n"     << std::endl;
  std::cout << "Rising time for ECAL shape (timePhase) = " << parameterMap.simParameters(barrel).timePhase()    << std::endl;
  std::cout << "Bin of maximum = "                         << parameterMap.simParameters(barrel).binOfMaximum() << std::endl;

  double ToM        = theShape.computeTimeOfMaximum();
  double T0         = theShape.computeT0();
  double risingTime = theShape.computeRisingTime();

  std::cout << "\n Maximum time from tabulated values = " << ToM        << std::endl;
  std::cout << "\n Tzero from tabulated values        = " << T0         << std::endl;
  std::cout << "\n Rising time from tabulated values  = " << risingTime << std::endl;
  std::cout << std::endl;
  
  //time of maximum
  double tMax      = ToM/25.0;
  double tMax_gain = tMax + 1;    
  std::cout << "TIME OF MAXIMUM =" << tMax << std::endl;

  //DETERMINING THE GROUP-ID ACCORDING TO THE SIGNAL TIMING ///////////////////////////////////////////////////////
  int    nGroupId   = 1;
  double meanTiming = 5.5;

  if(gen_groupID_){

    //looking for the corresponding group id file
    char* grpid_file = new char[20];
    std::sprintf (grpid_file,"GroupId_SM%02u.out",SuperModule_);
    std::cout << "LOOKING FOR GROUPID FILE=" << grpid_file << std::endl;

    std::ifstream groupid_in(grpid_file);
    if (groupid_in.is_open())
      {
	groupid_in >> nGroupId >> meanTiming;
	std::cout << std::endl;
	std::cout << "GENERATING WEIGHTS FOR GROUPIDs" << std::endl;
	std::cout << "NGroup = " << nGroupId << " Mean Timing = " << meanTiming << std::endl;
	groupid_in.close();
	//TMax is set to the mean timing 
	tMax = meanTiming;
      }//generating groupIds
    else
      {
	std::cout << "EcalWeightsBuilder : Error: "
		  << "can not open group id file file " << grpid_file
		  << std::endl;
	abort();
      }
    delete [] grpid_file;
  }// generating group id

  //mean timing for each group
  double meanTimingGroup[3]; //only 3 group Id for the moment.
  meanTimingGroup[0] = meanTiming - 0.08; 
  meanTimingGroup[1] = meanTiming;
  meanTimingGroup[2] = meanTiming + 0.08;
  if(gen_groupID_){for(int i=0; i < nGroupId ; ++i) std::cout << meanTimingGroup[i] << std::endl;
    std::cout << std::endl;}

  //
  //COMPUTING WEIGHTS ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //CREATE OBJECT TO COMPUTE GAIN 12 (low energy) WEIGHTS
  // For both the simulation and CMS: the amplitude weights are not jitter compensating
  // the timing weights are determined with an extra ComputeWeights object
  ComputeWeights weights(verbosity_, doFitBaseline_, doFitTime_, nPulseSamples_, nPrePulseSamples_);
  ComputeWeights weights_time(verbosity_, doFitBaseline_, true, nPulseSamples_, nPrePulseSamples_);

  // Create output weights files 
  std::ofstream fAmpWeights_TB, fPedWeights_TB, fTimeWeights_TB, fChi2_TB;
  std::ofstream fAmpWeights_CMS, fPedWeights_CMS, fTimeWeights_CMS, fChi2_CMS;
  if(debug_){
    //Amplitude weights:
    fAmpWeights_TB.open("ampWeights_TB.txt");
    //Pedestal weights:
    if (doFitBaseline_) fPedWeights_TB.open("pedWeights_TB.txt");
    //Timing weights:
    if (default_weights_ || doFitTime_) fTimeWeights_TB.open("timeWeights_TB.txt");
    //chi2 matrix
    fChi2_TB.open("chi2Matrix_TB.txt");

    //FOR CMS
    //Amplitude weights:
    fAmpWeights_CMS.open("ampWeights_CMS.txt");
    //Pedestal weights:
    if (doFitBaseline_) fPedWeights_CMS.open("pedWeights_CMS.txt");
    //Timing weights:
    if (default_weights_ || doFitTime_) fTimeWeights_CMS.open("timeWeights_CMS.txt");
    //chi2 matrix
    fChi2_CMS.open("chi2Matrix_CMS.txt");
  }//debug files

  // CREATE OBJECT TO COMPUTE GAIN 6 and GAIN 1 (high energy) WEIGHTS
  ComputeWeights weights_gain(verbosity_, doFitBaseline_gain_, doFitTime_gain_, nPulseSamples_gain_, nPrePulseSamples_gain_);
  ComputeWeights weights_time_gain(verbosity_, doFitBaseline_gain_, true, nPulseSamples_gain_, nPrePulseSamples_gain_);

//   // CREATE OBJECT TO COMPUTE GAIN 6 and GAIN 1 (high energy) EXTRA WEIGHTS
//   //For resonctruction at high energy
//   ComputeWeights weights_gain_extra(verbosity_, doFitBaseline_gain_, doFitTime_gain_, nPulseSamples_gain_, nPrePulseSamples_gain_);
//   ComputeWeights weights_time_gain_extra(verbosity_, doFitBaseline_gain_, true, nPulseSamples_gain_, nPrePulseSamples_gain_);

  // Create output files
  std::ofstream fAmpWeightsGain_TB, fPedWeightsGain_TB, fTimeWeightsGain_TB, fChi2Gain_TB; 
  std::ofstream fAmpWeightsGain_CMS, fPedWeightsGain_CMS, fTimeWeightsGain_CMS, fChi2Gain_CMS; 
  if(debug_){
    //Amplitude weights:
    fAmpWeightsGain_TB.open("ampWeightsAfterGainSwitch_TB.txt");
    //Pedestal weights:
    if (doFitBaseline_gain_) fPedWeightsGain_TB.open("pedWeightsAfterGainSwitch_TB.txt");
    //Timing weights:
    if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_TB.open("timeWeightsAfterGainSwitch_TB.txt"); 
    //chi2 matrix
    fChi2Gain_TB.open("chi2MatrixAfterGainSwitch_TB.txt");

    //FOR CMS
    //Amplitude weights:
    fAmpWeightsGain_CMS.open("ampWeightsAfterGainSwitch_CMS.txt");
    //Pedestal weights:
    if (doFitBaseline_gain_) fPedWeightsGain_CMS.open("pedWeightsAfterGainSwitch_CMS.txt");
    //Timing weights:
    if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_CMS.open("timeWeightsAfterGainSwitch_CMS.txt"); 
    //chi2 matrix
    fChi2Gain_CMS.open("chi2MatrixAfterGainSwitch_CMS.txt");
  }//debug file

  //LOOPING OVER GROUPIDs /////////////////////////////////
  for(int igroupId = 0; igroupId < nGroupId; ++igroupId){
    if(debug_) std::cout << "%%%%%%%%%%%%%%%%%%%%%%%GROUPID " << igroupId << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    std::vector<double> pulseShape(nSamples_);
    std::vector<double> pulseShapeDerivative(nSamples_);

    //vector to save the weights for high energy
    std::vector<double> extraWeights_amp; 
    std::vector<double> extraWeights_ped; 
    std::vector<double> extraWeights_time;
    std::vector<double> extraWeights_chi2;
 
    //Shift the shape (in ns)
    // > 0  -> shift left  (earlier)
    // < 0  -> shift right (later)
    // 5.5 is the default time of maximum. The default shape from EcalShape has 
    // its peak set in the following weights producer right on 5.5 = 137.5 ns. 

    double shiftTime = 0.0; //;2.0; //shift to Tmax of Xtal 704
    if(gen_groupID_){
      double diff_time = (tMaxRef_ - meanTimingGroup[igroupId])*25.0;
      shiftTime = (diff_time-(int)diff_time) <= 0.5 ? (int)diff_time : ((int)diff_time + 1); 
      std::cout << "Shifting shape by " << shiftTime << " ns" << std::endl;
    }// generating groupIds

    for (int unsigned iTdcBin = 0; iTdcBin < nTdcBins_; iTdcBin++) {
      std::vector<double> sample_check;//for debugging purposes

      //Determination of tdcbin
      double tdcbin = -0.5 + (iTdcBin * 0.04) + 0.02; //return value in bin center

      //Determination of tzero
      //  double tzero = risingTime-(parameterMap.simParameters(barrel).binOfMaximum()-tdcbin)*25.;
      //  risingTime-(parameterMap.simParameters(barrel).binOfMaximum()-(tdcbin+1))*25.0
      double tzero = risingTime -(BinOfMax_-tdcbin)*25.;

      if(debug_) std::cout << "BIN=" << iTdcBin << " TDCbin=" << tdcbin << " TZERO=" << tzero << std::endl;
    
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {

	if(debug_)
	  std::cout << iTdcBin << " " << tzero + iSample*nTdcBins_ << " Value=" << (theShape)(tzero + shiftTime + iSample*nTdcBins_) 
		    << " DERVIATIVE=" << theShape.derivative(tzero + shiftTime + iSample*nTdcBins_) << std::endl;
	
	pulseShape[iSample]           = (theShape)(tzero + shiftTime + iSample*nTdcBins_);
	pulseShapeDerivative[iSample] = theShape.derivative(tzero + shiftTime + iSample*nTdcBins_);
	
	//saving sampling for sum checks
	sample_check.push_back(pulseShape[iSample]); 
      }//loop sample
      
      //COMPUTING WEIGHTS FOR LOW ENERGY
      if (!weights.compute(pulseShape, pulseShapeDerivative, tMax)) {
	std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		  << " Impossible to compute weights iTdcBin = " << iTdcBin
		  << std::endl;
	continue;
      }//compute weights.
      if (!weights_time.compute(pulseShape, pulseShapeDerivative, tMax)) {
	std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		  << " Impossible to compute time weights iTdcBin = " << iTdcBin
		  << std::endl;
	continue;
      }//compute time weights.
      
      if(debug_){
	//CHECKING SUMS
	double sum_weights = 0.0;
	double sum_wf      = 0.0;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	  sum_weights += weights.getAmpWeight(iSample);
	  sum_wf      += weights.getAmpWeight(iSample) * sample_check[iSample];	
	}//loop sample
	std::cout << "SUM WEIGHTS=" << sum_weights << " " << " SUM Wi*Fi=" <<  sum_wf << std::endl;
	
	//WRITING WEIGHTS INTO OUTPUT FILE
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	  fAmpWeights_TB << std::setw(10) << std::setprecision(7) << weights.getAmpWeight(iSample) << " "; 
	  if (doFitBaseline_)
	    fPedWeights_TB << std::setw(10) << std::setprecision(7) << weights.getPedWeight(iSample) << " ";      
	  if(default_weights_) fTimeWeights_TB << std::setw(10) << std::setprecision(7) << weights_time.getTimeWeight(iSample) << " "; 
	  else if (doFitTime_) { fTimeWeights_TB << std::setw(10) << std::setprecision(7) << weights.getTimeWeight(iSample) << " ";}      
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    fChi2_TB << std::setw(10) << std::setprecision(7) << weights.getChi2Matrix(iSample, iSample2) << " "; 
	  fChi2_TB << std::endl;
	}//loop sample
	fAmpWeights_TB << std::endl;
	if (doFitBaseline_) fPedWeights_TB << std::endl;
	if (default_weights_ || doFitTime_) fTimeWeights_TB << std::endl;  

	if (iTdcBin == 12 && !gen_groupID_){
	  for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	    fAmpWeights_CMS << std::setw(10) << std::setprecision(7) << weights.getAmpWeight(iSample) << " "; 
	    if (doFitBaseline_)
	      fPedWeights_CMS << std::setw(10) << std::setprecision(7) << weights.getPedWeight(iSample) << " ";      
	    if(default_weights_) fTimeWeights_CMS << std::setw(10) << std::setprecision(7) << weights_time.getTimeWeight(iSample) << " "; 
	    else if (doFitTime_) { fTimeWeights_CMS << std::setw(10) << std::setprecision(7) << weights.getTimeWeight(iSample) << " ";}      
	    for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	      fChi2_CMS << std::setw(10) << std::setprecision(7) << weights.getChi2Matrix(iSample, iSample2) << " "; 
	    fChi2_CMS << std::endl;
	  }//loop sample
	  fAmpWeights_CMS << std::endl;
	  if (doFitBaseline_) fPedWeights_CMS << std::endl;
	  if (default_weights_ || doFitTime_) fTimeWeights_CMS << std::endl;  
	}//CMS
      }//debug files
      
      //COMPUTING WEIGHTS FOR HIGH ENERGY
      //1st CONFIGURATION SAMPLES : samples 5 ,6 ,7 ,8 ,9 are in gain 6
      //In that case Tmax = Tmax + 1;
      if (!weights_gain.compute(pulseShape, pulseShapeDerivative, tMax_gain)) { 
	if (verbosity_)
	  std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		    << " Impossible to compute weights iTdcBin =" << iTdcBin
		    << std::endl;
	continue;
      }//compute weights.
      if (!weights_time_gain.compute(pulseShape, pulseShapeDerivative, tMax_gain)) { 
	if (verbosity_)
	  std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		    << " Impossible to compute time weights iTdcBin =" << iTdcBin
		    << std::endl;
	continue;
      }//compute weights.
    
      if(debug_){
	//WRITING HIGH GAIN WEIGHTS INTO OUTPUT FILE
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	  fAmpWeightsGain_TB << std::setw(10) << std::setprecision(7) << weights_gain.getAmpWeight(iSample) << " "; 
	  if (doFitBaseline_gain_)
	    fPedWeightsGain_TB << std::setw(10) << std::setprecision(7) << weights_gain.getPedWeight(iSample) << " "; 
	  if(default_weights_) fTimeWeightsGain_TB << std::setw(10) << std::setprecision(7) << weights_time_gain.getTimeWeight(iSample) << " "; 
	  else if (doFitTime_gain_) { fTimeWeightsGain_TB << std::setw(10) << std::setprecision(7) << weights_gain.getTimeWeight(iSample) << " "; }
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    fChi2Gain_TB << std::setw(10) << std::setprecision(7) << weights_gain.getChi2Matrix(iSample, iSample2) << " "; 
	  fChi2Gain_TB << std::endl;
	}
	fAmpWeightsGain_TB << std::endl;
	if (doFitBaseline_gain_) fPedWeightsGain_TB << std::endl;
	if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_TB << std::endl;

	if (iTdcBin == 12 && !gen_groupID_){
	  if(default_weights_) {
	    for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	      //if(iSample == unsigned(parameterMap.simParameters(barrel).binOfMaximum()-1)) 
	      if(iSample == BinOfMax_)
		fAmpWeightsGain_CMS << std::setw(10) << std::setprecision(7) << 1.0 << " ";
	      else fAmpWeightsGain_CMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	    }//loop sample
	    fAmpWeightsGain_CMS << std::endl;
	    for (int unsigned iSample = 0; iSample < nSamples_; iSample++) fPedWeightsGain_CMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	    fPedWeightsGain_CMS << std::endl;
	    for (int unsigned iSample = 0; iSample < nSamples_; iSample++) fTimeWeightsGain_CMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	    fTimeWeightsGain_CMS << std::endl;
	    for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	      for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
		fChi2Gain_CMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	      fChi2Gain_CMS << std::endl;
	    }//loop sample
	  }//default weights
	  else
	    {
	      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
		fAmpWeightsGain_CMS << std::setw(10) << std::setprecision(7) << weights_gain.getAmpWeight(iSample) << " "; 
		if (doFitBaseline_gain_)
		  fPedWeightsGain_CMS  << std::setw(10) << std::setprecision(7) << weights_gain.getPedWeight(iSample) << " "; 
		if (doFitTime_gain_)
		  fTimeWeightsGain_CMS << std::setw(10) << std::setprecision(7) << weights_gain.getTimeWeight(iSample) << " "; 
		for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
		  fChi2Gain_CMS << std::setw(10) << std::setprecision(7) << weights_gain.getChi2Matrix(iSample, iSample2) << " "; 
		fChi2Gain_CMS << std::endl;
	      }//loop sample
	      fAmpWeightsGain_CMS << std::endl;
	      if (doFitBaseline_gain_) fPedWeightsGain_CMS  << std::endl;
	      if (doFitTime_gain_)     fTimeWeightsGain_CMS << std::endl;
	    }
	}//CMS
      }//debug files
      
      //CREATING THE WEIGHTS FILE:
      if (!iTdcBin)
	WeightsFileTB << igroupId << " " << nSamples_ << " " << nTdcBins_*2 << std::endl;
      
      //low energy
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << weights.getAmpWeight(iSample) << " "; }
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	if (doFitBaseline_) {WeightsFileTB << std::setw(10) << std::setprecision(7) << weights.getPedWeight(iSample) << " ";}
	else               {WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";}}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	if(default_weights_) WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_time.getTimeWeight(iSample) << " "; 
	else { if (doFitTime_) WeightsFileTB << std::setw(10) << std::setprecision(7) << weights.getTimeWeight(iSample) << " ";
	  else WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " "; }}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	  WeightsFileTB << std::setw(10) << std::setprecision(7) << weights.getChi2Matrix(iSample, iSample2) << " "; 
	WeightsFileTB << std::endl;
      }//loop sample
      
      //high energy
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_gain.getAmpWeight(iSample) << " "; }
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	if (doFitBaseline_gain_) {WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_gain.getPedWeight(iSample) << " "; }
	else                    {WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";}}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	if(default_weights_)      WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_time_gain.getTimeWeight(iSample) << " "; 
	else { if (doFitTime_gain_) WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_gain.getTimeWeight(iSample) << " ";
	  else WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";}}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	  WeightsFileTB << std::setw(10) << std::setprecision(7) << weights_gain.getChi2Matrix(iSample, iSample2) << " "; 
	WeightsFileTB << std::endl;
      }//loop sample      
      
      // WEIGHTS FOR CMS /////////////////////////////////////////////////////////////////////////////
      if (iTdcBin == 12 && default_weights_ && !gen_groupID_){
	WeightsFileCMS << 0 << " " << nSamples_ << " " << 1 << std::endl;
	
	//low energy
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	  WeightsFileCMS << std::setw(10) << std::setprecision(7) << weights.getAmpWeight(iSample) << " ";
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	  WeightsFileCMS << std::setw(10) << std::setprecision(7) << weights.getPedWeight(iSample) << " ";
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++)
	  WeightsFileCMS << std::setw(10) << std::setprecision(7) << weights_time.getTimeWeight(iSample) << " "; 
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    WeightsFileCMS << std::setw(10) << std::setprecision(7) << weights.getChi2Matrix(iSample, iSample2) << " "; 
	  WeightsFileCMS << std::endl;
	}//loop sample
	
	//high energy
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	  //if(iSample == unsigned(parameterMap.simParameters(barrel).binOfMaximum()-1)) 
	  if(iSample == BinOfMax_) 
	    WeightsFileCMS << std::setw(10) << std::setprecision(7) << 1.0 << " ";
	  else WeightsFileCMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	}
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) WeightsFileCMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) WeightsFileCMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	WeightsFileCMS << std::endl;
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    WeightsFileCMS << std::setw(10) << std::setprecision(7) << 0.0 << " ";
	  WeightsFileCMS << std::endl;
	}//loop sample
      }//CMS WEIGHTS
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////  
      
      //EXTRA SETS OF WEIGHTS FOR RECONSTRUCTION AT HIGH ENERGY
      //2nd CONFIGURATION SAMPLES : samples 4, 5 ,6 ,7 ,8 are in gain 6
      //In that case Tmax = Tmax;
      if (!weights_gain.compute(pulseShape, pulseShapeDerivative, tMax)) { 
	if (verbosity_)
	  std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		    << " Impossible to compute weights iTdcBin =" << iTdcBin
		    << std::endl;
	continue;
      }//compute weights.
      if (!weights_time_gain.compute(pulseShape, pulseShapeDerivative, tMax)) { 
	if (verbosity_)
	  std::cout << "EcalWeightsBuilder::~EcalWeightsBuilder: Warning:"
		    << " Impossible to compute time weights iTdcBin =" << iTdcBin
		    << std::endl;
	continue;
      }//compute weights.

      //saving these extra weights into a vector
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++)   extraWeights_amp.push_back(weights_gain.getAmpWeight(iSample));
      if (doFitBaseline_gain_) 
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) extraWeights_ped.push_back(weights_gain.getPedWeight(iSample));
      if(default_weights_) 
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) extraWeights_time.push_back(weights_time_gain.getTimeWeight(iSample));
      else if (doFitTime_gain_) 
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) extraWeights_time.push_back(weights_gain.getTimeWeight(iSample));
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	  extraWeights_chi2.push_back(weights_gain.getChi2Matrix(iSample, iSample2));
      }//loop sample

    }//loop TDC 1st

    //WRITING EXTRA WEIGHTS
    // for the reconstruction a low energy, extra weights are all 0 
    // only at high energy weights are used.
    for (int unsigned iTdcBin = 0; iTdcBin < nTdcBins_; iTdcBin++) {
      if(debug_){
	//WRITING ZEROS FOR LOW ENERGY
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	  fAmpWeights_TB << std::setw(10) << std::setprecision(7) << 0.0 << " "; 
	  if (doFitBaseline_)  fPedWeights_TB  << std::setw(10) << std::setprecision(7) << 0.0 << " ";      
	  if(default_weights_ || doFitTime_) fTimeWeights_TB << std::setw(10) << std::setprecision(7) << 0.0 << " "; 
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    fChi2_TB << std::setw(10) << std::setprecision(7) << 0.0 << " "; 
	  fChi2_TB << std::endl;
	}//loop sample
	fAmpWeights_TB << std::endl;
	if (doFitBaseline_) fPedWeights_TB << std::endl;
	if (default_weights_ || doFitTime_) fTimeWeights_TB << std::endl;  
	
	//WRITING HIGH GAIN WEIGHTS INTO OUTPUT FILE
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++)
	  fAmpWeightsGain_TB << std::setw(10) << std::setprecision(7)   << extraWeights_amp[iSample+iTdcBin*nSamples_] << " "; 	
	if (doFitBaseline_gain_) {
	  for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	    fPedWeightsGain_TB << std::setw(10) << std::setprecision(7) << extraWeights_ped[iSample+iTdcBin*nSamples_] << " "; 
	}//baseline	
	if(default_weights_ || doFitTime_gain_) 
	  for (int unsigned iSample = 0; iSample < nSamples_; iSample++)
	    fTimeWeightsGain_TB << std::setw(10) << std::setprecision(7) << extraWeights_time[iSample+iTdcBin*nSamples_] << " "; 	
	for (int unsigned iSample = 0; iSample < nSamples_; iSample++){
	  for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	    fChi2Gain_TB << std::setw(10) << std::setprecision(7) << extraWeights_chi2[iSample2+iSample*nSamples_+iTdcBin*nSamples_*nSamples_] << " "; 
	  fChi2Gain_TB << std::endl;
	}//loop sample	
	fAmpWeightsGain_TB << std::endl;
	if (doFitBaseline_gain_) fPedWeightsGain_TB << std::endl;
	if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_TB << std::endl;	
      }//debug

     //low energy
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";	
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " "; 
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	  WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " "; 
	WeightsFileTB << std::endl;
      }//loop sample
      
      //high energy
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	WeightsFileTB << std::setw(10) << std::setprecision(7) << extraWeights_amp[iSample+iTdcBin*nSamples_] << " "; }
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	if (doFitBaseline_gain_) {WeightsFileTB << std::setw(10) << std::setprecision(7) << extraWeights_ped[iSample+iTdcBin*nSamples_] << " "; }
	else                     {WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";}}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) {
	if(default_weights_ || doFitTime_gain_) WeightsFileTB << std::setw(10) << std::setprecision(7) << extraWeights_time[iSample+iTdcBin*nSamples_] << " "; 
	else WeightsFileTB << std::setw(10) << std::setprecision(7) << 0.0 << " ";}
      WeightsFileTB << std::endl;
      for (int unsigned iSample = 0; iSample < nSamples_; iSample++) { 
	for (int unsigned iSample2 = 0; iSample2 < nSamples_; iSample2++)
	  WeightsFileTB << std::setw(10) << std::setprecision(7) << extraWeights_chi2[iSample2+iSample*nSamples_+iTdcBin*nSamples_*nSamples_] << " "; 
	WeightsFileTB << std::endl;
      }//loop sample
    }//loop TDC 2nd

  }//loop group IDs

  //CLOSING FILES
  if(debug_){
    fAmpWeights_TB.close();
    if (doFitBaseline_) fPedWeights_TB.close();
    if (default_weights_ || doFitTime_) fTimeWeights_TB.close();
    fChi2_TB.close();
    
    fAmpWeightsGain_TB.close();
    if (doFitBaseline_gain_) fPedWeightsGain_TB.close();
    if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_TB.close();
    fChi2Gain_TB.close();

    fAmpWeights_CMS.close();
    if (doFitBaseline_) fPedWeights_CMS.close();
    if (default_weights_ || doFitTime_) fTimeWeights_CMS.close();
    fChi2_CMS.close();
    
    fAmpWeightsGain_CMS.close();
    if (doFitBaseline_gain_) fPedWeightsGain_CMS.close();
    if (default_weights_ || doFitTime_gain_) fTimeWeightsGain_CMS.close();
    fChi2Gain_CMS.close();
  }//debug files

  WeightsFileCMS.close();
  WeightsFileTB.close();

}//ANALYZE
