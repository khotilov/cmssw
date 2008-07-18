#include "DQM/SiStripCommissioningAnalysis/interface/SamplingAlgorithm.h"
#include "CondFormats/SiStripObjects/interface/SamplingAnalysis.h"
#include "DataFormats/SiStripCommon/interface/SiStripHistoTitle.h"
#include "DataFormats/SiStripCommon/interface/SiStripEnumsAndStrings.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQM/SiStripCommissioningAnalysis/interface/SiStripPulseShape.h"
#include "TProfile.h"
#include "TF1.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace sistrip;

// ----------------------------------------------------------------------------
// 
SamplingAlgorithm::SamplingAlgorithm( SamplingAnalysis* const anal ) 
  : CommissioningAlgorithm(anal),
    histo_(0,""),
    deconv_fitter_(0),
    peak_fitter_(0),
    samp_(0)
{
   peak_fitter_ = new TF1("peak_fitter",fpeak_convoluted,-4800,0,5);
   peak_fitter_->SetNpx(2000);
   peak_fitter_->FixParameter(0,0);
   peak_fitter_->SetParLimits(1,0,4800);
   peak_fitter_->SetParLimits(2,0,20);
   peak_fitter_->FixParameter(3,50);
   peak_fitter_->SetParLimits(4,0,50);
   peak_fitter_->SetParameters(0.,1250,10,50,10);
   deconv_fitter_ = new TF1("deconv_fitter",fdeconv_convoluted,-50,50,5);
   deconv_fitter_->SetNpx(1000);
   deconv_fitter_->FixParameter(0,0);
   deconv_fitter_->SetParLimits(1,-10,10);
   deconv_fitter_->SetParLimits(2,0,200);
   deconv_fitter_->SetParLimits(3,5,100);
   deconv_fitter_->FixParameter(3,50);
   deconv_fitter_->SetParLimits(4,0,50);
   deconv_fitter_->SetParameters(0.,-2.82,0.96,50,20);
}

// ----------------------------------------------------------------------------
// 
void SamplingAlgorithm::extract( const std::vector<TH1*>& histos) {

  if ( !anal() ) {
    edm::LogWarning(mlCommissioning_)
      << "[SamplingAlgorithm::" << __func__ << "]"
      << " NULL pointer to Analysis object!";
    return; 
  }
  
  CommissioningAnalysis* tmp = const_cast<CommissioningAnalysis*>( anal() );
  samp_ = dynamic_cast<SamplingAnalysis*>( tmp );
  if ( !samp_ ) {
    edm::LogWarning(mlCommissioning_)
      << "[SamplingAlgorithm::" << __func__ << "]"
      << " NULL pointer to derived Analysis object!";
    return; 
  }

  // Check
  if ( histos.size() != 1 && histos.size() != 2 ) {
    samp_->addErrorCode(sistrip::numberOfHistos_);
  }
  
  // Extract FED key from histo title
  if ( !histos.empty() ) { samp_->fedKey( extractFedKey( histos.front() ) ); }

  // Extract
  std::vector<TH1*>::const_iterator ihis = histos.begin();
  for ( ; ihis != histos.end(); ihis++ ) {
    
    // Check pointer
    if ( !(*ihis) ) {
      edm::LogWarning(mlCommissioning_) << " NULL pointer to histogram!";
      continue;
    }
    
    // Check name
    SiStripHistoTitle title( (*ihis)->GetName() );
    if ( title.runType() != sistrip::APV_LATENCY && title.runType() != sistrip::FINE_DELAY) {
      samp_->addErrorCode(sistrip::unexpectedTask_);
      continue;
    }
    // Set the mode for later fits
    samp_->runType_ = title.runType();
    
    // Set the granularity
    samp_->granularity_ = title.granularity();

    // Extract timing histo
    if ( title.extraInfo().find(sistrip::extrainfo::clusterCharge_) != std::string::npos ) {
      histo_.first = *ihis;
      histo_.second = (*ihis)->GetName();
    }

  }
  
}

// ----------------------------------------------------------------------------
// 
void SamplingAlgorithm::analyse() { 

  if ( !samp_ ) {
    edm::LogWarning(mlCommissioning_)
      << "[SamplingAlgorithm::" << __func__ << "]"
      << " NULL pointer to derived Analysis object!";
    return; 
  }

  TProfile* prof = (TProfile*)(histo_.first);
  if ( !prof ) {
    edm::LogWarning(mlCommissioning_)
      << " NULL pointer to histogram!" ;
    return;
  }

  // prune the profile
  pruneProfile(prof);

  // correct for the binning
  correctBinning(prof);

  // correct for clustering effects
  correctProfile(prof,samp_->sOnCut_);
  
  // fit depending on the mode
  if(samp_->runType_==sistrip::APV_LATENCY) {

    // initialize  the fit (overal latency)
    float max = prof->GetBinCenter(prof->GetMaximumBin());
    peak_fitter_->SetParameters(0.,50-max,10,50,10);

    // fit
    prof->Fit(peak_fitter_,"Q");
    prof->Fit(peak_fitter_,"QEWM");

    // Set monitorables
    samp_->max_   = peak_fitter_->GetMaximumX();
    samp_->error_ = peak_fitter_->GetParError(1);

  } else { // sistrip::FINE_DELAY

    // fit
    prof->Fit(deconv_fitter_,"Q");
    prof->Fit(deconv_fitter_,"QEWM");
    // Set monitorables
    samp_->max_   = deconv_fitter_->GetMaximumX();
    samp_->error_ = deconv_fitter_->GetParError(1);

  }

}

// ----------------------------------------------------------------------------
//
void SamplingAlgorithm::pruneProfile(TProfile* profile) const
{
  // loop over bins to find the max stat
  uint32_t nbins=profile->GetNbinsX();
  uint32_t max = 0;
  for(uint32_t bin=1;bin<=nbins;++bin) max = max < profile->GetBinEntries(bin) ? uint32_t(profile->GetBinEntries(bin)) : max;
  // loop over bins to clean
  uint32_t min = max/10;
  for(uint32_t bin=1;bin<=nbins;++bin)
    if(profile->GetBinEntries(bin)<min) {
      profile->SetBinContent(bin,0.);
      profile->SetBinError(bin,0.);
    }
}

// ----------------------------------------------------------------------------
//
void SamplingAlgorithm::correctBinning(TProfile* prof) const
{
  prof->GetXaxis()->SetLimits(prof->GetXaxis()->GetXmin()-prof->GetBinWidth(1)/2.,
                              prof->GetXaxis()->GetXmax()-prof->GetBinWidth(1)/2.);
}

// ----------------------------------------------------------------------------
//
void SamplingAlgorithm::correctProfile(TProfile* profile, float SoNcut) const
{
  if ( !samp_ ) { return; }
  uint32_t nbins=profile->GetNbinsX();
  float min = samp_->limit(SoNcut);
  for(uint32_t bin=1;bin<=nbins;++bin)
    if(profile->GetBinContent(bin)<min) {
      profile->SetBinContent(bin,0.);
      profile->SetBinError(bin,0.);
      profile->SetBinEntries(bin,0);
    }
    else {
      profile->SetBinContent(bin,
                             profile->GetBinEntries(bin)*
                             samp_->correctMeasurement(profile->GetBinContent(bin),SoNcut));
    }
}
