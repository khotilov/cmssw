#include "TauAnalysis/CandidateTools/plugins/NSVfitAlgorithmByIntegration2.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TArrayF.h>
#include <TH1D.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>

using namespace SVfit_namespace;

enum { kMax, kMedian };

namespace 
{
  class Integrand : public ROOT::Math::Functor
  {
   public:

    Integrand(NSVfitAlgorithmByIntegration2* algorithm)
      : algorithm_(algorithm)
    {}

    unsigned int NDim() const 
    {
      return algorithm_->getNumDimensions();
    }

   private:

    virtual double DoEval(const double* x) const
    {
      double nll = algorithm_->nll(x, 0);
      double retVal = TMath::Exp(-nll);
      //static long callCounter = 0;
      //if ( (callCounter % 10000) == 0 ) {
      //  std::cout << "<Integrand> (call = " << callCounter << "):" << std::endl;
      //  unsigned numDimensions = NDim();
      //  for ( unsigned iDimension = 0; iDimension < numDimensions; ++iDimension ) {
      //    std::cout << " x[" << iDimension << "] = " << x[iDimension] << std::endl;
      //  }
      //  std::cout << "nll = " << nll << " --> returning retVal = " << retVal << std::endl;
      //}
      //++callCounter;
      return retVal;
    } 

    NSVfitAlgorithmByIntegration2* algorithm_;
  };

  class AuxPhysicalSolutionFinder : public ROOT::Math::Functor
  {
   public:

    AuxPhysicalSolutionFinder(NSVfitAlgorithmByIntegration2* algorithm)
      : algorithm_(algorithm)
    {}

   private:

    virtual double DoEval(const double* x) const
    {
      bool isPhysicalSolution = algorithm_->update(x, 0);
      return ( isPhysicalSolution ) ? 1. : 0.;
    } 

    NSVfitAlgorithmByIntegration2* algorithm_;
  };

  class AuxFillProbHistograms : public ROOT::Math::Functor
  {
   public:

    AuxFillProbHistograms(NSVfitAlgorithmByIntegration2* algorithm)
      : algorithm_(algorithm)
    {}

   private:

    virtual double DoEval(const double* x) const
    {
      algorithm_->fillProbHistograms(x);
      return 0.;
    } 

    NSVfitAlgorithmByIntegration2* algorithm_;
  };

  class AuxResonanceMassValue : public ROOT::Math::Functor
  {
   public:

    AuxResonanceMassValue(NSVfitAlgorithmByIntegration2* algorithm)
      : algorithm_(algorithm)
    {}

   private:

    virtual double DoEval(const double* x) const
    {
      const NSVfitResonanceHypothesis* resonance = algorithm_->currentEventHypothesis()->resonance(0);
      assert(resonance);
      return resonance->p4_fitted().mass();
    } 

    NSVfitAlgorithmByIntegration2* algorithm_;
  };
}

NSVfitAlgorithmByIntegration2::NSVfitAlgorithmByIntegration2(const edm::ParameterSet& cfg)
  : NSVfitAlgorithmBase(cfg),
    integrand_(0),
    integrator_(0),
    auxPhysicalSolutionFinder_(0),
    auxResonanceMassValue_(0),
    fitParameterValues_(0),
    probHistEventMass_(0),        
    auxFillProbHistograms_(0)
{
  edm::ParameterSet cfgMarkovChainOptions = cfg.getParameter<edm::ParameterSet>("markovChainOptions");
  cfgMarkovChainOptions.addParameter<std::string>("name", pluginName_);
  cfgMarkovChainOptions.addParameter<int>("verbosity", verbosity_);
  integrator_ = new MarkovChainIntegrator(cfgMarkovChainOptions);

  monitorMarkovChain_ = ( cfg.exists("monitorMarkovChain") ) ?
    cfg.getParameter<bool>("monitorMarkovChain") : false;

  std::string max_or_median_string = cfg.getParameter<std::string>("max_or_median");
  if      ( max_or_median_string == "max"    ) max_or_median_ = kMax;
  else if ( max_or_median_string == "median" ) max_or_median_ = kMedian;
  else throw cms::Exception("NSVfitAlgorithmByIntegration2")
    << " Invalid Configuration Parameter 'max_or_median' = " << max_or_median_string << " !!\n";
}

NSVfitAlgorithmByIntegration2::~NSVfitAlgorithmByIntegration2() 
{
  delete integrand_;
  delete integrator_;
  delete auxPhysicalSolutionFinder_;

  delete [] fitParameterValues_;

  for ( std::vector<TH1*>::iterator it = probHistFitParameter_.begin();
	it != probHistFitParameter_.end(); ++it ) {
    delete (*it);
  }

  for ( std::map<std::string, TH1*>::iterator it = probHistResonanceMass_.begin();
	it != probHistResonanceMass_.end(); ++it ) {
    delete it->second;
  }
  delete probHistEventMass_;
  delete auxFillProbHistograms_;

  delete auxResonanceMassValue_;
}

void NSVfitAlgorithmByIntegration2::beginJob()
{
  NSVfitAlgorithmBase::beginJob();

  numDimensions_ = 0;
  numConstParameters_ = 0;

  for ( std::vector<NSVfitParameter>::const_iterator fitParameter = fitParameters_.begin();
	fitParameter != fitParameters_.end(); ++fitParameter ) {
    bool isFixed = fitParameter->IsFixed();
    
    if ( !isFixed ) {
      NSVfitParameterMappingType fitParameterMapping(&(*fitParameter));
      fitParameterMapping.idxByIntegration_ = numDimensions_;
      fitParameterMappings_.push_back(fitParameterMapping);
      ++numDimensions_;
    }

    if ( isFixed ) {
      NSVfitParameterMappingType fitParameterMapping(&(*fitParameter));
      fitParameterMapping.idxByIntegration_ = numConstParameters_;
      constParameterMappings_.push_back(fitParameterMapping);
      ++numConstParameters_;
    }
  }

  fitParameterValues_ = new double[fitParameters_.size()];

  integrand_ = new Integrand(this);
  integrator_->setIntegrand(*integrand_);

  auxPhysicalSolutionFinder_ = new AuxPhysicalSolutionFinder(this);
  integrator_->setStartPosition_and_MomentumFinder(*auxPhysicalSolutionFinder_);

  startPosition_.resize(numDimensions_);
  intBoundaryLower_.resize(numDimensions_);
  intBoundaryUpper_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const NSVfitParameter* fitParameter_ref = fitParameterMappings_[iDimension].base_;
    intBoundaryLower_[iDimension] = fitParameter_ref->LowerLimit();
    intBoundaryUpper_[iDimension] = fitParameter_ref->UpperLimit();
  }

  probHistFitParameter_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const NSVfitParameter* fitParameter_ref = fitParameterMappings_[iDimension].base_;
    const std::string& fitParameterName = fitParameter_ref->UniqueName();
    std::string histogramName = std::string(pluginName_).append("_").append(fitParameterName);
    TH1* histogram = new TH1D(histogramName.data(), histogramName.data(), 1000, fitParameter_ref->LowerLimit(), fitParameter_ref->UpperLimit());
    probHistFitParameter_[iDimension] = histogram;
  }
  for ( std::vector<resonanceModelType*>::const_iterator resonance = eventModel_->resonances_.begin();
	resonance != eventModel_->resonances_.end(); ++resonance ) {
    const std::string& resonanceName = (*resonance)->resonanceName_;
    std::string histogramName = std::string("probHistResonanceMass").append("_").append(resonanceName);
    TH1* histogram = bookMassHistogram(histogramName.data());
    probHistResonanceMass_.insert(std::pair<std::string, TH1*>(resonanceName, histogram));
  }
  probHistEventMass_ = bookMassHistogram("probHistEventMass");
  auxFillProbHistograms_ = new AuxFillProbHistograms(this);
  integrator_->registerCallBackFunction(*auxFillProbHistograms_);

  if ( monitorMarkovChain_ ) {
    auxResonanceMassValue_ = new AuxResonanceMassValue(this);
    integrator_->setF(*auxResonanceMassValue_);
  }
}

void NSVfitAlgorithmByIntegration2::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  NSVfitAlgorithmBase::beginEvent(evt, es);

  for ( std::vector<TH1*>::iterator histogram = probHistFitParameter_.begin();
	histogram != probHistFitParameter_.end(); ++histogram ) {
    (*histogram)->Reset();
  }

  for ( std::map<std::string, TH1*>::iterator histogram = probHistResonanceMass_.begin();
	histogram != probHistResonanceMass_.end(); ++histogram ) {
    histogram->second->Reset();
  }

  probHistEventMass_->Reset();
  probListEventMass_.clear();

  currentRunNumber_ = evt.id().run();
  currentLumiSectionNumber_ = evt.luminosityBlock();
  currentEventNumber_ = evt.id().event();
}

TH1* NSVfitAlgorithmByIntegration2::bookMassHistogram(const std::string& histogramName)
{
  double xMin = 1.e+1;
  double xMax = 1.e+4;
  double logBinWidth = 1.01;
  int numBins = TMath::Log(xMax/xMin)/TMath::Log(logBinWidth);
  TArrayF binning(numBins + 1);
  double x = xMin;  
  for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    binning[iBin] = x;
    //std::cout << "binning[" << iBin << "] = " << binning[iBin] << std::endl;
    x *= logBinWidth;
  }  
  std::string histogramName_full = std::string(pluginName_).append("_").append(histogramName);
  TH1* histogram = new TH1D(histogramName_full.data(), histogramName_full.data(), numBins, binning.GetArray());
  return histogram;
}

void NSVfitAlgorithmByIntegration2::fitImp() const
{
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const NSVfitParameter* fitParameter_ref = fitParameterMappings_[iDimension].base_; 
    startPosition_[iDimension] = fitParameter_ref->InitialValue();
    intBoundaryUpper_[iDimension] = fitParameter_ref->UpperLimit(); 
    intBoundaryLower_[iDimension] = fitParameter_ref->LowerLimit();
  }
  if ( verbosity_ >= 2 ) std::cout << "startPosition = " << format_vdouble(startPosition_) << std::endl;

  // CV: transform startPosition into interval ]-1..+1[
  //     expected by MarkovChainIntegrator class
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    double x_i = startPosition_[iDimension];
    double xMax = intBoundaryUpper_[iDimension];
    double xMin = intBoundaryLower_[iDimension];
    assert(xMax > xMin);
    startPosition_[iDimension] = (x_i - xMin)/(xMax - xMin);
  }
  if ( verbosity_ >= 2 ) std::cout << "startPosition (mapped into interval ]-1..+1[) = " << format_vdouble(startPosition_) << std::endl;
  integrator_->initializeStartPosition_and_Momentum(startPosition_);

  double integral, integralErr;
  int errorFlag = 0;
  std::string monitorFileName;
  if ( monitorMarkovChain_ ) {
    monitorFileName = Form("%s_run%i_ls%i_ev%i_mc.root", pluginName_.data(), currentRunNumber_, currentLumiSectionNumber_, currentEventNumber_);
  } else {
    monitorFileName = "";
  }
  integrator_->integrate(intBoundaryLower_, intBoundaryUpper_, integral, integralErr, errorFlag, monitorFileName);

//--- set central values and uncertainties on reconstructed masses
  for ( std::vector<resonanceModelType*>::const_iterator resonance = eventModel_->resonances_.begin();
	resonance != eventModel_->resonances_.end(); ++resonance ) {
    const std::string& resonanceName = (*resonance)->resonanceName_;
    NSVfitResonanceHypothesis* resonance = const_cast<NSVfitResonanceHypothesis*>(currentEventHypothesis_->resonance(resonanceName));
    assert(resonance);
    std::map<std::string, TH1*>::const_iterator histogram = probHistResonanceMass_.find(resonanceName);
    assert(histogram != probHistResonanceMass_.end());
    if ( errorFlag == 0 ) setMassResults(resonance, histogram->second);
    else resonance->isValidSolution_ = false;
  }

  if ( verbosity_ >= 2 ) {
    currentEventHypothesis_->print(std::cout);
    std::vector<double> quantiles;
    quantiles.push_back(0.16);
    quantiles.push_back(0.50);
    quantiles.push_back(0.84);
    std::sort(probListEventMass_.begin(), probListEventMass_.end());
    for ( std::vector<double>::const_iterator quantile = quantiles.begin();
	  quantile != quantiles.end(); ++quantile ) {
      int idx = TMath::Nint((*quantile)*probListEventMass_.size());
      std::cout << "probListEventMass[" << idx << "] (" << (*quantile)*100. << "% quantile) = " << probListEventMass_[idx] << std::endl;
    }
  }

  fittedEventHypothesis_ = currentEventHypothesis_;  
  if ( errorFlag == 0 ) {
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      double valueMaximum, valueMaximum_interpol, valueMean, valueQuantile016, valueQuantile050, valueQuantile084;
      extractHistogramProperties(
        probHistFitParameter_[iDimension], probHistFitParameter_[iDimension],
	valueMaximum, valueMaximum_interpol, valueMean, valueQuantile016, valueQuantile050, valueQuantile084);
      if      ( max_or_median_ == kMax    ) fitParameterValues_[iDimension] = valueMaximum_interpol;
      else if ( max_or_median_ == kMedian ) fitParameterValues_[iDimension] = valueQuantile050;
      else assert(0);
    }
    fittedEventHypothesis_nll_ = this->nll(fitParameterValues_, 0);
  } else {
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      fitParameterValues_[iDimension] = 0.;
    }
    fittedEventHypothesis_nll_ = -1.;
  }
}

void NSVfitAlgorithmByIntegration2::setMassResults(NSVfitResonanceHypothesisBase* resonance, const TH1* histMassResult) const
{
  std::string histMassResultName_density = std::string(histMassResult->GetName()).append("_cloned");
  TH1* histMassResult_density = (TH1*)histMassResult->Clone(histMassResultName_density.data());
  for ( int iBin = 1; iBin <= histMassResult->GetNbinsX(); ++iBin ) {
    double binContent = histMassResult->GetBinContent(iBin);
    double binError = histMassResult->GetBinError(iBin);
    double binWidth = histMassResult->GetBinWidth(iBin);
    assert(binWidth > 0.);
    histMassResult_density->SetBinContent(iBin, binContent/binWidth);
    histMassResult_density->SetBinError(iBin, binError/binWidth);
  }

  if ( histMassResult_density->Integral() > 0. ) {
    double massMaximum, massMaximum_interpol, massMean, massQuantile016, massQuantile050, massQuantile084;
    extractHistogramProperties(
      histMassResult, histMassResult_density,
      massMaximum, massMaximum_interpol, massMean, massQuantile016, massQuantile050, massQuantile084);
    
    double mass;
    if      ( max_or_median_ == kMax    ) mass = massMaximum_interpol;
    else if ( max_or_median_ == kMedian ) mass = massQuantile050;
    else assert(0);
    double massErrUp   = TMath::Abs(massQuantile084 - mass);
    double massErrDown = TMath::Abs(mass - massQuantile016);
    NSVfitAlgorithmBase::setMassResults(resonance, mass, massErrUp, massErrDown);
  
    resonance->isValidSolution_ = true;
    
    if ( verbosity_ >= 1 ) {
      std::cout << "<NSVfitAlgorithmByIntegration2::setMassResults>:" << std::endl;
      std::cout << " pluginName = " << pluginName_ << std::endl;
      std::cout << "--> mass = " << resonance->mass_ << std::endl;
      std::cout << " (mean = " << massMean << ", median = " << massQuantile050 << ", max = " << massMaximum << ")" << std::endl;
    }
  } else {
    edm::LogWarning("NSVfitAlgorithmByIntegration2::setMassResults")
      << "Likelihood functions returned Probability zero for all tested mass hypotheses --> no valid solution found !!";
    resonance->isValidSolution_ = false;
  }
  
  delete histMassResult_density;
}

bool NSVfitAlgorithmByIntegration2::update(const double* x, const double* param) const
{
//--- copy fitParameter
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const NSVfitParameter* fitParameter_ref = fitParameterMappings_[iDimension].base_;
    fitParameterValues_[fitParameter_ref->index()] = x[iDimension];
  }

//--- copy constant parameters
  for ( unsigned iConstParameter = 0; iConstParameter < numConstParameters_; ++iConstParameter ) {
    const NSVfitParameter* fitParameter_ref = constParameterMappings_[iConstParameter].base_;
    fitParameterValues_[fitParameter_ref->index()] = fitParameter_ref->Value();    
  }

//--- build event, resonance and particle hypotheses
//    and check if hypothesis corresponds to a "valid" (physically allowed) solution
  currentEventHypothesis_isValidSolution_ = eventModel_->builder_->applyFitParameter(currentEventHypothesis_, fitParameterValues_);
  if ( verbosity_ >= 2 ) {
    currentEventHypothesis_->print(std::cout);
    std::cout << "isValidSolution = " << currentEventHypothesis_isValidSolution_ << std::endl;
  }
  return currentEventHypothesis_isValidSolution_;
}

double NSVfitAlgorithmByIntegration2::nll(const double* x, const double* param) const
{
  bool isPhysicalSolution = update(x, param);

  double nll = std::numeric_limits<float>::max();
  if ( isPhysicalSolution ) {
    nll = eventModel_->nll(currentEventHypothesis_);
  }
  if ( TMath::IsNaN(nll) ) nll = std::numeric_limits<float>::max();

  return nll;
}

void NSVfitAlgorithmByIntegration2::fillProbHistograms(const double* x)
{
//--- fill histograms of fitParameter distributions
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    probHistFitParameter_[iDimension]->Fill(x[iDimension]);
  }

//--- fill mass distribution histograms
  for ( std::vector<resonanceModelType*>::const_iterator resonance = eventModel_->resonances_.begin();
	resonance != eventModel_->resonances_.end(); ++resonance ) {
    const std::string& resonanceName = (*resonance)->resonanceName_;
    const NSVfitResonanceHypothesis* resonance = currentEventHypothesis_->resonance(resonanceName);
    assert(resonance);
    TH1* histogram = probHistResonanceMass_[resonanceName];
    assert(histogram);
    histogram->Fill(resonance->p4_fitted().mass());
  }
  probHistEventMass_->Fill(currentEventHypothesis_->p4_fitted().mass());  
  probListEventMass_.push_back(currentEventHypothesis_->p4_fitted().mass());  
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitAlgorithmPluginFactory, NSVfitAlgorithmByIntegration2, "NSVfitAlgorithmByIntegration2");


