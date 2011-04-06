#include "TauAnalysis/CandidateTools/plugins/NSVfitAlgorithmByIntegration.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/CandidateTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include <limits>

using namespace SVfit_namespace;

double g(double* x, size_t dim, void* param)
{
  double nll = NSVfitAlgorithmBase::gNSVfitAlgorithm->nll(x, (double*)param);
  double retVal = TMath::Exp(-nll);
  static long callCounter = 0;
  //if ( (callCounter % 10000) == 0 ) 
  //  std::cout << "<g> (call = " << callCounter << "):" 
  //	        << " nll = " << nll << " --> returning retVal = " << retVal << std::endl;
  ++callCounter;
  return retVal;
}

NSVfitAlgorithmByIntegration::NSVfitAlgorithmByIntegration(const edm::ParameterSet& cfg)
  : NSVfitAlgorithmBase(cfg),
    fitParameterValues_(0),
    xl_(0),
    xu_(0),
    integrand_(0),
    workspace_(0),
    rnd_(0),
    numMassParameters_(0),
    massParForReplacements_(0)
{
  edm::ParameterSet cfg_replacements = cfg.getParameter<edm::ParameterSet>("parameters");
  std::vector<std::string> replacementNames = cfg_replacements.getParameterNamesForType<edm::ParameterSet>();
  for ( std::vector<std::string>::const_iterator replacementName = replacementNames.begin();
	replacementName != replacementNames.end(); ++replacementName ) {
    edm::ParameterSet cfg_replacement = cfg_replacements.getParameter<edm::ParameterSet>(*replacementName);
    
    fitParameterReplacementType* replacement = new fitParameterReplacementType();
    replacement->name_ = (*replacementName);
    replacement->iterLowerLimit_ = cfg_replacement.getParameter<double>("min");
    replacement->iterUpperLimit_ = cfg_replacement.getParameter<double>("max");
    replacement->iterStepSize_ = cfg_replacement.getParameter<double>("stepSize");

    std::string toReplace_string = cfg_replacement.getParameter<std::string>("replace");
    replacement->toReplace_ = toReplace_string;

    replacement->idxMassParameter_ = numMassParameters_;
    ++numMassParameters_;

    std::string replaceBy_string = cfg_replacement.getParameter<std::string>("by");
    size_t pos_token0 = -1;
    size_t pos = 0;
    std::set<std::string> tokens;
    while ( pos < replaceBy_string.length() ) {
      bool isSymbol = (replaceBy_string[pos] == '(' || replaceBy_string[pos] == ')' ||
		       replaceBy_string[pos] == '*' || replaceBy_string[pos] == '/' ||
		       replaceBy_string[pos] == '+' || replaceBy_string[pos] == '-');      
      if ( (isSymbol || pos == (replaceBy_string.length() - 1)) && (pos - pos_token0) > 1 ) {
	size_t num = ( pos != (replaceBy_string.length() - 1) ) ? pos - (pos_token0 + 1) : pos - pos_token0;
	std::string token = std::string(replaceBy_string, pos_token0 + 1, num);	
	tokens.insert(token);
      }
      if ( isSymbol ) pos_token0 = pos;
      ++pos;
    }

    std::string formula_string = replaceBy_string;
    int errorFlag = 0;    
    replacement->numParForReplacements_ = 0;
    for ( std::set<std::string>::const_iterator token = tokens.begin();
	  token != tokens.end(); ++token ) {
      if        ( (*token) == (*replacementName) ) {
	formula_string = replace_string(formula_string, *token, "x", 0, 1000, errorFlag);
      } else {
	size_t posSeparator = token->find(".");
	if ( posSeparator == std::string::npos ) {
	  throw cms::Exception("NSVfitAlgorithmByIntegration::NSVfitAlgorithmByIntegration")
	    << " Parameter token = " << (*token) << " has invalid format;" 
	    << " expected format is 'daughter.type' !!\n";
	}
	
	std::string particleName = std::string(*token, 0, posSeparator);
	std::string value_string = std::string(*token, posSeparator + 1);
	
	if        ( isDaughter(particleName)  ) {
	  replaceParByFitParameter* parForReplacement = new replaceParByFitParameter();	  
	  parForReplacement->fitParameterName_ = (*token);
	  parForReplacement->iPar_ = replacement->numParForReplacements_;
	  std::ostringstream par_string;
	  par_string << "[" << parForReplacement->iPar_ << "]";
	  formula_string = replace_string(formula_string, *token, par_string.str(), 0, 1000, errorFlag);
	  replacement->parForReplacements_.push_back(parForReplacement);
	  ++replacement->numParForReplacements_;
	} else if ( isResonance(particleName) ) {
	  replaceParByResonanceHypothesis* parForReplacement = new replaceParByResonanceHypothesis();	  
	  parForReplacement->resonanceName_ = particleName;
	  parForReplacement->valueExtractor_ = new StringObjectFunction<NSVfitResonanceHypothesis>(value_string);
	  parForReplacement->iPar_ = replacement->numParForReplacements_;
	  std::ostringstream par_string;
	  par_string << "[" << parForReplacement->iPar_ << "]";
	  formula_string = replace_string(formula_string, *token, par_string.str(), 0, 1000, errorFlag);
	  replacement->parForReplacements_.push_back(parForReplacement);
	  ++replacement->numParForReplacements_;
	} else {
	  throw cms::Exception("NSVfitAlgorithmByIntegration::NSVfitAlgorithmByIntegration")
	    << " No resonance/daughter of name = " << particleName << " defined !!\n";
	}
      } 
    }

    std::string formulaName = std::string(*replacementName).append("_formula");    
    replacement->replaceBy_ = new TFormula(formulaName.data(), formula_string.data());

    fitParameterReplacements_.push_back(replacement);
  }

  edm::ParameterSet cfg_vegas = cfg.getParameter<edm::ParameterSet>("vegasOptions");
  numCalls_ = cfg_vegas.getParameter<unsigned>("numCalls");

  massParForReplacements_ = new IndepCombinatoricsGeneratorT<double>(numMassParameters_);
  for ( unsigned iMassParameter = 0; iMassParameter < numMassParameters_; ++iMassParameter ) {
    const fitParameterReplacementType* fitParameterReplacement = fitParameterReplacements_[iMassParameter];
    massParForReplacements_->setLowerLimit(iMassParameter, fitParameterReplacement->iterLowerLimit_);
    massParForReplacements_->setUpperLimit(iMassParameter, fitParameterReplacement->iterUpperLimit_);
    massParForReplacements_->setStepSize(iMassParameter,   fitParameterReplacement->iterStepSize_);
  }
}

NSVfitAlgorithmByIntegration::~NSVfitAlgorithmByIntegration() 
{
  for( std::vector<fitParameterReplacementType*>::iterator it = fitParameterReplacements_.begin();
       it != fitParameterReplacements_.end(); ++it ) {
    delete (*it);
  }

  delete fitParameterValues_;

  delete [] xl_;
  delete [] xu_;

  if ( integrand_ ) delete [] (double*)integrand_->params;
  delete integrand_;
  if ( workspace_ ) gsl_monte_vegas_free(workspace_);
  if ( rnd_       ) gsl_rng_free(rnd_);

  delete massParForReplacements_;
}

void NSVfitAlgorithmByIntegration::beginJob()
{
  NSVfitAlgorithmBase::beginJob();
  
  for ( std::vector<fitParameterReplacementType*>::iterator fitParameterReplacement = fitParameterReplacements_.begin();
	fitParameterReplacement != fitParameterReplacements_.end(); ++fitParameterReplacement ) {
    (*fitParameterReplacement)->beginJob(this);
  } 

  numDimensions_ = 0;

  for ( std::vector<fitParameterType>::const_iterator fitParameter = fitParameters_.begin();
	fitParameter != fitParameters_.end(); ++fitParameter ) {
    bool isReplaced = false;
    for ( std::vector<fitParameterReplacementType*>::const_iterator fitParameterReplacement = fitParameterReplacements_.begin();
	  fitParameterReplacement != fitParameterReplacements_.end(); ++fitParameterReplacement ) {
      if ( fitParameter->idx_ == (*fitParameterReplacement)->idxToReplace_ ) isReplaced = true;
    }
    
    if ( !isReplaced ) {
      fitParameterByIntegrationType fitParameterByIntegration(&(*fitParameter));
      fitParameterByIntegration.idxByIntegration_ = numDimensions_;
      fitParametersByIntegration_.push_back(fitParameterByIntegration);
      ++numDimensions_;
    }
  }

  fitParameterValues_ = new double[fitParameters_.size()];

  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];

  integrand_ = new gsl_monte_function;
  integrand_->f = &g;
  integrand_->dim = numDimensions_;
  integrand_->params = new double[numMassParameters_];
  workspace_ = gsl_monte_vegas_alloc(numDimensions_);
  gsl_rng_env_setup();
  rnd_ = gsl_rng_alloc(gsl_rng_default);
}

void NSVfitAlgorithmByIntegration::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  NSVfitAlgorithmBase::beginEvent(evt, es);

  currentRunNumber_ = evt.id().run();
  currentLumiSectionNumber_ = evt.luminosityBlock();
  currentEventNumber_ = evt.id().event();
}

void getBinning(const IndepCombinatoricsGeneratorT<double>& grid, unsigned iDimension, int& numBins, double& min, double& max)
{
  double lowerLimit = grid.lowerLimit(iDimension);
  double upperLimit = grid.upperLimit(iDimension);
  double stepSize   = grid.stepSize(iDimension);
  numBins = 1 + (upperLimit - lowerLimit)/grid.stepSize(iDimension);
  min = lowerLimit - 0.5*(stepSize*numBins - (upperLimit - lowerLimit));
  max = min + numBins*stepSize;
}

void NSVfitAlgorithmByIntegration::fitImp() const
{
  //std::cout << "<NSVfitAlgorithmByIntegration::fitImp>:" << std::endl;

  for ( std::vector<fitParameterReplacementType*>::const_iterator fitParameterReplacement = fitParameterReplacements_.begin();
	fitParameterReplacement != fitParameterReplacements_.end(); ++fitParameterReplacement ) {
    for ( std::vector<replaceParBase*>::const_iterator par = (*fitParameterReplacement)->parForReplacements_.begin();
	  par != (*fitParameterReplacement)->parForReplacements_.end(); ++par ) {
      replaceParByResonanceHypothesis* par_resonance = dynamic_cast<replaceParByResonanceHypothesis*>(*par);
      if ( par_resonance ) {
	const NSVfitResonanceHypothesis* resonance = currentEventHypothesis_->resonance(par_resonance->resonanceName_);
	par_resonance->value_ = (*par_resonance->valueExtractor_)(*resonance);
      }
    }
  }

  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    xl_[iDimension] = fitParametersByIntegration_[iDimension].base_->lowerLimit_; 
    xu_[iDimension] = fitParametersByIntegration_[iDimension].base_->upperLimit_;
    //std::cout << " fitParameter #" << iDimension << " (" << fitParametersByIntegration_[iDimension].base_->name_ << ":" 
    //	        << fitParametersByIntegration_[iDimension].base_->type_ << "):"
    //	        << " xl = " << xl_[iDimension] << ", xu = " << xu_[iDimension] << std::endl;
  }

  gsl_monte_vegas_init(workspace_);
  massParForReplacements_->reset();

  TH1* histResults = 0;
  std::ostringstream histResultsName;
  histResultsName << pluginName_;
  histResultsName << "@" << currentRunNumber_ << ":" << currentLumiSectionNumber_ << ":" << currentEventNumber_;
  if ( numMassParameters_ == 1 ) {
    int numBinsX;
    double minX, maxX;
    getBinning(*massParForReplacements_, 0, numBinsX, minX, maxX);
    histResults = new TH1F(histResultsName.str().data(), histResultsName.str().data(), numBinsX, minX, maxX);
  } else if ( numMassParameters_ == 2 ) {
    int numBinsX, numBinsY;
    double minX, maxX, minY, maxY;
    getBinning(*massParForReplacements_, 0, numBinsX, minX, maxX);
    getBinning(*massParForReplacements_, 1, numBinsY, minY, maxY);
    histResults = new TH2F(histResultsName.str().data(), histResultsName.str().data(), numBinsX, minX, maxX, numBinsY, minY, maxY);    
  } else {
    throw cms::Exception("NSVfitAlgorithmByIntegration::fitImp")
      << " Only fits in one or two dimensions supported yet "
      << " --> please contact the developers Evan Friis (friis@physics.ucdavis.edu) and Christian Veelken (christian.veelken@cern.ch)"
      << " and request support for more dimensions !!\n";
  }

  while ( massParForReplacements_->isValid() ) {
//--- set mass parameters
    for ( unsigned iMassParameter = 0; iMassParameter < numMassParameters_; ++iMassParameter ) {
      ((double*)integrand_->params)[iMassParameter] = (*massParForReplacements_)[iMassParameter];
    }

//--- call VEGAS routine (part of GNU scientific library)
//    to perform actual integration
    double p, pErr;
    gsl_monte_vegas_integrate(integrand_, xl_, xu_, numDimensions_, numCalls_, rnd_, workspace_, &p, &pErr);
    //std::cout << "--> M = " << (*massParForReplacements_) << ": p = " << p << " +/- " << pErr << std::endl;

    if      ( numMassParameters_ == 1 ) histResults->Fill((*massParForReplacements_)[0], p);
    else if ( numMassParameters_ == 2 ) {
      TH2* histResults2d = dynamic_cast<TH2*>(histResults);
      assert(histResults2d);
      histResults2d->Fill((*massParForReplacements_)[0], (*massParForReplacements_)[1], p);
    } else assert(0);

    massParForReplacements_->next();
  }

//--- set central values and uncertainties on reconstructed masses
  for ( unsigned iMassParameter = 0; iMassParameter < numMassParameters_; ++iMassParameter ) {  
    const std::string& resonanceName = eventModel_->resonances_[iMassParameter]->resonanceName_;
    NSVfitResonanceHypothesis* resonance = const_cast<NSVfitResonanceHypothesis*>(currentEventHypothesis_->resonance(resonanceName));
    setMassResults(resonance, histResults, iMassParameter);
  }

  currentEventHypothesis_->histMassResults_ = std::auto_ptr<TH1>(histResults);
}

void NSVfitAlgorithmByIntegration::setMassResults(
       NSVfitResonanceHypothesis* resonance, const TH1* histMassResults, unsigned iDimension) const
{
  const TH1* histMassResult1d = 0;
  if      ( histMassResults->GetDimension() == 1 ) histMassResult1d = histMassResults;
  else if ( histMassResults->GetDimension() == 2 ) {
    const TH2* histMassResults2d = dynamic_cast<const TH2*>(histMassResults);
    assert(histMassResults2d);
    if      ( iDimension == 0 ) histMassResult1d = histMassResults2d->ProjectionX();
    else if ( iDimension == 1 ) histMassResult1d = histMassResults2d->ProjectionY();
  }
  assert(histMassResult1d);

//--- compute median, -1 sigma and +1 sigma limits on reconstructed mass
  if ( histMassResult1d->Integral() > 0. ) {
    Double_t q[3];
    Double_t probSum[3];
    probSum[0] = 0.16;
    probSum[1] = 0.50;
    probSum[2] = 0.84;
    (const_cast<TH1*>(histMassResult1d))->GetQuantiles(3, q, probSum);

    int binMaximum = histMassResult1d->GetMaximumBin();
    double massMaxInterpol = 0.;
    if ( binMaximum > 1 && binMaximum < histMassResult1d->GetNbinsX() ) {
      double xMaximum = histMassResult1d->GetBinCenter(binMaximum);
      double yMaximum = histMassResult1d->GetBinContent(binMaximum);
      
      int binLeft     = binMaximum - 1;
      double xLeft    = histMassResult1d->GetBinCenter(binLeft);
      double yLeft    = histMassResult1d->GetBinContent(binLeft);    
      
      int binRight    = binMaximum + 1;
      double xRight   = histMassResult1d->GetBinCenter(binRight);
      double yRight   = histMassResult1d->GetBinContent(binRight); 
      
      double xMinus   = xLeft - xMaximum;
      double yMinus   = yLeft - yMaximum;
      double xPlus    = xRight - xMaximum;
      double yPlus    = yRight - yMaximum;
      
      massMaxInterpol = xMaximum + 0.5*(yPlus*square(xMinus) - yMinus*square(xPlus))/(yPlus*xMinus - yMinus*xPlus);
    } else {
      massMaxInterpol = histMassResult1d->GetBinCenter(binMaximum);
    }

    NSVfitAlgorithmBase::setMassResults(resonance, q[1], TMath::Abs(q[2] - q[1]), TMath::Abs(q[1] - q[0]));

    resonance->massMean_ = histMassResult1d->GetMean();
    resonance->massMedian_ = q[1];
    resonance->massMaximum_ = histMassResult1d->GetBinCenter(binMaximum);
    resonance->massMaxInterpol_ = massMaxInterpol;
    resonance->isValidSolution_ = false;
    
    //std::cout << "<NSVfitAlgorithmByIntegration::setMassResults>:" << std::endl;
    //std::cout << "--> mass = " << resonance->mass_ << std::endl;
  } else {
    edm::LogWarning("NSVfitAlgorithmByIntegration::setMassResults")
      << "Likelihood functions returned Probability zero for all tested mass hypotheses --> no valid solution found !!";
    resonance->isValidSolution_ = false;
  }
  
  if ( histMassResult1d != histMassResults ) delete histMassResult1d;
}

bool NSVfitAlgorithmByIntegration::isDaughter(const std::string& daughterName)
{
  bool isDaughter = false;

  for ( std::vector<resonanceModelType*>::const_iterator resonance = eventModel_->resonances_.begin();
	resonance != eventModel_->resonances_.end(); ++resonance ) {
    for ( std::vector<daughterModelType*>::const_iterator daughter = (*resonance)->daughters_.begin();
	  daughter != (*resonance)->daughters_.end(); ++daughter ) {
      if ( (*daughter)->daughterName_ == daughterName ) isDaughter = true;
    }
  }

  return isDaughter;
}
 
bool NSVfitAlgorithmByIntegration::isResonance(const std::string& resonanceName)
{
  bool isResonance = false;

  for ( std::vector<resonanceModelType*>::const_iterator resonance = eventModel_->resonances_.begin();
	resonance != eventModel_->resonances_.end(); ++resonance ) {
    if ( (*resonance)->resonanceName_ == resonanceName ) isResonance = true;
  }

  return isResonance;
}

NSVfitAlgorithmBase::fitParameterType* NSVfitAlgorithmByIntegration::getFitParameter(const std::string& token)
{
  size_t posSeparator = token.find(".");
  if ( posSeparator == std::string::npos || posSeparator == (token.length() - 1) ) {
    throw cms::Exception("NSVfitAlgorithmByIntegration::getFitParameter")
      << " Parameter token = " << token << " passed as function argument has invalid format;" 
      << " expected format is 'daughter.type' !!\n";
  }

  std::string name = std::string(token, 0, posSeparator);
  std::string type_string = std::string(token, posSeparator + 1);
  int type = -1;
  if ( type_string == "x" ) type = nSVfit_namespace::kTau_visEnFracX;
  else throw cms::Exception("NSVfitAlgorithmByIntegration::getFitParameter")
    << " Type = " << type << " not defined !!\n";

  return NSVfitAlgorithmBase::getFitParameter(name, type);
}

double NSVfitAlgorithmByIntegration::nll(double* x, double* param) const
{
//--- copy fitParameter
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    fitParameterValues_[fitParametersByIntegration_[iDimension].base_->idx_] = x[iDimension];
  }

//--- set additional fitParameters according to mass parameter values
  for ( std::vector<fitParameterReplacementType*>::const_iterator fitParameterReplacement = fitParameterReplacements_.begin();
	fitParameterReplacement != fitParameterReplacements_.end(); ++fitParameterReplacement ) {
    TFormula* formula = (*fitParameterReplacement)->replaceBy_;

    for ( int iPar = 0; iPar < (*fitParameterReplacement)->numParForReplacements_; ++iPar ) {
      formula->SetParameter(iPar, (*(*fitParameterReplacement)->parForReplacements_[iPar])(fitParameterValues_));
    }

//--- check if fitParameter is within limits;
//    return probability zero if not
    double fitParameterValue = formula->Eval(param[(*fitParameterReplacement)->idxMassParameter_]);
    if ( fitParameterValue >= fitParameterLimits_[(*fitParameterReplacement)->idxToReplace_].first  &&
	 fitParameterValue <= fitParameterLimits_[(*fitParameterReplacement)->idxToReplace_].second ) {
      fitParameterValues_[(*fitParameterReplacement)->idxToReplace_] = fitParameterValue;
    } else {
      return std::numeric_limits<float>::max();
    }
  }

//--- build event, resonance and particle hypotheses
  eventModel_->builder_->applyFitParameter(currentEventHypothesis_, fitParameterValues_);

//--- compute likelihood
  return eventModel_->nll(currentEventHypothesis_);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitAlgorithmPluginFactory, NSVfitAlgorithmByIntegration, "NSVfitAlgorithmByIntegration");


