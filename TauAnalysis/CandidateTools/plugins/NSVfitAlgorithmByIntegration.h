#ifndef TauAnalysis_CandidateTools_NSVfitAlgorithmByIntegration_h
#define TauAnalysis_CandidateTools_NSVfitAlgorithmByIntegration_h

/** \class SVfitAlgorithmByIntegration
 *
 * Concrete implementation of (n)SVfit algorithm
 * by integration of likelihood functions
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitAlgorithmByIntegration.h,v 1.1 2011/03/03 13:04:47 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/IndepCombinatoricsGeneratorT.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <TFormula.h>

#include <vector>
#include <string>

class NSVfitAlgorithmByIntegration : public NSVfitAlgorithmBase
{
 public:
  NSVfitAlgorithmByIntegration(const edm::ParameterSet&);
  ~NSVfitAlgorithmByIntegration();

  void beginJob();

  void print(std::ostream&) const {}

  double nll(double*, double*) const;

 protected:
  void fitImp() const;

  bool isDaughter(const std::string&);
  bool isResonance(const std::string&);
  
  NSVfitAlgorithmBase::fitParameterType* getFitParameter(const std::string&);

  struct fitParameterByIntegrationType
  {
    fitParameterByIntegrationType(const NSVfitAlgorithmBase::fitParameterType* base)
      : base_(base)
    {}
    const NSVfitAlgorithmBase::fitParameterType* base_;
    int idxByIntegration_;
  };

  std::vector<fitParameterByIntegrationType> fitParametersByIntegration_;

  struct replaceParBase
  {
    virtual void beginJob(NSVfitAlgorithmByIntegration*) {}
    virtual double operator()(double*) const = 0;
    int iPar_;
  };

  struct replaceParByFitParameter : replaceParBase
  {
    void beginJob(NSVfitAlgorithmByIntegration* algorithm)
    {
      idx_ = algorithm->getFitParameter(fitParameterName_)->idx_;
    }
    double operator()(double* param) const { return param[idx_]; }
    std::string fitParameterName_;
    int idx_;
  };

  struct replaceParByResonanceHypothesis : replaceParBase
  {
    replaceParByResonanceHypothesis()
      : valueExtractor_(0)
    {}
    ~replaceParByResonanceHypothesis()
    {
      delete valueExtractor_;
    }
    double operator()(double* param) const { return value_; }
    std::string resonanceName_;
    NSVfitResonanceHypothesis* resonanceHypothesis_;
    StringObjectFunction<NSVfitResonanceHypothesis>* valueExtractor_;
    mutable double value_;
  };

  struct fitParameterReplacementType
  {    
    fitParameterReplacementType()
      : replaceBy_(0)
    {}
    ~fitParameterReplacementType() 
    {
      delete replaceBy_;
      for ( std::vector<replaceParBase*>::iterator it = parForReplacements_.begin();
	    it != parForReplacements_.end(); ++it ) {
	delete (*it);
      }
    }
    void beginJob(NSVfitAlgorithmByIntegration* algorithm)
    {
      NSVfitAlgorithmBase::fitParameterType* fitParameterToReplace = algorithm->getFitParameter(toReplace_);
      if ( !fitParameterToReplace ) {
	throw cms::Exception("fitParameterReplacementType::beginJob")
	  << " No fitParameter of name = " << toReplace_ << " defined !!";
      }
      idxToReplace_ = fitParameterToReplace->idx_;

      for ( std::vector<replaceParBase*>::iterator par = parForReplacements_.begin();
	    par != parForReplacements_.end(); ++par ) {
	(*par)->beginJob(algorithm);
      }
    }
    std::string name_;
    double iterLowerLimit_;
    double iterUpperLimit_;
    double iterStepSize_;
    std::string toReplace_;
    int idxToReplace_;
    TFormula* replaceBy_;
    int idxMassParameter_;
    std::vector<replaceParBase*> parForReplacements_;
    int numParForReplacements_;
  };

  std::vector<fitParameterReplacementType*> fitParameterReplacements_;

  double* fitParameterValues_;

  double* xl_;
  double* xu_;

  gsl_monte_function* integrand_;
  gsl_monte_vegas_state* workspace_;
  gsl_rng* rnd_;
  unsigned numCalls_;
  size_t numDimensions_;

  size_t numMassParameters_;
  IndepCombinatoricsGeneratorT<double>* massParForReplacements_;
};

#endif

