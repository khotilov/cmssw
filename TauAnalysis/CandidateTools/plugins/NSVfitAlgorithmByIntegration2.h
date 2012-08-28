#ifndef TauAnalysis_CandidateTools_NSVfitAlgorithmByIntegration2_h
#define TauAnalysis_CandidateTools_NSVfitAlgorithmByIntegration2_h

/** \class SVfitAlgorithmByIntegration2
 *
 * Concrete implementation of (n)SVfit algorithm
 * by integration of likelihood functions
 * (based on Markov Chain integration algorithm)
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitAlgorithmByIntegration2.h,v 1.1 2012/04/07 15:44:43 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/MarkovChainIntegrator.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <Math/Functor.h>
#include <TH1.h>

#include <vector>
#include <string>

class NSVfitAlgorithmByIntegration2 : public NSVfitAlgorithmBase
{
 public:
  NSVfitAlgorithmByIntegration2(const edm::ParameterSet&);
  ~NSVfitAlgorithmByIntegration2();

  void beginJob();
  void beginEvent(const edm::Event&, const edm::EventSetup&);

  unsigned getNumDimensions() const { return numDimensions_; }

  void print(std::ostream&) const {}

  bool update(const double* x, const double* param) const;
  double nll(const double* x, const double* param) const;

  void fillProbHistograms(const double*);

 protected:
  TH1* bookMassHistogram(const std::string& name);

  void fitImp() const;

  void setMassResults(NSVfitResonanceHypothesisBase*, const TH1*) const;
    
  struct NSVfitParameterMappingType
  {
    NSVfitParameterMappingType(const NSVfitParameter* base)
      : base_(base)
    {}
    const NSVfitParameter* base_;
    int idxByIntegration_;
  };

  std::vector<NSVfitParameterMappingType> fitParameterMappings_;

  std::vector<NSVfitParameterMappingType> constParameterMappings_;
  unsigned numConstParameters_;

  ROOT::Math::Functor* integrand_;
  MarkovChainIntegrator* integrator_;
  ROOT::Math::Functor* auxPhysicalSolutionFinder_;
  std::vector<double> intBoundaryLower_;
  std::vector<double> intBoundaryUpper_;
  unsigned numDimensions_;

  bool monitorMarkovChain_;  
  ROOT::Math::Functor* auxResonanceMassValue_;

  double* fitParameterValues_;

  std::vector<TH1*> probHistFitParameter_;
  std::map<std::string, TH1*> probHistResonanceMass_;
  TH1* probHistEventMass_;
  mutable std::vector<double> probListEventMass_;
  ROOT::Math::Functor* auxFillProbHistograms_;

  edm::RunNumber_t currentRunNumber_;
  edm::LuminosityBlockNumber_t currentLumiSectionNumber_;
  edm::EventNumber_t currentEventNumber_;
};

#endif

