#ifndef TauAnalysis_CandidateTools_NSVfitEventLikelihoodMEt2_h
#define TauAnalysis_CandidateTools_NSVfitEventLikelihoodMEt2_h

/** \class NSVfitEventLikelihoodMEt2
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * New version using covariance matrix of (PF)MET significance calculation
 * (CMS AN-10/400) to compute the likehood
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: NSVfitEventLikelihoodMEt2.h,v 1.3 2011/05/26 13:06:36 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoMET/METAlgorithms/interface/SignAlgoResolutions.h"
#include "RecoMET/METAlgorithms/interface/SigInputObj.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitEventLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/PFMEtSignInterface.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"

#include <TMatrixD.h>

#include <list>

class NSVfitEventLikelihoodMEt2 : public NSVfitEventLikelihood
{
 public:
  NSVfitEventLikelihoodMEt2(const edm::ParameterSet&);
  ~NSVfitEventLikelihoodMEt2();

  void beginJob(NSVfitAlgorithmBase*);
  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const NSVfitEventHypothesis*) const;

  double operator()(const NSVfitEventHypothesis*) const;

 private:

  double power_;

  PFMEtSignInterface* pfMEtSign_;

  //mutable TMatrixD pfMEtCovInverse_;

  mutable ROOT::Math::SMatrix<double, 2, 2> pfMEtCovInverse_;

  // Vector which holds the current residual.  We make this a member to prevent
  // excessive memory allocations.
  mutable ROOT::Math::SVector<double, 2> residual_fitted_;
};

#endif
