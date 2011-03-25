#ifndef TauAnalysis_CandidateTools_NSVfitEventLikelihoodTrackInfo_h
#define TauAnalysis_CandidateTools_NSVfitEventLikelihoodTrackInfo_h

/** \class NSVfitEventLikelihoodTrackInfo
 *
 * Plugin for computing likelihood for position of primary event vertex 
 * refitted after excluding tracks of particles produced in tau lepton decays 
 * to be compatible with reconstructed tau lepton decay kinematics
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitEventLikelihoodTrackInfo.h,v 1.1 2011/03/23 17:46:39 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitEventLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesis.h"

class NSVfitEventLikelihoodTrackInfo : public NSVfitEventLikelihood
{
 public:
  NSVfitEventLikelihoodTrackInfo(const edm::ParameterSet&);
  ~NSVfitEventLikelihoodTrackInfo();

  void beginJob(NSVfitAlgorithmBase*) const; 

  double operator()(const NSVfitEventHypothesis*) const;
};

#endif
