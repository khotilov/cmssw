#include "TauAnalysis/CandidateTools/interface/NSVfitTauDecayBuilderBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToHadHypothesis.h"
#include "TauAnalysis/CandidateTools/interface/nSVfitParameter.h"

/** \class NSVfitTauToHadBuilder
 *
 * Auxiliary class reconstructing tau --> had decays and
 * building NSVfitTauToHadHypothesis objects;
 * used by NSVfit algorithm
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.16 $
 *
 * $Id: NSVfitTauToHadBuilder.cc,v 1.16 2011/03/31 09:29:25 veelken Exp $
 *
 */

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitTrackService.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

using namespace nSVfit_namespace;

class NSVfitTauToHadBuilder : public NSVfitTauDecayBuilderBase
{
 public:
  NSVfitTauToHadBuilder(const edm::ParameterSet& cfg)
    : NSVfitTauDecayBuilderBase(cfg)
  {}
  
  void beginJob(NSVfitAlgorithmBase* algorithm) 
  {
    NSVfitTauDecayBuilderBase::beginJob(algorithm);

    idxFitParameter_thetaVMrho_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_theta_rho, true);
    idxFitParameter_mass2VMrho_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_mass2_rho, true);
    idxFitParameter_thetaVMa1_  = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_theta_a1,  true);
    idxFitParameter_thetaVMa1r_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_theta_a1r, true);
    idxFitParameter_phiVMa1r_   = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_phi_a1r,   true);
    idxFitParameter_mass2VMa1_  = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTauVM_mass2_a1,  true);

    //if ( verbosity_ ) print(std::cout);
  }

  NSVfitSingleParticleHypothesisBase* build(const NSVfitTauDecayBuilderBase::inputParticleMap& inputParticles) const 
  {
    inputParticleMap::const_iterator particlePtr = inputParticles.find(prodParticleLabel_);
    assert(particlePtr != inputParticles.end());

    NSVfitTauToHadHypothesis* hypothesis = new NSVfitTauToHadHypothesis(particlePtr->second, prodParticleLabel_, barcodeCounter_);
    ++barcodeCounter_;

    NSVfitTauDecayBuilderBase::initialize(hypothesis, particlePtr->second.get());

    // Three prong case: check if we can fit a reconstructed vertex.
    const std::vector<reco::TrackBaseRef>& tracks = hypothesis->tracks();

    if ( tracks.size() > 1 ) {
      KalmanVertexFitter kvf(true);
      // Get the associated transient tracks
      std::vector<reco::TransientTrack> transTracks = trackService_->transientTracks(tracks.begin(), tracks.end());
      TransientVertex vertex = kvf.vertex(transTracks);
      hypothesis->fittedVertex_ = vertex;
    }

    return hypothesis;
  }

  void applyFitParameter(NSVfitSingleParticleHypothesisBase* hypothesis, double* param) const
  {
    NSVfitTauDecayBuilderBase::applyFitParameter(hypothesis, param);

    NSVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<NSVfitTauToHadHypothesis*>(hypothesis);
    assert(hypothesis_T);

    applyOptionalFitParameter(param, idxFitParameter_thetaVMrho_, hypothesis_T->decay_angle_VMrho_);
    applyOptionalFitParameter(param, idxFitParameter_mass2VMrho_, hypothesis_T->mass2_VMrho_);
    applyOptionalFitParameter(param, idxFitParameter_thetaVMa1_,  hypothesis_T->decay_angle_VMa1_);
    applyOptionalFitParameter(param, idxFitParameter_thetaVMa1r_, hypothesis_T->decay_angle_VMa1r_theta_);
    applyOptionalFitParameter(param, idxFitParameter_phiVMa1r_,   hypothesis_T->decay_angle_VMa1r_phi_);
    applyOptionalFitParameter(param, idxFitParameter_mass2VMa1_,  hypothesis_T->mass2_VMa1_);
  }

  // The two neutrion system in a leptonic decay can have mass.
  bool nuSystemIsMassless() const { return true; }

  virtual int getDecayMode(const reco::Candidate* candidate) const
  {
    const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(candidate);
    assert(tauPtr);
    return tauPtr->decayMode();
  }

  virtual std::vector<reco::TrackBaseRef> extractTracks(const reco::Candidate* candidate) const
  {
    const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(candidate);
    assert(tauPtr);
    return trackExtractor_(*tauPtr);
  }

  void print(std::ostream& stream) const 
  {
    NSVfitTauDecayBuilderBase::print(stream);
    stream << " idxFitParameter_thetaVMrho = " << idxFitParameter_thetaVMrho_ << std::endl;
    stream << " idxFitParameter_mass2VMrho = " << idxFitParameter_mass2VMrho_ << std::endl;
    stream << " idxFitParameter_thetaVMa1 = " << idxFitParameter_thetaVMa1_ << std::endl;
    stream << " idxFitParameter_thetaVMa1r = " << idxFitParameter_thetaVMa1r_ << std::endl;
    stream << " idxFitParameter_phiVMa1r = " << idxFitParameter_phiVMa1r_ << std::endl;
    stream << " idxFitParameter_mass2VMa1 = " << idxFitParameter_mass2VMa1_ << std::endl;
  }

  private:
    edm::Service<NSVfitTrackService> trackService_;
    SVfitLegTrackExtractor<pat::Tau> trackExtractor_;

    int idxFitParameter_thetaVMrho_;
    int idxFitParameter_mass2VMrho_;
    int idxFitParameter_thetaVMa1_;
    int idxFitParameter_thetaVMa1r_;
    int idxFitParameter_phiVMa1r_;
    int idxFitParameter_mass2VMa1_;
};

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitSingleParticleBuilderPluginFactory, NSVfitTauToHadBuilder, "NSVfitTauToHadBuilder");


