/*
 * RecoTauBuilderConePlugin
 *
 * Build a PFTau using cones defined in DeltaR.
 *
 * Original Authors: Ludovic Houchu, Simone Gennai
 * Modifications: Evan K. Friis
 *
 * $Id $
 */

#include <vector>
#include <algorithm>

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/RecoTauConstructor.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

namespace reco { namespace tau {

namespace {
// Predicate class that tests if a candidate lies within some deltaR (min,
// max) about a supplied axis
template<class CandType>
class DeltaRFilter {
  public:
    DeltaRFilter(const reco::Candidate::LorentzVector& axis,
                 double min, double max): axis_(axis), min_(min), max_(max) {}
    bool operator()(const CandType& b) {
      double deltaR = reco::deltaR<reco::Candidate::LorentzVector>
          (axis_, b.p4());
      return(deltaR >= min_ && deltaR < max_);
    }
  private:
    const reco::Candidate::LorentzVector& axis_;
    const double min_;
    const double max_;
};

// Wrapper around DeltaRFilter to support reference types like Ptr<>
template<class CandType>
class DeltaRPtrFitler {
  public:
    DeltaRPtrFitler(const reco::Candidate::LorentzVector& axis,
                    double min, double max): filter_(axis, min, max) {}
    bool operator()(const CandType& b) { return filter_(*b); }
  private:
    DeltaRFilter<typename CandType::value_type> filter_;
};
}


class RecoTauBuilderConePlugin : public RecoTauBuilderPlugin {
  public:
    explicit RecoTauBuilderConePlugin(const edm::ParameterSet& pset);
    ~RecoTauBuilderConePlugin() {}
    // Build a tau from a jet
    std::vector<reco::PFTau> operator()(
        const reco::PFJetRef& jet,
        const std::vector<RecoTauPiZero>& piZeros) const;
  private:
    double leadObjecPtThreshold_;
    /* String function to extract values from PFJets */
    typedef StringObjectFunction<reco::PFJet> JetFunc;
    // Cone defintions
    JetFunc matchingCone_;
    JetFunc signalConeChargedHadrons_;
    JetFunc isoConeChargedHadrons_;
    JetFunc signalConePiZeros_;
    JetFunc isoConePiZeros_;
    JetFunc signalConeNeutralHadrons_;
    JetFunc isoConeNeutralHadrons_;
};

// ctor - initialize all of our variables
RecoTauBuilderConePlugin::RecoTauBuilderConePlugin(
    const edm::ParameterSet& pset):RecoTauBuilderPlugin(pset),
    leadObjecPtThreshold_(pset.getParameter<double>("leadObjectPt")),
    matchingCone_(pset.getParameter<std::string>("matchingCone")),
    signalConeChargedHadrons_(pset.getParameter<std::string>(
            "signalConeChargedHadrons")),
    isoConeChargedHadrons_(
        pset.getParameter<std::string>("isoConeChargedHadrons")),
    signalConePiZeros_(
        pset.getParameter<std::string>("signalConePiZeros")),
    isoConePiZeros_(
        pset.getParameter<std::string>("isoConePiZeros")),
    signalConeNeutralHadrons_(
        pset.getParameter<std::string>("signalConeNeutralHadrons")),
    isoConeNeutralHadrons_(
        pset.getParameter<std::string>("isoConeNeutralHadrons"))
{}

std::vector<reco::PFTau> RecoTauBuilderConePlugin::operator()(
    const reco::PFJetRef& jet,
    const std::vector<RecoTauPiZero>& piZeros) const {
  /* Define our filters */
  typedef DeltaRPtrFitler<PFCandidatePtr> PFCandPtrDRFilter;
  typedef boost::filter_iterator< PFCandPtrDRFilter,
          std::vector<PFCandidatePtr>::const_iterator> PFCandPtrDRFilterIter;

  typedef DeltaRFilter<RecoTauPiZero> PiZeroDRFilter;
  typedef boost::filter_iterator< PiZeroDRFilter,
          std::vector<RecoTauPiZero>::const_iterator> PiZeroDRFilterIter;

  // Define output.  We only produce one tau per jet for the cone algo.
  std::vector<reco::PFTau> output;

  // Our tau builder - the true indicates to automatically copy gamma candidates
  // from the pizeros.
  RecoTauConstructor tau(jet, getPFCands(), true);

  typedef std::vector<PFCandidatePtr> PFCandPtrs;

  // Get the PF Charged hadrons
  PFCandPtrs pfchs = pfCandidates(*jet, reco::PFCandidate::h);
  // Get the PF gammas
  PFCandPtrs pfGammas = pfCandidates(*jet, reco::PFCandidate::gamma);
  // Neutral hadrons
  PFCandPtrs pfnhs = pfCandidates(*jet, reco::PFCandidate::h0);

  /***********************************************
   ******     Lead Candidate Finding    **********
   ***********************************************/

  // Define our matching cone and filters
  double matchingCone = matchingCone_(*jet);
  PFCandPtrDRFilter matchingConeFilter(jet->p4(), 0, matchingCone);

  // Find the maximum PFCharged hadron in the matching cone.  The call to
  // PFCandidates always a sorted list, so we can just take the first if it
  // if it exists.
  PFCandidatePtr leadPFCH;
  PFCandPtrs::iterator leadPFCH_iter =
      std::find_if(pfchs.begin(), pfchs.end(), matchingConeFilter);

  if (leadPFCH_iter != pfchs.end()) {
    leadPFCH = *leadPFCH_iter;
    // Set leading candidate
    tau.setleadPFChargedHadrCand(leadPFCH);
  }

  // Find the leading neutral candidate
  PFCandidatePtr leadPFGamma;
  PFCandPtrs::iterator leadPFGamma_iter =
      std::find_if(pfGammas.begin(), pfGammas.end(), matchingConeFilter);

  if (leadPFGamma_iter != pfGammas.end()) {
    leadPFGamma = *leadPFGamma_iter;
    // Set leading neutral candidate
    tau.setleadPFNeutralCand(leadPFGamma);
  }

  PFCandidatePtr leadPFCand;
  // Always use the leadPFCH if it is above our threshold
  if (leadPFCH.isNonnull() && leadPFCH->pt() > leadObjecPtThreshold_) {
    leadPFCand = leadPFCH;
  } else if (leadPFGamma.isNonnull() &&
             leadPFGamma->pt() > leadObjecPtThreshold_) {
    // Otherwise use the lead Gamma if it is above threshold
    leadPFCand = leadPFGamma;
  } else {
    // If both are too low PT, just take the charged one
    leadPFCand = leadPFCH;
  }

  tau.setleadPFCand(leadPFCand);

  // If there is no leading charged candidate at all, return a blank tau
  if (leadPFCH.isNull()) {
    // Set the P4 as that of the jet.
    PFTau tau(0, jet->p4());
    tau.setjetRef(jet);
    output.push_back(tau);
    return output;
  }

  // Our cone axis is defined about the lead charged hadron
  reco::Candidate::LorentzVector coneAxis = leadPFCH->p4();

  /***********************************************
   ******     Cone Construction         **********
   ***********************************************/

  // Define the signal and isolation cone sizes for this jet and build filters
  // to select elements in the given DeltaR regions

  PFCandPtrDRFilter signalConePFCHFilter(
      coneAxis, -0.1, signalConeChargedHadrons_(*jet));
  PFCandPtrDRFilter signalConePFNHFilter(
      coneAxis, -0.1, signalConeNeutralHadrons_(*jet));
  PiZeroDRFilter signalConePiZeroFilter(
      coneAxis, -0.1, signalConePiZeros_(*jet));

  PFCandPtrDRFilter isoConePFCHFilter(
      coneAxis, signalConeChargedHadrons_(*jet), isoConeChargedHadrons_(*jet));
  PFCandPtrDRFilter isoConePFNHFilter(
      coneAxis, signalConeNeutralHadrons_(*jet), isoConeNeutralHadrons_(*jet));
  PiZeroDRFilter isoConePiZeroFilter(
      coneAxis, signalConePiZeros_(*jet), isoConePiZeros_(*jet));

  // Build signal charged hadrons
  tau.addPFCands(RecoTauConstructor::kSignal,
                 RecoTauConstructor::kChargedHadron,
                 PFCandPtrDRFilterIter(signalConePFCHFilter, pfchs.begin(),
                                       pfchs.end()),
                 PFCandPtrDRFilterIter(signalConePFCHFilter, pfchs.end(),
                                       pfchs.end()));

  // Build signal neutral hadrons
  tau.addPFCands(RecoTauConstructor::kSignal,
                 RecoTauConstructor::kNeutralHadron,
                 PFCandPtrDRFilterIter(signalConePFNHFilter, pfnhs.begin(),
                                       pfnhs.end()),
                 PFCandPtrDRFilterIter(signalConePFNHFilter, pfnhs.end(),
                                       pfnhs.end()));

  // Build signal PiZeros
  tau.addPiZeros(RecoTauConstructor::kSignal,
                 PiZeroDRFilterIter(signalConePiZeroFilter,
                                    piZeros.begin(), piZeros.end()),
                 PiZeroDRFilterIter(signalConePiZeroFilter,
                                    piZeros.end(), piZeros.end()));

  // Build isolation charged hadrons
  tau.addPFCands(RecoTauConstructor::kIsolation,
                 RecoTauConstructor::kChargedHadron,
                 PFCandPtrDRFilterIter(isoConePFCHFilter, pfchs.begin(),
                                       pfchs.end()),
                 PFCandPtrDRFilterIter(isoConePFCHFilter, pfchs.end(),
                                       pfchs.end()));

  // Build isolation neutral hadrons
  tau.addPFCands(RecoTauConstructor::kIsolation,
                 RecoTauConstructor::kNeutralHadron,
                 PFCandPtrDRFilterIter(isoConePFNHFilter, pfnhs.begin(),
                                       pfnhs.end()),
                 PFCandPtrDRFilterIter(isoConePFNHFilter, pfnhs.end(),
                                       pfnhs.end()));

  // Build isolation PiZeros
  tau.addPiZeros(RecoTauConstructor::kIsolation,
                 PiZeroDRFilterIter(isoConePiZeroFilter, piZeros.begin(),
                                    piZeros.end()),
                 PiZeroDRFilterIter(isoConePiZeroFilter, piZeros.end(),
                                    piZeros.end()));

  // Put our built tau in the output - 'false' indicates don't build the
  // leading candidtes, we already did that explicitly above.

  output.push_back(tau.get(false));
  return output;
}
}}  // end namespace reco::tauk

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(RecoTauBuilderPluginFactory,
                  reco::tau::RecoTauBuilderConePlugin,
                  "RecoTauBuilderConePlugin");
