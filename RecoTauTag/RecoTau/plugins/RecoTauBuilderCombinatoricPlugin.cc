#include <vector>

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "RecoTauTag/RecoTau/interface/CombinatoricGenerator.h"
#include "RecoTauTag/RecoTau/interface/RecoTauConstructor.h"

namespace reco { namespace tau {
class RecoTauBuilderCombinatoricPlugin : public RecoTauBuilderPlugin {
 public:
  explicit RecoTauBuilderCombinatoricPlugin(const edm::ParameterSet& pset);
  virtual ~RecoTauBuilderCombinatoricPlugin() {}
  virtual std::vector<reco::PFTau> operator()
      (const reco::PFJetRef& jet,
       const std::vector<RecoTauPiZero>& piZeros) const;
 private:
  struct decayModeInfo {
    uint32_t maxPiZeros_;
    uint32_t maxPFCHs_;
    uint32_t nCharged_;
    uint32_t nPiZeros_;
  };
  std::vector<decayModeInfo> decayModesToBuild_;
};

RecoTauBuilderCombinatoricPlugin::RecoTauBuilderCombinatoricPlugin(
    const edm::ParameterSet& pset): RecoTauBuilderPlugin(pset) {
  typedef std::vector<edm::ParameterSet> VPSet;
  const VPSet& decayModes = pset.getParameter<VPSet>("decayModes");
  for (VPSet::const_iterator dm = decayModes.begin();
       dm != decayModes.end(); ++dm) {
    decayModeInfo info;
    info.nCharged_ = dm->getParameter<uint32_t>("nCharged");
    info.nPiZeros_ = dm->getParameter<uint32_t>("nPiZeros");
    info.maxPFCHs_ = dm->getParameter<uint32_t>("maxTracks");
    info.maxPiZeros_ = dm->getParameter<uint32_t>("maxPiZeros");
    decayModesToBuild_.push_back(info);
  }
}

std::vector<reco::PFTau> RecoTauBuilderCombinatoricPlugin::operator()(
    const reco::PFJetRef& jet,
    const std::vector<RecoTauPiZero>& piZeros) const {
  typedef std::vector<PFCandidatePtr> PFCandPtrs;
  typedef std::vector<RecoTauPiZero> PiZeroList;

  std::vector<PFTau> output;

  // Get PFCHs from this jet.  They are already sorted by descending Pt
  PFCandPtrs pfchs = pfCandidates(*jet, reco::PFCandidate::h);

  // Loop over the decay modes we want to build
  for (std::vector<decayModeInfo>::const_iterator
       decayMode = decayModesToBuild_.begin();
       decayMode != decayModesToBuild_.end(); ++decayMode) {
    // Find how many piZeros are in this decay mode
    size_t piZerosToBuild = decayMode->nPiZeros_;
    // Find how many tracks are in this decay mode
    size_t tracksToBuild = decayMode->nCharged_;

    // Skip decay mode if jet doesn't have the multiplicity to support it
    if (pfchs.size() < tracksToBuild || piZeros.size() < piZerosToBuild)
      continue;

    // Find the start and end of potential signal tracks
    PFCandPtrs::iterator pfch_begin = pfchs.begin();
    PFCandPtrs::iterator pfch_end =  pfchs.end();
    pfch_end = takeNElements(pfch_begin, pfch_end, decayMode->maxPFCHs_);

    // Build our track combo generator
    typedef tau::CombinatoricGenerator<PFCandPtrs> PFCombo;
    PFCombo trackCombos(pfch_begin, pfch_end, tracksToBuild);

    // Find the start and end of potential signal tracks
    PiZeroList::const_iterator piZero_begin = piZeros.begin();
    PiZeroList::const_iterator piZero_end = piZeros.end();
    piZero_end = takeNElements(piZero_begin, piZero_end,
                               decayMode->maxPiZeros_);

    // Build our piZero combo generator
    typedef tau::CombinatoricGenerator<PiZeroList> PiZeroCombo;
    PiZeroCombo piZeroCombos(piZero_begin, piZero_end, piZerosToBuild);

    /*
     * Begin combinatoric loop for this decay mode
     */

    // Loop over the different combinations of tracks
    for (PFCombo::iterator trackCombo = trackCombos.begin();
         trackCombo != trackCombos.end(); ++trackCombo) {
      // Loop over the different combinations of PiZeros
      for (PiZeroCombo::iterator piZeroCombo = piZeroCombos.begin();
           piZeroCombo != piZeroCombos.end(); ++piZeroCombo) {
        // Output tau
        RecoTauConstructor tau(jet, getPFCands(), true);
        // Reserve space in our collections
        tau.reserve(RecoTauConstructor::kSignal,
                    RecoTauConstructor::kChargedHadron, tracksToBuild);
        tau.reserve(
            RecoTauConstructor::kSignal,
            RecoTauConstructor::kGamma, 2*piZerosToBuild);  // k-factor = 2
        tau.reservePiZero(RecoTauConstructor::kSignal, piZerosToBuild);

        tau.reserve(
            RecoTauConstructor::kIsolation,
            RecoTauConstructor::kChargedHadron, pfchs.size() - tracksToBuild);
        tau.reserve(RecoTauConstructor::kIsolation,
                    RecoTauConstructor::kGamma,
                    (piZeros.size() - piZerosToBuild)*2);
        tau.reservePiZero(RecoTauConstructor::kIsolation,
                          (piZeros.size() - piZerosToBuild));

        // Set signal and isolation components for charged hadrons, after
        // converting them to a PFCandidateRefVector
        tau.addPFCands(RecoTauConstructor::kSignal,
                       RecoTauConstructor::kChargedHadron,
                       trackCombo->combo_begin(),
                       trackCombo->combo_end());
        tau.addPFCands(RecoTauConstructor::kIsolation,
                       RecoTauConstructor::kChargedHadron,
                       trackCombo->remainder_begin(),
                       trackCombo->remainder_end());

        // Add all the candidates that weren't included in the combinatoric
        // generation
        tau.addPFCands(RecoTauConstructor::kIsolation,
                       RecoTauConstructor::kChargedHadron,
                       pfch_end, pfchs.end());

        // Get PiZero constituents and add them to the tau.
        // The sub-gammas are automatically added.
        tau.addPiZeros(RecoTauConstructor::kSignal,
                       piZeroCombo->combo_begin(),
                       piZeroCombo->combo_end());
        tau.addPiZeros(RecoTauConstructor::kIsolation,
                       piZeroCombo->remainder_begin(),
                       piZeroCombo->remainder_end());
        tau.addPiZeros(RecoTauConstructor::kIsolation,
                       piZero_end, piZeros.end());

        output.push_back(tau.get(true));
      }
    }
  }
  return output;
}
}}  // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(RecoTauBuilderPluginFactory,
                  reco::tau::RecoTauBuilderCombinatoricPlugin,
                  "RecoTauBuilderCombinatoricPlugin");
