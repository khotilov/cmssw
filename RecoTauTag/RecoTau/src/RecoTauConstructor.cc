#include "RecoTauTag/RecoTau/interface/RecoTauConstructor.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

namespace reco {
  namespace tau {

    RecoTauConstructor::RecoTauConstructor(const PFJetRef& jet, 
        const edm::Handle<PFCandidateCollection>& pfCands,
        bool copyGammasFromPiZeros):pfCands_(pfCands)
    {
      copyGammas_ = copyGammasFromPiZeros;
      // Initialize our Accessors
      collections_[std::make_pair(kSignal, kChargedHadron)] = &tau_.selectedSignalPFChargedHadrCands_;
      collections_[std::make_pair(kSignal, kGamma)] = &tau_.selectedSignalPFGammaCands_;
      collections_[std::make_pair(kSignal, kNeutralHadron)] = &tau_.selectedSignalPFNeutrHadrCands_;
      collections_[std::make_pair(kSignal, kAll)] = &tau_.selectedSignalPFCands_;

      collections_[std::make_pair(kIsolation, kChargedHadron)] = &tau_.selectedIsolationPFChargedHadrCands_;
      collections_[std::make_pair(kIsolation, kGamma)] = &tau_.selectedIsolationPFGammaCands_;
      collections_[std::make_pair(kIsolation, kNeutralHadron)] = &tau_.selectedIsolationPFNeutrHadrCands_;
      collections_[std::make_pair(kIsolation, kAll)] = &tau_.selectedIsolationPFCands_;

      tau_.setjetRef(jet);
    }

    void RecoTauConstructor::addPFCand(Region region, ParticleType type, const PFCandidateRef& ref)
    {
      getCollection(region, type)->push_back(ref);
      // Add to global collection
      getCollection(region, kAll)->push_back(ref);
    }

    void RecoTauConstructor::reserve(Region region, ParticleType type, size_t size)
    {
      getCollection(region, type)->reserve(size);
      // Reserve global collection as well
      getCollection(region, kAll)->reserve(
          getCollection(region, kAll)->size() + size);
    }

    void RecoTauConstructor::reservePiZero(Region region, size_t size)
    {
      if(region == kSignal)
      {
        tau_.signalPiZeroCandidates_.reserve(size);
        // If we are building the gammas with the pizeros, resize that 
        // vector as well
        if(copyGammas_)
          reserve(kSignal, kGamma, 2*size);
      }
      else
      {
        tau_.isolationPiZeroCandidates_.reserve(size);
        if(copyGammas_)
          reserve(kIsolation, kGamma, 2*size);
      }
    }

    void RecoTauConstructor::addPiZero(Region region, const RecoTauPiZero& piZero)
    {
      if(region == kSignal) {
        tau_.signalPiZeroCandidates_.push_back(piZero);
        // Copy the daughter gammas into the gamma collection if desired
        if(copyGammas_) {
          addPFCands(kSignal, kGamma, piZero.daughterPtrVector().begin(), 
              piZero.daughterPtrVector().end());
        }
      }
      else
      {
        tau_.isolationPiZeroCandidates_.push_back(piZero);
        if(copyGammas_) {
          addPFCands(kIsolation, kGamma, piZero.daughterPtrVector().begin(), 
              piZero.daughterPtrVector().end());
        }
      }
    }

    PFCandidateRefVector * RecoTauConstructor::getCollection(Region region, ParticleType type) 
    {
      return collections_[std::make_pair(region, type)];
    }

    // Trivial converter needed for polymorphism
    PFCandidateRef RecoTauConstructor::convertToRef(const PFCandidateRef& pfRef) const
    {
      return pfRef;
    }

    // Convert from a Ptr to a Ref
    PFCandidateRef RecoTauConstructor::convertToRef(const PFCandidatePtr& pfPtr) const
    {
      assert(pfPtr.id() == pfCands_.id());
      return PFCandidateRef(pfCands_, pfPtr.key());
    }

    // Convert from a CandidatePtr to a Ref
    PFCandidateRef RecoTauConstructor::convertToRef(const CandidatePtr& candPtr) const
    {
      assert(candPtr.id() == pfCands_.id());
      return PFCandidateRef(pfCands_, candPtr.key());
    }

    const PFTau& RecoTauConstructor::get(bool setupLeadingObjects) {
      // Setup all the important member variables of the tau

      // Set charge of tau
      tau_.setCharge(
          sumPFCandCharge(
            getCollection(kSignal, kChargedHadron)->begin(),
            getCollection(kSignal, kChargedHadron)->end()
            )
          );

      // Set PDG id
      tau_.setPdgId(tau_.charge() > 0 ? 15 : -15);

      // Set P4
      tau_.setP4(
          sumPFCandP4(
            getCollection(kSignal, kAll)->begin(),
            getCollection(kSignal, kAll)->end()
            )
          );

      // Set charged isolation quantities
      tau_.setisolationPFChargedHadrCandsPtSum(
          sumPFCandPt(
            getCollection(kIsolation, kChargedHadron)->begin(),
            getCollection(kIsolation, kChargedHadron)->end()
            )
          );

      // Set gamma isolation quantities
      tau_.setisolationPFGammaCandsEtSum(
          sumPFCandPt(
            getCollection(kIsolation, kGamma)->begin(),
            getCollection(kIsolation, kGamma)->end()
            )
          );

      if(setupLeadingObjects)
      {
        typedef PFCandidateRefVector::const_iterator Iter;
        // Find the highest PT object in the signal cone
        Iter leadingCand = leadPFCand(
            getCollection(kSignal, kAll)->begin(),
            getCollection(kSignal, kAll)->end()
            );

        if(leadingCand != getCollection(kSignal, kAll)->end())
          tau_.setleadPFCand(*leadingCand);

        // Hardest charged object in signal cone
        Iter leadingChargedCand = leadPFCand(
            getCollection(kSignal, kChargedHadron)->begin(),
            getCollection(kSignal, kChargedHadron)->end()
            );

        if(leadingChargedCand != getCollection(kSignal, kChargedHadron)->end())
          tau_.setleadPFChargedHadrCand(*leadingChargedCand);

        // Hardest gamma object in signal cone
        Iter leadingGammaCand = leadPFCand(
            getCollection(kSignal, kGamma)->begin(),
            getCollection(kSignal, kGamma)->end()
            );

        if(leadingGammaCand != getCollection(kSignal, kGamma)->end())
          tau_.setleadPFNeutralCand(*leadingGammaCand);
      }
      return tau_;
    }
  }
}






