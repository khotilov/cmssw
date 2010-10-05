/*
 * RecoTauPiZeroProducer
 *
 * Author: Evan K. Friis, UC Davis
 *
 * Associates reconstructed PiZeros to PFJets.  The PiZeros are built using one
 * or more RecoTauBuilder plugins.  Any overlaps (PiZeros sharing constituents)
 * are removed, with the best PiZero candidates taken.  The 'best' are defined
 * via the input list of RecoTauPiZeroQualityPlugins, which form a
 * lexicograpical ranking.
 *
 * $Id $
 */

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/foreach.hpp>
#include <algorithm>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoTauTag/RecoTau/interface/RecoTauPiZeroPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCleaningTools.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/JetPiZeroAssociation.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

class RecoTauPiZeroProducer : public edm::EDProducer {
  public:
    typedef reco::tau::RecoTauPiZeroBuilderPlugin Builder;
    typedef reco::tau::RecoTauPiZeroQualityPlugin Ranker;

    explicit RecoTauPiZeroProducer(const edm::ParameterSet& pset);
    ~RecoTauPiZeroProducer() {}
    void produce(edm::Event& evt, const edm::EventSetup& es);
    void print(const std::vector<reco::RecoTauPiZero>& piZeros,
               std::ostream& out);

  private:
    typedef boost::ptr_vector<Builder> builderList;
    typedef boost::ptr_vector<Ranker> rankerList;
    typedef reco::tau::RecoTauLexicographicalRanking<rankerList,
            reco::RecoTauPiZero> PiZeroPredicate;
    edm::InputTag jetSrc_;
    builderList builders_;
    rankerList rankers_;
    std::auto_ptr<PiZeroPredicate> predicate_;
};

RecoTauPiZeroProducer::RecoTauPiZeroProducer(const edm::ParameterSet& pset) {
  jetSrc_ = pset.getParameter<edm::InputTag>("jetSrc");

  typedef std::vector<edm::ParameterSet> VPSet;
  // Get each of our PiZero builders
  const VPSet& builders = pset.getParameter<VPSet>("builders");

  for (VPSet::const_iterator builderPSet = builders.begin();
      builderPSet != builders.end(); ++builderPSet) {
    // Get plugin name
    const std::string& pluginType =
      builderPSet->getParameter<std::string>("plugin");
    // Build the plugin
    builders_.push_back(RecoTauPiZeroBuilderPluginFactory::get()->create(
          pluginType, *builderPSet));
  }

  // Get each of our quality rankers
  const VPSet& rankers = pset.getParameter<VPSet>("ranking");
  for (VPSet::const_iterator rankerPSet = rankers.begin();
      rankerPSet != rankers.end(); ++rankerPSet) {
    const std::string& pluginType =
      rankerPSet->getParameter<std::string>("plugin");
    rankers_.push_back(RecoTauPiZeroQualityPluginFactory::get()->create(
          pluginType, *rankerPSet));
  }

  // Build the sorting predicate
  predicate_ = std::auto_ptr<PiZeroPredicate>(new PiZeroPredicate(rankers_));

  produces<reco::JetPiZeroAssociation>();
}

void RecoTauPiZeroProducer::produce(edm::Event& evt,
                                    const edm::EventSetup& es) {
  edm::Handle<reco::PFJetCollection> pfJets;
  evt.getByLabel(jetSrc_, pfJets);

  // Make our association
  std::auto_ptr<reco::JetPiZeroAssociation> association(
      new reco::JetPiZeroAssociation(reco::PFJetRefProd(pfJets)));

  size_t iJet = 0;
  // Loop over our jets
  for (reco::PFJetCollection::const_iterator jet = pfJets->begin();
      jet != pfJets->end(); ++jet, ++iJet) {
    size_t numberOfGammas = reco::tau::pfCandidates(
        *jet, reco::PFCandidate::gamma).size();
    typedef boost::ptr_vector<reco::RecoTauPiZero> PiZeroVector;
    typedef boost::ptr_list<reco::RecoTauPiZero> PiZeroList;
    // Build our global list of RecoTauPiZero
    PiZeroList dirtyPiZeros;

    // Compute the pi zeros from this jet for all the desired algorithms
    BOOST_FOREACH(const Builder& builder, builders_) {
      PiZeroVector result(builder(*jet));
      dirtyPiZeros.transfer(dirtyPiZeros.end(), result);
    }
    // Rank the candidates according to our quality plugins
    dirtyPiZeros.sort(*predicate_);

    // Keep track of the photons in the clean collection
    std::vector<reco::RecoTauPiZero> cleanPiZeros;
    std::set<reco::CandidatePtr> photonsInCleanCollection;
    while (dirtyPiZeros.size() &&
           numberOfGammas > photonsInCleanCollection.size()) {
      // Pull our candidate pi zero from the front of the list
      std::auto_ptr<reco::RecoTauPiZero> toAdd(
          dirtyPiZeros.pop_front().release());
      // Find the sub-gammas that are not already in the cleaned collection
      std::vector<reco::CandidatePtr> uniqueGammas;
      std::set_difference(toAdd->daughterPtrVector().begin(),
                          toAdd->daughterPtrVector().end(),
                          photonsInCleanCollection.begin(),
                          photonsInCleanCollection.end(),
                          std::back_inserter(uniqueGammas));
      // If the pi zero has no unique gammas, discard it.  Note toAdd is deleted
      // when it goes out of scope.
      if (!uniqueGammas.size()) {
        continue;
      } else if (uniqueGammas.size() == toAdd->daughterPtrVector().size()) {
        // Check if it is composed entirely of unique gammas.  In this case
        // immediately add it to the clean collection.
        photonsInCleanCollection.insert(toAdd->daughterPtrVector().begin(),
                                        toAdd->daughterPtrVector().end());
        cleanPiZeros.push_back(*toAdd);
      } else {
        // Otherwise update the pizero that contains only the unique gammas and
        // add it back into the sorted list of dirty PiZeros
        toAdd->clearDaughters();
        // Add each of the unique daughters back to the pizero
        BOOST_FOREACH(const reco::CandidatePtr& gamma, uniqueGammas) {
          toAdd->addDaughter(gamma);
        }
        // Update the four vector
        AddFourMomenta p4Builder_;
        p4Builder_.set(*toAdd);
        // Put this pi zero back into the collection of sorted dirty pizeros
        PiZeroList::iterator insertionPoint = std::lower_bound(
            dirtyPiZeros.begin(), dirtyPiZeros.end(), *toAdd, *predicate_);
        dirtyPiZeros.insert(insertionPoint, toAdd);
      }
    }
    // Add to association
    association->setValue(iJet, cleanPiZeros);
  }
  evt.put(association);
}

// Print some helpful information
void RecoTauPiZeroProducer::print(
    const std::vector<reco::RecoTauPiZero>& piZeros, std::ostream& out) {
  const unsigned int width = 25;
  BOOST_FOREACH(const reco::RecoTauPiZero& piZero, piZeros) {
    out << piZero;
    out << "* Rankers:" << std::endl;
    for (rankerList::const_iterator ranker = rankers_.begin();
        ranker != rankers_.end(); ++ranker) {
      out << "* " << std::setiosflags(std::ios::left)
        << std::setw(width) << ranker->name()
        << " " << std::resetiosflags(std::ios::left)
        << std::setprecision(3) << (*ranker)(piZero);
      out << std::endl;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoTauPiZeroProducer);

