#ifndef RecoTauTag_RecoTau_RecoTauMVAHelper_h
#define RecoTauTag_RecoTau_RecoTauMVAHelper_h

/*
 * RecoTauMVAHelper
 *
 * Author: Evan K. Friis (UC Davis)
 *
 * Manages DB retrieval and applictionat of MVAComputers to reco::PFTaus.
 *
 * Based on code by Christophe Saoute in
 *  - PhysicsTools/MVAComputer/interface/MVAModuleHelper.h
 */

#include <boost/ptr_container/ptr_map.hpp>
#include <string>

#include "RecoTauTag/RecoTau/interface/RecoTauDiscriminantPlugins.h"

#include "PhysicsTools/MVAComputer/interface/MVAComputerCache.h"
#include "PhysicsTools/MVAComputer/interface/Variable.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

namespace reco { namespace tau {

class RecoTauMVAHelper {
  public:
    explicit RecoTauMVAHelper(const std::string &name,
                              const std::string eslabel="");
    ~RecoTauMVAHelper() {}

    // Setup event information and retrive MVA from DB
    void setEvent(const edm::Event& evt, const edm::EventSetup &es);

    // Apply MVA to tau and return result
    double operator()(const reco::PFTauRef &tau) const;

    // Add a training event of type <target> with given weight
    void train(const reco::PFTauRef &tau, bool target,
               double weight = 1.0) const;
  private:
    // Name of computer in the DB record
    std::string name_;
    // Name of label for event setup record
    std::string eslabel_;
    PhysicsTools::MVAComputerCache computer_;
    // Map our discriminant plugins to their "MVA name"
    typedef boost::ptr_map<PhysicsTools::AtomicId, RecoTauDiscriminantPlugin>
        PluginMap;
    PluginMap plugins_;
    // Helper function to load relevant plugins for a MVA computer
    void loadDiscriminantPlugins(
        const PhysicsTools::Calibration::MVAComputer &computer);
    // Load the plugin values for this tau
    void fillValues(const reco::PFTauRef& tau) const;
    mutable PhysicsTools::Variable::ValueList values_;
};
}}  // end reco::tau
#endif
