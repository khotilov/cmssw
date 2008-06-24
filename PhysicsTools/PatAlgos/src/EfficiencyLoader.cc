#include "PhysicsTools/PatAlgos/interface/EfficiencyLoader.h"

#include <algorithm>

using pat::helper::EfficiencyLoader;

EfficiencyLoader::EfficiencyLoader(const edm::ParameterSet &iConfig) 
{
    // Get the names (sorted)
    names_ = iConfig.getParameterNamesForType<edm::InputTag>();
    std::sort(names_.begin(), names_.end());
    
    // get the InputTags
    for (std::vector<std::string>::const_iterator it = names_.begin(), ed = names_.end(); it != ed; ++it) {
        tags_.push_back( iConfig.getParameter<edm::InputTag>(*it) );
    }

    // prepare the Handles
    handles_.resize(names_.size());
}

void
EfficiencyLoader::newEvent(const edm::Event &iEvent) const {
    for (size_t i = 0, n = names_.size(); i < n; ++i) {
        iEvent.getByLabel(tags_[i], handles_[i]);
    }    
}
