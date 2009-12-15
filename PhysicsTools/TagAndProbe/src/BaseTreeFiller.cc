#include "PhysicsTools/TagAndProbe/interface/BaseTreeFiller.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

tnp::ProbeVariable::~ProbeVariable() {}

tnp::ProbeFlag::~ProbeFlag() {}

void tnp::ProbeFlag::init(const edm::Event &iEvent) const {
    if (external_) {
        edm::Handle<edm::View<reco::Candidate> > view;
        iEvent.getByLabel(src_, view);
        passingProbes_.clear();
        for (size_t i = 0, n = view->size(); i < n; ++i) passingProbes_.push_back(view->refAt(i));
    }

}

void tnp::ProbeFlag::fill(const reco::CandidateBaseRef &probe) const {
    if (external_) {
        value_ = (std::find(passingProbes_.begin(), passingProbes_.end(), probe) != passingProbes_.end());
    } else {
        value_ = bool(cut_(*probe));
    }
}

tnp::BaseTreeFiller::BaseTreeFiller(const char *name, const edm::ParameterSet iConfig) {
    // make trees as requested
    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>(name,name);

    // set up variables
    edm::ParameterSet variables = iConfig.getParameter<edm::ParameterSet>("variables");
    //.. the ones that are strings
    std::vector<std::string> stringVars = variables.getParameterNamesForType<std::string>();
    for (std::vector<std::string>::const_iterator it = stringVars.begin(), ed = stringVars.end(); it != ed; ++it) {
        vars_.push_back(tnp::ProbeVariable(*it, variables.getParameter<std::string>(*it)));
    }
    //.. the ones that are InputTags
    std::vector<std::string> inputTagVars = variables.getParameterNamesForType<edm::InputTag>();
    for (std::vector<std::string>::const_iterator it = inputTagVars.begin(), ed = inputTagVars.end(); it != ed; ++it) {
        vars_.push_back(tnp::ProbeVariable(*it, variables.getParameter<edm::InputTag>(*it)));
    }
 
    // set up flags
    edm::ParameterSet flags = iConfig.getParameter<edm::ParameterSet>("flags");
    //.. the ones that are strings
    std::vector<std::string> stringFlags = flags.getParameterNamesForType<std::string>();
    for (std::vector<std::string>::const_iterator it = stringFlags.begin(), ed = stringFlags.end(); it != ed; ++it) {
        flags_.push_back(tnp::ProbeFlag(*it, flags.getParameter<std::string>(*it)));
    }
    //.. the ones that are InputTags
    std::vector<std::string> inputTagFlags = flags.getParameterNamesForType<edm::InputTag>();
    for (std::vector<std::string>::const_iterator it = inputTagFlags.begin(), ed = inputTagFlags.end(); it != ed; ++it) {
        flags_.push_back(tnp::ProbeFlag(*it, flags.getParameter<edm::InputTag>(*it)));
    }

    // then make all the variables in the trees
    for (std::vector<tnp::ProbeVariable>::iterator it = vars_.begin(), ed = vars_.end(); it != ed; ++it) {
        tree_->Branch(it->name().c_str(), it->address(), (it->name()+"/F").c_str());
    }
    for (std::vector<tnp::ProbeFlag>::iterator it = flags_.begin(), ed = flags_.end(); it != ed; ++it) {
        tree_->Branch(it->name().c_str(), it->address(), (it->name()+"/B").c_str());
    }
    
    // set up weights, if needed
    if (iConfig.existsAs<double>("eventWeight")) { 
        weightMode_ = Fixed;
        weight_ = iConfig.getParameter<double>("eventWeight");
    } else if (iConfig.existsAs<edm::InputTag>("eventWeight")) { 
        weightMode_ = External;
        weightSrc_ = iConfig.getParameter<edm::InputTag>("eventWeight");
    } else {
        weightMode_ = None;
    }
    if (weightMode_ != None) {
        tree_->Branch("weight", &weight_, "weight/F");
    }
}

tnp::BaseTreeFiller::~BaseTreeFiller() { }

void tnp::BaseTreeFiller::init(const edm::Event &iEvent) const {
    for (std::vector<tnp::ProbeVariable>::const_iterator it = vars_.begin(), ed = vars_.end(); it != ed; ++it) {
        it->init(iEvent);
    }
    for (std::vector<tnp::ProbeFlag>::const_iterator it = flags_.begin(), ed = flags_.end(); it != ed; ++it) {
        it->init(iEvent);
    }
    if (weightMode_ == External) {
        edm::Handle<double> weight;
        iEvent.getByLabel(weightSrc_, weight);
        weight_ = *weight;
    }
}

void tnp::BaseTreeFiller::fill(const reco::CandidateBaseRef &probe) const {
    for (std::vector<tnp::ProbeVariable>::const_iterator it = vars_.begin(), ed = vars_.end(); it != ed; ++it) {
        it->fill(probe);
    }
    for (std::vector<tnp::ProbeFlag>::const_iterator it = flags_.begin(), ed = flags_.end(); it != ed; ++it) {
        it->fill(probe);
    }
    tree_->Fill();
}

