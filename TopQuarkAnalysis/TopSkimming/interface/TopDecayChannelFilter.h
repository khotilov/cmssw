#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

template <typename S>
class TopDecayChannelFilter : public edm::EDFilter {
 public:
  TopDecayChannelFilter(const edm::ParameterSet&);
  ~TopDecayChannelFilter();
  
 private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);  
  edm::InputTag src_;    
  S sel_;  
};

template<typename S>
TopDecayChannelFilter<S>::TopDecayChannelFilter(const edm::ParameterSet& cfg):
  src_( cfg.template getParameter<edm::InputTag>( "src" ) ),
  sel_( cfg )
{ }

template<typename S>
TopDecayChannelFilter<S>::~TopDecayChannelFilter()
{ }

template<typename S>
bool
TopDecayChannelFilter<S>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection> parts;
  iEvent.getByLabel( src_,parts );
 
  return sel_( *parts );
}
