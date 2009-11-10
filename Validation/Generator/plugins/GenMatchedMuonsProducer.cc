#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "UserCode/TopMCValidation/plugins/GenMatchedCandsProducerBase.h"
typedef GenMatchedCandsProducerBase<reco::Muon, reco::GenParticle> GenMatchedMuonsProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMatchedMuonsProducer);
