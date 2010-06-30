//
// $Id: NearbyMuonsInfo.cc,v 1.2 2010/05/27 21:43:57 gpetrucc Exp $
//

/**
  \class    pat::NearbyMuonsInfo NearbyMuonsInfo.h "PhysicsTools/PatAlgos/interface/NearbyMuonsInfo.h"
  \brief    Matcher of reconstructed objects to L1 Muons 
            
  \author   Giovanni Petrucciani
  \version  $Id: NearbyMuonsInfo.cc,v 1.2 2010/05/27 21:43:57 gpetrucc Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include <DataFormats/Candidate/interface/CompositeCandidate.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h>
#include <functional>

class NearbyMuonsInfo : public edm::EDProducer {
    public:
      explicit NearbyMuonsInfo(const edm::ParameterSet & iConfig);
      virtual ~NearbyMuonsInfo() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
      virtual void beginRun(edm::Run & iRun, const edm::EventSetup & iSetup);
    private:
        edm::InputTag src_;
        PropagateToMuon prop_;

        /// Write a ValueMap<float> in the event
        void writeValueMap(edm::Event &iEvent,
                const edm::Handle<edm::View<reco::Candidate> > & handle,
                const std::vector<float> & values,
                const std::string    & label) const ;

};

NearbyMuonsInfo::NearbyMuonsInfo(const edm::ParameterSet & iConfig) :
    src_(iConfig.getParameter<edm::InputTag>("src")),
    prop_(iConfig)
{
    produces<edm::ValueMap<float> >("dphiVtxTimesQ");
    produces<edm::ValueMap<float> >("drM2");
    produces<edm::ValueMap<float> >("dphiM2");
    produces<edm::ValueMap<float> >("distM2");
    produces<edm::ValueMap<float> >("drStaIn");
    produces<edm::ValueMap<float> >("dphiStaIn");
}

void 
NearbyMuonsInfo::beginRun(edm::Run & iRun, const edm::EventSetup & iSetup) {
    prop_.init(iSetup);
}

void 
NearbyMuonsInfo::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
    using namespace edm;
    using namespace std;

    Handle<View<reco::Candidate> > src;
    iEvent.getByLabel(src_, src);

    size_t n = src->size();
    std::vector<float> dphiVtxTimesQ(n);
    std::vector<float> drStaIn(n, -999), dphiStaIn(n, -999);
    std::vector<float> drM2(n,    -999),    dphiM2(n, -999), distM2(n, -999);
    for (size_t i = 0; i < n; ++i) {
        const reco::Candidate & ci = (*src)[i];
        if (ci.numberOfDaughters() != 2) throw cms::Exception("CorruptData") << 
            "NearbyMuonsInfo should be used on composite candidates with two daughters, this one has " << ci.numberOfDaughters() << "\n";
        const reco::Candidate &d1 = *ci.daughter(0), &d2 = *ci.daughter(1);
        dphiVtxTimesQ[i] = d1.charge() * deltaPhi(d1.phi(), d2.phi());
        const reco::RecoCandidate *mu1 = dynamic_cast<const reco::RecoCandidate *>(&*d1.masterClone());
        const reco::RecoCandidate *mu2 = dynamic_cast<const reco::RecoCandidate *>(&*d2.masterClone());
        if (mu1 == 0) throw cms::Exception("CorruptData") << "First daughter of candidate is not a ShallowClone of a reco::RecoCandidate\n";
        if (mu2 == 0) throw cms::Exception("CorruptData") << "Second daughter of candidate is not a ShallowClone of a reco::RecoCandidate\n";
        if (mu1->standAloneMuon().isNonnull()            && mu2->standAloneMuon().isNonnull() &&
            mu1->standAloneMuon().isAvailable()          && mu2->standAloneMuon().isAvailable() &&
            mu1->standAloneMuon()->extra().isAvailable() && mu2->standAloneMuon()->extra().isAvailable()  ) {
            dphiStaIn[i] = deltaPhi(mu1->standAloneMuon()->innerPosition().Phi(), mu2->standAloneMuon()->innerPosition().Phi());
            drStaIn[i]   = hypot(dphiStaIn[i], std::abs(mu1->standAloneMuon()->innerPosition().Eta() - mu2->standAloneMuon()->innerPosition().Eta()));
        }
        TrajectoryStateOnSurface prop1 = prop_.extrapolate(*mu1);
        TrajectoryStateOnSurface prop2 = prop_.extrapolate(*mu2);
        if (prop1.isValid() && prop2.isValid()) {
            dphiM2[i] = deltaPhi<float>(prop1.globalPosition().phi(), prop2.globalPosition().phi());
            drM2[i]   = hypot(dphiM2[i], std::abs<float>(prop1.globalPosition().eta() - prop2.globalPosition().eta()));
            distM2[i] = (prop1.globalPosition()-prop2.globalPosition()).mag();
        }
    }

    writeValueMap(iEvent, src, dphiVtxTimesQ,  "dphiVtxTimesQ");
    writeValueMap(iEvent, src, drStaIn,        "drStaIn");
    writeValueMap(iEvent, src, dphiStaIn,      "dphiStaIn");
    writeValueMap(iEvent, src, drM2,           "drM2");
    writeValueMap(iEvent, src, dphiM2,         "dphiM2");
    writeValueMap(iEvent, src, distM2,         "distM2");
}

void
NearbyMuonsInfo::writeValueMap(edm::Event &iEvent,
        const edm::Handle<edm::View<reco::Candidate> > & handle,
        const std::vector<float> & values,
        const std::string    & label) const 
{
    using namespace edm; 
    using namespace std;
    auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*valMap);
    filler.insert(handle, values.begin(), values.end());
    filler.fill();
    iEvent.put(valMap, label);
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(NearbyMuonsInfo);
