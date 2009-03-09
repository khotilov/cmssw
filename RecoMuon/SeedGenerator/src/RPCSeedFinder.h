#ifndef RecoMuon_SeedGenerator_RPCSeedFinder_H
#define RecoMuon_SeedGenerator_RPCSeedFinder_H


/** \class RPCSeedFinder
 *  
 *   \author Haiyun.Teng - Peking University
 *
 *  
 */


#include "RecoMuon/SeedGenerator/src/RPCSeedPattern.h"
#include <DataFormats/TrajectorySeed/interface/TrajectorySeed.h>
#include <RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <vector>
#include <algorithm>

namespace edm {class EventSetup;}

class RPCSeedFinder 
{
    typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
    typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
    typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;

    public:
        RPCSeedFinder();
        ~RPCSeedFinder();
        void configure(const edm::ParameterSet& iConfig);
        void setOutput(std::vector<TrajectorySeed> *goodRef, std::vector<TrajectorySeed> *candidateRef);
        void setrecHits(ConstMuonRecHitContainer &recHits);
        void setEventSetup(const edm::EventSetup& iSetup);
        void seed();

    private:
        // Signal for call fillLayers()
        bool isrecHitsset;
        bool isConfigured;
        bool isOutputset;
        bool isEventSetupset;
        const edm::EventSetup *eSetup;
        RPCSeedPattern oneSeed;
        //ConstMuonRecHitContainer theRecHits;
        std::vector<TrajectorySeed> *goodSeedsRef;
        std::vector<TrajectorySeed> *candidateSeedsRef;
};
#endif
