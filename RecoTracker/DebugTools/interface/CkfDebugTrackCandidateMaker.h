#ifndef CkfDebugTrackCandidateMaker_h
#define CkfDebugTrackCandidateMaker_h

#include "RecoTracker/CkfPattern/interface/CkfTrackCandidateMakerBase.h"
#include "RecoTracker/DebugTools/interface/CkfDebugTrajectoryBuilder.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Common/interface/EDProduct.h"

namespace cms {
  class CkfDebugTrackCandidateMaker : public CkfTrackCandidateMakerBase, public edm::EDProducer {
  public:
    CkfDebugTrackCandidateMaker(const edm::ParameterSet& conf) : CkfTrackCandidateMakerBase(conf) {
      produces<TrackCandidateCollection>();
    }

    virtual void beginJob (edm::EventSetup const & es){
      beginJobBase(es); 
      initDebugger(es);
    }

    virtual void produce(edm::Event& e, const edm::EventSetup& es){produceBase(e,es);}
    virtual void endJob() {delete dbg; }

  private:
    virtual TrajectorySeedCollection::const_iterator 
      lastSeed(TrajectorySeedCollection& theSeedColl){return theSeedColl.begin()+1;}

    void initDebugger(edm::EventSetup const & es){
      dbg = new CkfDebugger(es);
/*       myTrajectoryBuilder = dynamic_cast<const CkfDebugTrajectoryBuilder*>(theTrajectoryBuilder); */
/*       if (myTrajectoryBuilder) myTrajectoryBuilder->setDebugger( dbg); */
/* 	else */
	  theTrajectoryBuilder->setDebugger( dbg);
    };
    
    void printHitsDebugger(edm::Event& e){dbg->printSimHits(e);};
    void countSeedsDebugger(){dbg->countSeed();};
    void deleteAssocDebugger(){dbg->deleteHitAssociator();};
    void deleteDebugger(){delete dbg;};
    CkfDebugger *  dbg;
/*     const CkfDebugTrajectoryBuilder* myTrajectoryBuilder; */
  };
}

#endif
