// -*- C++ -*-
//
// Package:    RPCSeedGenerator
// Class:      RPCSeedGenerator
// 
/**\class RPCSeedGenerator RecoMuon/SeedGenerator/src/RPCSeedGenerator.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Haiyun Teng
//         Created:  Wed Oct 29 17:24:36 CET 2008
// $Id$
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// special include files
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include <vector>
// Using other classes
#include "RecoMuon/SeedGenerator/src/RPCSeedPattern.h"
#include "RecoMuon/SeedGenerator/src/RPCSeedFinder.h"
#include "RecoMuon/SeedGenerator/src/RPCSeedrecHitFinder.h"
#include "RecoMuon/SeedGenerator/src/RPCSeedLayerFinder.h"
#include "RecoMuon/SeedGenerator/src/RPCSeedOverlapper.h"
// Geometry
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
// Math
#include "math.h"
// C++
#include <vector>

//
// constants, enums and typedefs
//
using namespace std;
using namespace edm;

typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitPointer ConstMuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;
typedef MuonTransientTrackingRecHit::ConstMuonRecHitContainer ConstMuonRecHitContainer;


#ifndef RPCLayerNumber
#define RPCLayerNumber 12
#endif

#ifndef BarrelLayerNumber
#define BarrelLayerNumber 6
#endif

#ifndef EachEndcapLayerNumber
#define EachEndcapLayerNumber 3
#endif

//
// class decleration
//

class RPCSeedFinder;

class RPCSeedGenerator : public edm::EDProducer {
    public:
        explicit RPCSeedGenerator(const edm::ParameterSet& iConfig);
        ~RPCSeedGenerator();

    private:
        virtual void beginJob(const edm::EventSetup& iSetup);
        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
        virtual void endJob();

        void CheckOverlap(const edm::EventSetup& iSetup);
        void CheckcandidateOverlap(const edm::EventSetup& iSetup);
        bool isShareHit(const edm::OwnVector<TrackingRecHit> RecHits, const TrackingRecHit& hit, edm::ESHandle<RPCGeometry> rpcGeometry);

        // ----------member data ---------------------------
        RPCSeedFinder Finder;
        RPCSeedrecHitFinder recHitFinder;
        RPCSeedLayerFinder LayerFinder;
        RPCSeedOverlapper Overlapper;
        std::vector<TrajectorySeed> candidateSeeds;
        std::vector<TrajectorySeed> goodSeeds;
        edm::InputTag theRPCRecHits;
};


//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
RPCSeedGenerator::RPCSeedGenerator(const edm::ParameterSet& iConfig)
{
    //register your products
    /* Examples
       produces<ExampleData2>();

    //if do put with a label
    produces<ExampleData2>("label");
    */

    // Now do what ever other initialization is needed
    // Configure other modules
    Finder.configure(iConfig);
    recHitFinder.configure(iConfig);
    LayerFinder.configure(iConfig);
    Overlapper.configure(iConfig);
    // Register the production
    produces<TrajectorySeedCollection>("goodSeeds");
    produces<TrajectorySeedCollection>("candidateSeeds");
    // Get event data Tag
    theRPCRecHits = iConfig.getParameter<edm::InputTag>("RPCRecHitsLabel");
    
    cout << endl << "[RPCSeedGenerator] --> Constructor called" << endl;
}


RPCSeedGenerator::~RPCSeedGenerator()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    cout << "[RPCSeedGenerator] --> Destructor called" << endl;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
RPCSeedGenerator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    /* This is an event example
    //Read 'ExampleData' from the Event
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);

    //Use the ExampleData to create an ExampleData2 which 
    // is put into the Event
    std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
    iEvent.put(pOut);
    */

    /* this is an EventSetup example
    //Read SetupData from the SetupRecord in the EventSetup
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    */

    // clear goodSeeds from last reconstruction  
    goodSeeds.clear();
    candidateSeeds.clear();

    // Create the pointer to the Seed container
    auto_ptr<TrajectorySeedCollection> goodCollection(new TrajectorySeedCollection());
    auto_ptr<TrajectorySeedCollection> candidateCollection(new TrajectorySeedCollection());

    // Muon Geometry - DT, CSC and RPC 
    edm::ESHandle<MuonDetLayerGeometry> muonLayers;
    iSetup.get<MuonRecoGeometryRecord>().get(muonLayers);

    // Get the RPC layers
    vector<DetLayer*> RPCBarrelLayers = muonLayers->barrelRPCLayers();
    const DetLayer* RB4L  = RPCBarrelLayers[5];
    const DetLayer* RB3L  = RPCBarrelLayers[4];
    const DetLayer* RB22L = RPCBarrelLayers[3];
    const DetLayer* RB21L = RPCBarrelLayers[2];
    const DetLayer* RB12L = RPCBarrelLayers[1];
    const DetLayer* RB11L = RPCBarrelLayers[0];
    vector<DetLayer*> RPCEndcapLayers = muonLayers->endcapRPCLayers();
    const DetLayer* REM3L = RPCEndcapLayers[0];
    const DetLayer* REM2L = RPCEndcapLayers[1];
    const DetLayer* REM1L = RPCEndcapLayers[2];
    const DetLayer* REP1L = RPCEndcapLayers[3];
    const DetLayer* REP2L = RPCEndcapLayers[4];
    const DetLayer* REP3L = RPCEndcapLayers[5];

    // Get RPC recHits by MuonDetLayerMeasurements, while CSC and DT is set to false and with empty InputTag
    MuonDetLayerMeasurements muonMeasurements(edm::InputTag(), edm::InputTag(), theRPCRecHits, false, false, true);

    // Dispatch RPC recHits to the corresponding DetLayer, 6 layers for barrel and 3 layers for each endcap
    MuonRecHitContainer recHitsRPC[RPCLayerNumber];
    recHitsRPC[0] = muonMeasurements.recHits(RB11L, iEvent);
    recHitsRPC[1] = muonMeasurements.recHits(RB12L, iEvent);
    recHitsRPC[2] = muonMeasurements.recHits(RB21L, iEvent);
    recHitsRPC[3] = muonMeasurements.recHits(RB22L, iEvent);
    recHitsRPC[4] = muonMeasurements.recHits(RB3L, iEvent);
    recHitsRPC[5] = muonMeasurements.recHits(RB4L, iEvent); 
    recHitsRPC[6] = muonMeasurements.recHits(REM1L, iEvent);
    recHitsRPC[7] = muonMeasurements.recHits(REM2L, iEvent);
    recHitsRPC[8] = muonMeasurements.recHits(REM3L, iEvent);
    recHitsRPC[9] = muonMeasurements.recHits(REP1L, iEvent);
    recHitsRPC[10] = muonMeasurements.recHits(REP2L, iEvent);
    recHitsRPC[11] = muonMeasurements.recHits(REP3L, iEvent);

    // Print the size of recHits in each DetLayer
    cout << "RB1in "  << recHitsRPC[0].size()  << " recHits" << endl;
    cout << "RB1out " << recHitsRPC[1].size()  << " recHits" << endl;
    cout << "RB2in "  << recHitsRPC[2].size()  << " recHits" << endl;
    cout << "RB2out " << recHitsRPC[3].size()  << " recHits" << endl;
    cout << "RB3 "    << recHitsRPC[4].size()  << " recHits" << endl;
    cout << "RB4 "    << recHitsRPC[5].size()  << " recHits" << endl;
    cout << "REM1 "   << recHitsRPC[6].size()  << " recHits" << endl;
    cout << "REM2 "   << recHitsRPC[7].size()  << " recHits" << endl;
    cout << "REM3 "   << recHitsRPC[8].size()  << " recHits" << endl;
    cout << "REP1 "   << recHitsRPC[9].size()  << " recHits" << endl;
    cout << "REP2 "   << recHitsRPC[10].size() << " recHits" << endl;
    cout << "REP3 "   << recHitsRPC[11].size() << " recHits" << endl;

    // Set Input of RPCSeedFinder, PCSeedrecHitFinder, RPCSeedLayerFinder
    recHitFinder.setInput(recHitsRPC);
    LayerFinder.setInput(recHitsRPC);
    Overlapper.setIO(&goodSeeds, &candidateSeeds);
    // Start from filling layers to filling seeds
    LayerFinder.fillLayers();
    Overlapper.run();

    // Save seeds to event
    for(vector<TrajectorySeed>::iterator seed = goodSeeds.begin(); seed != goodSeeds.end(); ++seed)
        goodCollection->push_back(*seed);
    for(vector<TrajectorySeed>::iterator seed = candidateSeeds.begin(); seed != candidateSeeds.end(); ++seed)
        candidateCollection->push_back(*seed);

    // Put the seed to event
    iEvent.put(goodCollection, "goodSeeds");
    iEvent.put(candidateCollection, "candidateSeeds");

    // Unset the input of RPCSeedFinder, PCSeedrecHitFinder, RPCSeedLayerFinder
    recHitFinder.unsetInput();
    LayerFinder.unsetInput();
    Overlapper.unsetIO();
}

void RPCSeedGenerator::beginJob(const edm::EventSetup& iSetup) {

    // Set link and EventSetup of RPCSeedFinder, PCSeedrecHitFinder, RPCSeedLayerFinder
    cout << "set link and EventSetup of RPCSeedFinder, PCSeedrecHitFinder, RPCSeedLayerFinder" << endl;
    Finder.setEventSetup(iSetup);
    Finder.setOutput(&goodSeeds, &candidateSeeds);
    recHitFinder.setOutput(&Finder);
    LayerFinder.setOutput(&recHitFinder);
    Overlapper.setEventSetup(iSetup);
}

void RPCSeedGenerator::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(RPCSeedGenerator);
