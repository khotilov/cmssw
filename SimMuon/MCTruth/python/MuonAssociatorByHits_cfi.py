import FWCore.ParameterSet.Config as cms

muonAssociatorByHits = cms.EDProducer("MuonAssociatorEDProducer",
    # for Muon Track association
    #
    #     input collections
    #
    # ... reco::Track collection
    tracksTag = cms.InputTag("standAloneMuons"),
    # tracksTag = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    # tracksTag = cms.InputTag("globalMuons"),
    # tracksTag = cms.InputTag("generalTracks"),
    #
    # ... TrackingParticle collection
    tpTag = cms.InputTag("mergedtruth","MergedTrackTruth"),
    #
    dumpInputCollections = cms.bool(False),
    #
    #....... general input parameters
    #
    AbsoluteNumberOfHits_track = cms.bool(False),
    MinHitCut_track = cms.uint32(1),
    AbsoluteNumberOfHits_muon = cms.bool(False),
    MinHitCut_muon = cms.uint32(1),
    #
    PurityCut_track = cms.double(0.5),
    PurityCut_muon = cms.double(0.5),
    #
    SimToReco_useTracker = cms.bool(False),
    EfficiencyCut_track = cms.double(0.5),
    #
    SimToReco_useMuon = cms.bool(True),
    EfficiencyCut_muon = cms.double(0.5),
    #
    #........(for inner tracker stub of Global Muons)...
    UsePixels = cms.bool(True),
    UseGrouped = cms.bool(True),
    UseSplitting = cms.bool(True),
    ThreeHitTracksAreSpecial = cms.bool(True),
    #
    # for DT Hit associator
    crossingframe = cms.bool(True),
    simtracksTag = cms.InputTag("g4SimHits"),
    simtracksXFTag = cms.InputTag("mix","g4SimHits"),
    #
    DTsimhitsTag = cms.InputTag("g4SimHits","MuonDTHits"),
    DTsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonDTHits"),
    DTdigiTag = cms.InputTag("simMuonDTDigis"),
    DTdigisimlinkTag = cms.InputTag("simMuonDTDigis"),
    DTrechitTag = cms.InputTag("dt1DRecHits"),
    #
    dumpDT = cms.bool(False),
    links_exist = cms.bool(True),
    associatorByWire = cms.bool(False),
    #
    # for CSC Hit associator
    CSCsimHitsXFTag = cms.InputTag("mix","g4SimHitsMuonCSCHits"),
    CSClinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigiSimLinks"),
    CSCwireLinksTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigiSimLinks"),
    #
    # for RPC Hit associator
    RPCsimhitsXFTag = cms.InputTag("mix","g4SimHitsMuonRPCHits"),
    RPCdigisimlinkTag = cms.InputTag("simMuonRPCDigis","RPCDigiSimLink"),
    #
    # for Tracker Hit associator
    #
    associatePixel = cms.bool(True),
    associateStrip = cms.bool(True),
    associateRecoTracks = cms.bool(True),
    #                                
    ROUList = cms.vstring('TrackerHitsTIBLowTof', 
        'TrackerHitsTIBHighTof', 
        'TrackerHitsTIDLowTof', 
        'TrackerHitsTIDHighTof', 
        'TrackerHitsTOBLowTof', 
        'TrackerHitsTOBHighTof', 
        'TrackerHitsTECLowTof', 
        'TrackerHitsTECHighTof', 
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof')
)



