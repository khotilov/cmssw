import FWCore.ParameterSet.Config as cms

ecalTimeTree = cms.EDAnalyzer("EcalTimePi0TreeMaker",
    barrelEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE"),

    useRaw = cms.untracked.bool(False),

    barrelEcalUncalibratedRecHitCollection = cms.InputTag("ecalRatioUncalibRecHit","EcalUncalibRecHitsEB"),
    endcapEcalUncalibratedRecHitCollection = cms.InputTag("ecalRatioUncalibRecHit","EcalUncalibRecHitsEE"),

    # gf set correct cluster producrs
    # gf here SC are used, switch to BC for us
    barrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters",""),
    endcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower",""),
    barrelBasicClusterCollection = cms.InputTag("correctedHybridSuperClusters",""),
    endcapBasicClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower",""),
    barrelClusterShapeAssociationCollection = cms.InputTag("hybridSuperClusters","hybridShapeAssoc"),
    endcapClusterShapeAssociationCollection = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapShapeAssoc"),

    vertexCollection  = cms.InputTag("offlinePrimaryVertices",""),
                               
    muonCollection = cms.InputTag("GLBMuons"),
    hbTreshold = cms.double(1.),                               
    l1GlobalReadoutRecord = cms.string('gtDigis'),
    GTRecordCollection = cms.untracked.string('gtDigis'),
    runNum = cms.untracked.int32(-1),
    fileName = cms.untracked.string('EcalTimePi0Tree'),
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dREcal = cms.double(9999.0),
        dRPreshowerPreselection = cms.double(9999.0),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        usePreshower = cms.bool(True),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        trajectoryUncertaintyTolerance = cms.double(-1),
        accountForTrajectoryChangeCalo = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    )
)


