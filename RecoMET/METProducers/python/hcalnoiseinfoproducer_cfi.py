import FWCore.ParameterSet.Config as cms

hcalnoise = cms.EDProducer(
    'HcalNoiseInfoProducer',

    # what to fill
    fillDigis = cms.bool(True),
    fillRecHits = cms.bool(True),
    fillCaloTowers = cms.bool(True),
    fillTracks = cms.bool(True),

    # maximum number of RBXs to fill
    # if you want to record all RBXs above some energy threshold,
    # change maxProblemRBXs to 999 and pMinE to the threshold you want
    maxProblemRBXs  = cms.int32(999),

    # parameters for calculating summary variables
    maxCaloTowerIEta = cms.int32(20),
    maxTrackEta = cms.double(2.0),
    minTrackPt = cms.double(1.0),

    # collection names
    digiCollName = cms.string('hcalDigis'),
    recHitCollName = cms.string('hbhereco'),
    caloTowerCollName = cms.string('towerMaker'),
    trackCollName = cms.string('generalTracks'),

    # define hit energy thesholds
    minRecHitE = cms.double(1.5),
    minLowHitE = cms.double(10.0),
    minHighHitE = cms.double(25.0),

    # define energy threshold for "problematic" cuts
    pMinERatio = cms.double(25.0),
    pMinEZeros = cms.double(5.0),
    pMinEEMF = cms.double(10.0),

    # define energy threshold for loose/tight/high level cuts
    minERatio = cms.double(50.0),
    minEZeros = cms.double(10.0),
    minEEMF = cms.double(20.0),

    # define problematic RBX
    pMinE = cms.double(5.0),
    pMinRatio = cms.double(0.75),
    pMaxRatio = cms.double(0.85),
    pMinHPDHits = cms.int32(10),
    pMinRBXHits = cms.int32(25),
    pMinHPDNoOtherHits = cms.int32(7),
    pMinZeros = cms.int32(5),
    pMinLowEHitTime = cms.double(-6.0),
    pMaxLowEHitTime = cms.double(6.0),
    pMinHighEHitTime = cms.double(-4.0),
    pMaxHighEHitTime = cms.double(5.0),
    pMaxHPDEMF = cms.double(0.02),
    pMaxRBXEMF = cms.double(0.02),

    # define loose noise cuts
    lMinRatio = cms.double(0.70),
    lMaxRatio = cms.double(0.90),
    lMinHPDHits = cms.int32(17),
    lMinRBXHits = cms.int32(999),
    lMinHPDNoOtherHits = cms.int32(10),
    lMinZeros = cms.int32(9),
    lMinLowEHitTime = cms.double(-9999.0),
    lMaxLowEHitTime = cms.double(9999.0),
    lMinHighEHitTime = cms.double(-8.0),
    lMaxHighEHitTime = cms.double(7.0),

    # define tight noise cuts
    tMinRatio = cms.double(0.73),
    tMaxRatio = cms.double(0.88),
    tMinHPDHits = cms.int32(13),
    tMinRBXHits = cms.int32(40),
    tMinHPDNoOtherHits = cms.int32(8),
    tMinZeros = cms.int32(6),
    tMinLowEHitTime = cms.double(-9999.0),
    tMaxLowEHitTime = cms.double(9999.0),
    tMinHighEHitTime = cms.double(-6.0),
    tMaxHighEHitTime = cms.double(5.0),

    # define high level noise cuts
    hlMaxHPDEMF = cms.double(-999.0),
    hlMaxRBXEMF = cms.double(0.01),

    )
