# The following comments couldn't be translated into the new config version:

# All/OuterSurface/InnerSurface/ImpactPoint/default(track)
#

import FWCore.ParameterSet.Config as cms

# MonitorTrackGlobal
MonitorTrackSTAMuons = cms.EDFilter("TrackingMonitor",
    OutputMEsInRootFile = cms.bool(False),
    beamSpot = cms.InputTag("offlineBeamSpot"),                
    doTrackerSpecific = cms.bool(False),
    doAllPlots = cms.bool(False),
    doSeedParameterHistos = cms.bool(False), 
    SeedProducer = cms.InputTag("combinedP5SeedsForCTF"),
    TCProducer = cms.InputTag("ckfTrackCandidatesP5"),
    phiErrMax = cms.double(1.0),
    MeasurementState = cms.string('default'),
    etaErrMax = cms.double(0.5),
    TrackPxBin = cms.int32(1000),
    VXMax = cms.double(20.0),
    VXBin = cms.int32(20),
    RecHitBin = cms.int32(120),
    VZBin = cms.int32(50),
    TrackPzMin = cms.double(-500.0),
    Chi2Max = cms.double(500.0),
    Chi2Bin = cms.int32(250),
    TrackPzBin = cms.int32(1000),
    pxErrBin = cms.int32(100),
    VYMin = cms.double(-20.0),
    etaErrMin = cms.double(0.0),
    TrackPxMax = cms.double(500.0),
    TrackPzMax = cms.double(500.0),
    VYMax = cms.double(20.0),
    ThetaBin = cms.int32(100),
    RecHitMin = cms.double(0.0),
    OutputFileName = cms.string('monitortrackparameters_stamuons.root'),
    #
    FolderName = cms.string('Muons/standAloneMuonsUpdatedAtVtx'),
    pzErrMin = cms.double(0.0),
    EtaMin = cms.double(-3.0),
    pErrBin = cms.int32(100),
    pyErrMin = cms.double(0.0),
    VZMax = cms.double(100.0),
    phiErrBin = cms.int32(100),
    EtaMax = cms.double(3.0),
    etaErrBin = cms.int32(100),
    Chi2Min = cms.double(0.0),
    ThetaMin = cms.double(0.0),
    PhiMin = cms.double(-3.2),
    TrackPtMax = cms.double(1000),
    RecHitMax = cms.double(120.0),
    TrackPyBin = cms.int32(1000),
    ptErrMin = cms.double(0.0),
    TkSizeMin = cms.double(-0.5),
    TrackPxMin = cms.double(-500.0),
    VXMin = cms.double(-20.0),
    pyErrMax = cms.double(10.0),
    ThetaMax = cms.double(3.2),
    pzErrMax = cms.double(10.0),
    pzErrBin = cms.int32(100),
    pErrMin = cms.double(0.0),
    EtaBin = cms.int32(100),
    VZMin = cms.double(-100.0),
    pErrMax = cms.double(10.0),
    VYBin = cms.int32(20),
    pxErrMax = cms.double(10.0),
    TkSizeBin = cms.int32(11),
    AlgoName = cms.string('sta'),
    TrackPyMin = cms.double(-500.0),
    D0Max = cms.double(0.5),
    D0Bin = cms.int32(100),
    TrackProducer = cms.InputTag("standAloneMuons","UpdatedAtVtx"),
    TkSizeMax = cms.double(10.5),
    TrackPtBin = cms.int32(1000),
    TrackPyMax = cms.double(500.0),
    D0Min = cms.double(-0.5),
    phiErrMin = cms.double(0.0),
    pyErrBin = cms.int32(100),
    pxErrMin = cms.double(0.0),
    PhiBin = cms.int32(36),
    ptErrMax = cms.double(10.0),
    PhiMax = cms.double(3.2),
    TrackPtMin = cms.double(0),
    ptErrBin = cms.int32(100),
    RecLayBin = cms.int32(120),
    RecLayMax = cms.double(120.0),
    RecLayMin = cms.double(0.0),
    RecLostBin = cms.int32(120),
    RecLostMin = cms.double(0.0),
    RecLostMax = cms.double(20),
    TTRHBuilder = cms.string('WithTrackAngle'),
    TkSeedSizeBin = cms.int32(20),
    TkSeedSizeMin = cms.double(-0.5),
    TkSeedSizeMax = cms.double(19.5),
    TrackPBin = cms.int32(1000),
    TrackPMin = cms.double(0),
    TrackPMax = cms.double(1000),
    Chi2ProbMax = cms.double(1.0),
    Chi2ProbBin = cms.int32(100),
    Chi2ProbMin = cms.double(0.0),

    BSFolderName = cms.string('Muons/standAloneMuonsUpdatedAtVtx/BeamSpotParameters'),
    doBeamSpotPlots = cms.bool(False),
    X0Bin = cms.int32(100),
    X0Min = cms.double(-3.0),
    X0Max = cms.double(3.0),

    Y0Bin = cms.int32(100),
    Y0Min = cms.double(-3.0),
    Y0Max = cms.double(3.0),

    Z0Bin = cms.int32(60),
    Z0Min = cms.double(-30.0),
    Z0Max = cms.double(30.0)
)


