import FWCore.ParameterSet.Config as cms

from RecoMuon.GlobalTrackingTools.MuonTrackingRegionCommon_cff import *
from RecoMuon.TransientTrackingRecHit.MuonTransientTrackingRecHitBuilder_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
GlobalTrajectoryBuilderCommon = cms.PSet(
    MuonTrackingRegionCommon,
    Chi2ProbabilityCut = cms.double(30.0),
    Direction = cms.int32(0),
    Chi2CutCSC = cms.double(150.0),
    HitThreshold = cms.int32(1),
    MuonHitsOption = cms.int32(1),
    ScaleTECxFactor = cms.double(-1.0),
    ScaleTECyFactor = cms.double(-1.0),
    TrackRecHitBuilder = cms.string('WithTrackAngle'),
    DTRecSegmentLabel = cms.InputTag("dt4DSegments"),
    Chi2CutRPC = cms.double(1.0),
    CSCRecSegmentLabel = cms.InputTag("cscSegments"),
    Chi2CutDT = cms.double(10.0),
    TrackerRecHitBuilder = cms.string('WithTrackAngle'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    RefitRPCHits = cms.bool(True),
    TrackTransformer = cms.PSet(
        Fitter = cms.string('KFFitterForRefitInsideOut'),
        TrackerRecHitBuilder = cms.string('WithTrackAngle'),
        Smoother = cms.string('KFSmootherForRefitInsideOut'),
        MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
        RefitDirection = cms.string('insideOut'),
        RefitRPCHits = cms.bool(True),
        DoPredictionsOnly = cms.bool(False)
    ),
    GlobalMuonTrackMatcher = cms.PSet(
        MinP = cms.double(2.5),
        Chi2Cut = cms.double(50.0),
        MinPt = cms.double(1.0),
        DeltaDCut = cms.double(10.0),
        DeltaRCut = cms.double(0.2)
    ),
    PtCut = cms.double(1.0),
    TrackerPropagator = cms.string('SteppingHelixPropagatorAny'),
    RPCRecSegmentLabel = cms.InputTag("rpcRecHits")
)

