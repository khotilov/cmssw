import FWCore.ParameterSet.Config as cms

from TrackingTools.KalmanUpdators.KFUpdatorESProducer_cfi import *
from TrackingTools.GeomPropagators.SmartPropagator_cff import *
from RecoMuon.TrackingTools.MuonUpdatorAtVertex_cff import *
Chi2EstimatorForMuonTrackLoader = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2EstimatorForMuonTrackLoader'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(100000.0)
)

KFSmootherForMuonTrackLoader = cms.ESProducer("KFTrajectorySmootherESProducer",
    errorRescaling = cms.double(10.0),
    minHits = cms.int32(3),
    ComponentName = cms.string('KFSmootherForMuonTrackLoader'),
    Estimator = cms.string('Chi2EstimatorForMuonTrackLoader'),
    Updator = cms.string('KFUpdator'),
    Propagator = cms.string('SmartPropagatorAnyRK')
)

MuonTrackLoaderForSTA = cms.PSet(
    TrackLoaderParameters = cms.PSet(
        MuonUpdatorAtVertex,
        Smoother = cms.string('KFSmootherForMuonTrackLoader'),
        DoSmoothing = cms.bool(False),
        VertexConstraint = cms.bool(True)
    )
)
MuonTrackLoaderForGLB = cms.PSet(
    TrackLoaderParameters = cms.PSet(
        MuonUpdatorAtVertex,
        Smoother = cms.string('KFSmootherForMuonTrackLoader'),
        DoSmoothing = cms.bool(True),
        VertexConstraint = cms.bool(False)
    )
)
MuonTrackLoaderForL3 = cms.PSet(
    TrackLoaderParameters = cms.PSet(
        MuonUpdatorAtVertex,
        PutTkTrackIntoEvent = cms.untracked.bool(True),
        Smoother = cms.string('KFSmootherForMuonTrackLoader'),
        SmoothTkTrack = cms.untracked.bool(False),
        MuonSeededTracksInstance = cms.untracked.string('L2Seeded'),
        VertexConstraint = cms.bool(False),
        DoSmoothing = cms.bool(True)
    )
)
MuonTrackLoaderForCosmic = cms.PSet(
    TrackLoaderParameters = cms.PSet(
        MuonUpdatorAtVertexAnyDirection,
        PutTrajectoryIntoEvent = cms.untracked.bool(False),
        VertexConstraint = cms.bool(False),
        AllowNoVertex = cms.untracked.bool(True),
        Smoother = cms.string('KFSmootherForMuonTrackLoader'),
        DoSmoothing = cms.bool(False)
    )
)


