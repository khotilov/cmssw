import FWCore.ParameterSet.Config as cms

gsfTrackProducer = cms.EDProducer("GsfTrackProducer",
    src = cms.InputTag("CkfElectronCandidates"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    producer = cms.string(''),
    Fitter = cms.string('GsfElectronFittingSmoother'),
    useHitsSplitting = cms.bool(False),
    TrajectoryInEvent = cms.bool(False),
    TTRHBuilder = cms.string('WithTrackAngle'),
    Propagator = cms.string('fwdElectronPropagator'),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    AlgorithmName = cms.string('gsf')
)


