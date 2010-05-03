import FWCore.ParameterSet.Config as cms

trajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('TrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    fractionShared = cms.double('0.19'),
    allowSharedFirstHit = cms.bool(True)
)
