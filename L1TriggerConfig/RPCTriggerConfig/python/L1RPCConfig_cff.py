import FWCore.ParameterSet.Config as cms

rpcconf = cms.ESProducer("RPCTriggerConfig",
    filedir = cms.untracked.string('L1Trigger/RPCTrigger/data/Eff90PPT12/'),
    PACsPerTower = cms.untracked.int32(12)
)

rpcconfsrc = cms.ESSource("EmptyESSource",
    recordName = cms.string('L1RPCConfigRcd'),
    iovIsRunNotTime = cms.bool(True),
    firstValid = cms.vuint32(1)
)


