import FWCore.ParameterSet.Config as cms

# Module to generate response of L1Trigger/RPCTrigger
l1RpcEmulDigis = cms.EDFilter("RPCTrigger",
    RPCTriggerDebug = cms.untracked.int32(0),
    label = cms.string('muonRPCDigis')
)


