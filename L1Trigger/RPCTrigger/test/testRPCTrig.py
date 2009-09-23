import FWCore.ParameterSet.Config as cms

process = cms.Process("rpctest")


#process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
    log = cms.untracked.PSet( threshold = cms.untracked.string("DEBUG") ),
    debugModules = cms.untracked.vstring("l1RpcEmulDigis"),
    destinations = cms.untracked.vstring('log')
)

# rpc geometry
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

process.load("L1TriggerConfig.RPCTriggerConfig.RPCConeDefinition_cff")
# emulation
process.load("L1TriggerConfig.RPCTriggerConfig.L1RPCConfig_cff")

process.load("L1Trigger.RPCTrigger.RPCConeConfig_cff")
process.load("L1TriggerConfig.RPCTriggerConfig.RPCHwConfig_cff")
process.load("L1Trigger.RPCTrigger.l1RpcEmulDigis_cfi")
process.l1RpcEmulDigis.label = cms.string('simMuonRPCDigis')
process.l1RpcEmulDigis.RPCTriggerDebug = 1

# rpc r2d
#process.load("EventFilter.RPCRawToDigi.RPCSQLiteCabling_cfi")
#process.load("EventFilter.RPCRawToDigi.rpcUnpacker_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

#process.source = cms.Source("NewEventStreamFileReader",
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

            'file:digi.root'

    )
)

#process.a = cms.Path(process.rpcunpacker*process.l1RpcEmulDigis)
process.a = cms.Path(process.l1RpcEmulDigis)
