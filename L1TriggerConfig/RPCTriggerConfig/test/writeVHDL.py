import FWCore.ParameterSet.Config as cms

process = cms.Process("writeVHD")




process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

process.load("L1TriggerConfig.RPCTriggerConfig.L1RPCConfig_cff")
process.rpcconf.filedir = cms.untracked.string('L1TriggerConfig/RPCTriggerConfig/test/pats/') 

process.load("L1Trigger.RPCTrigger.RPCConeConfig_cff")
process.load("L1TriggerConfig.RPCTriggerConfig.RPCHwConfig_cff")
process.load("L1TriggerConfig.RPCTriggerConfig.RPCConeDefinition_cff")


process.load("EventFilter.RPCRawToDigi.RPCSQLiteCabling_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)


process.write = cms.EDAnalyzer("WriteVHDL",
          minTower = cms.int32(-12),
          maxTower = cms.int32(12),
          minSector = cms.int32(9),
          maxSector = cms.int32(9),
          templateName = cms.string("pacTemplate.vhd"),
          outDir =  cms.string("./vhd/")
)


process.p1 = cms.Path(process.write)
