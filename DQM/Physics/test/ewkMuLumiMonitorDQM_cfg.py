import FWCore.ParameterSet.Config as cms

process = cms.Process("EwkDQM")
process.load("DQM.Physics.ewkMuLumiMonitorDQM_cfi")


process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")





process.DQM.collectorHost = ''

process.dqmSaver.workflow = cms.untracked.string('/Physics/EWK/LumiMonitor')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(        
#          'file:/ciet3b/data4/MUSKIM2010/MUSKIM_143958-144114_6.root'
            "file:/ciet3b/data4/MUSKIM2010B/MUAODRED_SKIM_147117-148058_5.root"
    )
)


process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('detailedInfo'),
    detailedInfo = cms.untracked.PSet(
            default = cms.untracked.PSet( limit = cms.untracked.int32(1) ),
            threshold = cms.untracked.string('INFO')
            #threshold = cms.untracked.string('ERROR'),
     )
)


process.p = cms.Path(
                     process.ewkMuLumiMonitorDQM* 
                     process.dqmSaver
                     )

