import FWCore.ParameterSet.Config as cms

process = cms.Process("DTDQMOfflineSources")

# the source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/data/Commissioning08/Monitor/RAW/v1/000/067/818/E4AEAE98-B0A4-DD11-8414-0019B9F707D8.root',
      '/store/data/Commissioning08/Monitor/RAW/v1/000/067/818/E4EED1F7-DEA4-DD11-BC31-001D09F2AD4D.root',
      '/store/data/Commissioning08/Monitor/RAW/v1/000/067/818/E4FD1F4A-CBA4-DD11-9E99-001D09F23174.root',
      '/store/data/Commissioning08/Monitor/RAW/v1/000/067/818/E67F79EA-E3A4-DD11-94BD-001D09F2A465.root',
      '/store/data/Commissioning08/Monitor/RAW/v1/000/067/818/E6D73E61-A0A4-DD11-A769-000423D98844.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
    )


#process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("CondCore.DBCommon.CondDBSetup_cfi")


# Conditions (Global Tag is used here):

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cfi")
process.GlobalTag.globaltag = "CRAFT_ALL_V4::All"
#process.prefer("GlobalTag")

# Magnetic fiuld: force mag field to be 3.8 tesla
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

#Geometry
process.load("Configuration.StandardSequences.Geometry_cff")


# Real data raw to digi
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")

# reconstruction sequence for Cosmics
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
#process.load("RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff")
#process.load("RecoMuon.Configuration.RecoMuonCosmics_cff")



# offline DQM
process.load("DQM.DTMonitorModule.dtDQMOfflineSources_cff")
process.load("DQMServices.Components.MEtoEDMConverter_cff")



process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *', 
                                                                      'keep *_MEtoEDMConverter_*_*'),
                               fileName = cms.untracked.string('DTDQMOffline.root')
                               )








# message logger
process.MessageLogger = cms.Service("MessageLogger",
                                    debugModules = cms.untracked.vstring('*'),
                                    destinations = cms.untracked.vstring('cout'),
                                    categories = cms.untracked.vstring('DTChamberEfficiency'), 
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string('WARNING'),
                                                              noLineBreaks = cms.untracked.bool(False),
                                                              DEBUG = cms.untracked.PSet(
                                                                      limit = cms.untracked.int32(0)),
                                                              INFO = cms.untracked.PSet(
                                                                      limit = cms.untracked.int32(0)),
                                                              DTSegmentAnalysisTest = cms.untracked.PSet(
                                                                                 limit = cms.untracked.int32(-1)),
                                                              DTChamberEfficiency = cms.untracked.PSet(
                                                                                 limit = cms.untracked.int32(-1))
                                                              )
                                    )




# raw to digi
process.unpackers = cms.Sequence(process.muonDTDigis + process.muonCSCDigis + process.muonRPCDigis)
# reco
#process.reco = cms.Sequence(process.dt1DRecHits + process.dt4DSegments + process.muonRecoGR)
process.reco = cms.Sequence(process.offlineBeamSpot + process.muonsLocalRecoCosmics + process.STAmuontrackingforcosmics1Leg)

process.DTDQMOfflineCosmics = cms.Sequence(process.dtSources)


#Paths
process.allPath = cms.Path(process.unpackers *
                           process.reco *
                           process.DTDQMOfflineCosmics *
                           process.MEtoEDMConverter)


process.outpath = cms.EndPath(process.out)




# f = file('aNewconfigurationFile.cfg', 'w')
# f.write(process.dumpConfig())
# f.close()


