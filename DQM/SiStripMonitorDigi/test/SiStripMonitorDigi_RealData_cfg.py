import FWCore.ParameterSet.Config as cms

process = cms.Process("MonitorDigiRealData")

#--------------------------
# Event Source
#--------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#      '/store/data/Commissioning08/Cosmics/RAW/v1/000/066/993/2A8E332E-769F-DD11-8E22-001617DBD224.root'
#      '/store/data/Commissioning08/Cosmics/RAW/v1/000/066/993/4A64178D-799F-DD11-9C47-000423D996B4.root'
#       '/store/data/Commissioning08/Cosmics/RAW/v1/000/066/668/48F6B6B5-519C-DD11-B32A-000423D6B444.root'
#        '/store/data/Commissioning08/Cosmics/RAW/v1/000/067/647/0000721C-35A3-DD11-9132-001D09F291D7.root'
       '/store/data/Commissioning08/Cosmics/RAW/v1/000/067/647/22CBBD11-07A3-DD11-9DFB-001D09F2447F.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

#-------------------------------------------------
# Message Logger
#-------------------------------------------------
process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siStripDigis'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    ),
    destinations = cms.untracked.vstring('cout')
)

#-------------------------------------------------
# Magnetic Field
#-------------------------------------------------
process.load("Configuration.GlobalRuns.ForceZeroTeslaField_cff")

#-------------------------------------------------
# Geometry
#-------------------------------------------------
process.load("Configuration.StandardSequences.Geometry_cff")

#-------------------------------------------------
# Calibration
#-------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.connect = "frontier://FrontierProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "CRAFT_V3P::All"
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

#-----------------------
#  Reconstruction Modules
#-----------------------
process.load("EventFilter.SiStripRawToDigi.SiStripDigis_cfi")
process.siStripDigis.ProductLabel = 'source'

#--------------------------
# DQM Services
#--------------------------
process.DQMStore = cms.Service("DQMStore",
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0)
)

#--------------------------
# SiStrip MonitorDigi
#--------------------------
process.load("DQM.SiStripMonitorDigi.SiStripMonitorDigi_cfi")
process.SiStripMonitorDigi.CreateTrendMEs = True
process.SiStripMonitorDigi.OutputMEsInRootFile = True
process.SiStripMonitorDigi.OutputFileName = 'SiStripMonitorDigi_RealData.root'
process.SiStripMonitorDigi.SelectAllDetectors = True

process.outP = cms.OutputModule("AsciiOutputModule")
process.AdaptorConfig = cms.Service("AdaptorConfig")

#--------------------------
# Sequences 
#--------------------------

process.RecoForDQM = cms.Sequence(process.siStripDigis)

process.p = cms.Path(process.RecoForDQM*process.SiStripMonitorDigi)
process.ep = cms.EndPath(process.outP)


