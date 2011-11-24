import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalCreateTimeCalibrations")

# Global Tag -- for original timing calibrations
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.GlobalTag.globaltag = 'GR_R_44_V5::All'


process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("ecalCreateTimeCalibs.root"),
    closeFileFast = cms.untracked.bool(True)
    )


process.createTimeCalibs = cms.EDAnalyzer("EcalCreateTimeCalibrations",
  OutputFileName = cms.string("file:converted1.root"),##Not Used
  OutputTimeCalibFileName = cms.string("EcalTimeCalibConstants.xml"),##Name of output xml file with new constants
  FileNameStart = cms.string("ecalCreateTimeCalibs"),
  ZeroGlobalOffset = cms.bool(True),
  SubtractDBcalibs = cms.bool(True),
  BxIncludeExclude = cms.string("-1"),
  OrbitIncludeExclude = cms.string("-1"),
  TriggerBitIncludeExclude = cms.string("-1"),
  TechTriggerBitIncludeExclude = cms.string("-1"),
  LumiIncludeExclude = cms.string("-1"),
  RunIncludeExclude = cms.string("-1"),
  AvgTimeMin = cms.double(-5),
  AvgTimeMax = cms.double(5),
  MinHitAmpEB = cms.double(26),  
  MinHitAmpEE = cms.double(47), 
  MaxSwissCross = cms.double(0.95),
  MaxHitTimeEB = cms.double(5),
  MinHitTimeEB = cms.double(-5),
  MaxHitTimeEE = cms.double(7),
  MinHitTimeEE = cms.double(-7),
  EventsUsedFractionNum = cms.double(1),
  EventsUsedFractionDen = cms.double(1),
  InputFileNames = cms.vstring('file:/data/jared/data/EcalTiming/DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_RAW-RECO/rh-DoubleElectron_Run2011A-ZElectron-PromptSkim-v4_RAW-RECO-jul15.HADDED.root')
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource",
       numberEventsInRun = cms.untracked.uint32(1),
       firstRun = cms.untracked.uint32(888888)
)
process.p = cms.Path(process.createTimeCalibs)
