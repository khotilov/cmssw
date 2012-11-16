import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("SiStrpDQMQTestTuning")

#prepare options

options = VarParsing.VarParsing("analysis")

options.register ('globalTag',
                                    "DONOTEXIST::All",
                                    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                                    VarParsing.VarParsing.varType.string,          # string, int, or float
                                    "GlobalTag")
options.register ('dqmFile',
                                    "",
                                    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                                    VarParsing.VarParsing.varType.string,          # string, int, or float
                                    "DQM root file")
options.register ('runNumber',
                                    0,
                                    VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                                    VarParsing.VarParsing.varType.int,          # string, int, or float
                                    "run number")

options.parseArguments()

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('siStripDigis', 
                                         'siStripClusters', 
                                         'siStripZeroSuppression', 
                                         'SiStripClusterizer',
                                         'siStripOfflineAnalyser'),
    cout = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
    destinations = cms.untracked.vstring('cout')
)


## Empty Event Source
process.source = cms.Source("EmptyIOVSource",
                              timetype = cms.string('runnumber'),
                              firstValue= cms.uint64(options.runNumber),
                              lastValue= cms.uint64(options.runNumber),
                              interval = cms.uint64(1)
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = options.globalTag

# DQM Environment
process.load("DQMServices.Core.DQMStore_cfg")

# SiStrip Offline DQM Client
process.siStripOfflineAnalyser = cms.EDAnalyzer("SiStripOfflineDQM",
       GlobalStatusFilling      = cms.untracked.int32(-1),
#        GlobalStatusFilling      = cms.untracked.int32(2),
        SummaryCreationFrequency  = cms.untracked.int32(-1),                                              
#       CreateSummary            = cms.untracked.bool(False),
       SummaryConfigPath        = cms.untracked.string("DQM/SiStripMonitorClient/data/sistrip_monitorelement_config.xml"),
       UsedWithEDMtoMEConverter = cms.untracked.bool(False),
       PrintFaultyModuleList    = cms.untracked.bool(False),

      InputFileName            = cms.untracked.string(options.dqmFile),
       OutputFileName           = cms.untracked.string("/tmp/testRunNum.root"), 
       CreateTkMap              = cms.untracked.bool(True),
       TkmapParameters          = cms.untracked.PSet(
          loadFedCabling    = cms.untracked.bool(True),
          trackerdatPath    = cms.untracked.string('CommonTools/TrackerMap/data/'),
          trackermaptxtPath = cms.untracked.string('CommonTools/TrackerMap/data/'),
          mapMin            = cms.untracked.double(0.)
       ),
       TkMapOptions             = cms.untracked.VPSet(
    cms.PSet(mapName=cms.untracked.string('QTestAlarm')),
    cms.PSet(mapName=cms.untracked.string('FractionOfBadChannels'),mapMax=cms.untracked.double(-1.),logScale=cms.untracked.bool(True)),
    cms.PSet(mapName=cms.untracked.string('NumberOfCluster')),
    cms.PSet(mapName=cms.untracked.string('NumberOfDigi')),
    cms.PSet(mapName=cms.untracked.string('NumberOfOfffTrackCluster')),
    cms.PSet(mapName=cms.untracked.string('NumberOfOnTrackCluster')),
    cms.PSet(mapName=cms.untracked.string('StoNCorrOnTrack')),
    cms.PSet(mapName=cms.untracked.string('NApvShots'),mapMax=cms.untracked.double(-1.),logScale=cms.untracked.bool(True)),
#    cms.PSet(mapName=cms.untracked.string('NApvShots'),mapMax=cms.untracked.double(-1.),psuMap=cms.untracked.bool(True),loadLVCabling=cms.untracked.bool(True)),
    cms.PSet(mapName=cms.untracked.string('MedianChargeApvShots'),mapMax=cms.untracked.double(-1.))
    )
)

# Services needed for TkHistoMap
process.TkDetMap = cms.Service("TkDetMap")
process.SiStripDetInfoFileReader = cms.Service("SiStripDetInfoFileReader")


process.p1 = cms.Path(process.siStripOfflineAnalyser)
