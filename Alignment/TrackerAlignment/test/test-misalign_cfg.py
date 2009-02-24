import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# -- Load default module/services configurations -- //
# Message logger service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )
)
#replace MessageLogger.debugModules = { "*" }

# service = Tracer {}
# Ideal geometry producer
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

# Misalignment example scenario producer
process.load("Alignment.TrackerAlignment.MisalignedTracker_cfi")
process.MisalignedTracker.saveToDbase = True # to store to DB
import Alignment.TrackerAlignment.Scenarios_cff as Scenarios
process.MisalignedTracker.scenario = Scenarios.Tracker10pbScenario
#process.MisalignedTracker.scenario = Scenarios.SurveyLASOnlyScenario
#process.MisalignedTracker.scenario = Scenarios.SurveyLASCosmicsScenario
#process.MisalignedTracker.scenario = Scenarios.TrackerCRAFTScenario

# the module
process.prod = cms.EDAnalyzer("TestAnalyzer",
    fileName = cms.untracked.string('misaligned.root')
)

# data loop
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Database output service
import CondCore.DBCommon.CondDBSetup_cfi
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    CondCore.DBCommon.CondDBSetup_cfi.CondDBSetup,
    # Writing to oracle needs the following shell variable setting (in zsh):
    # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
    # connect = cms.string('oracle://cms_orcoff_prep/CMS_COND_ALIGNMENT'),  # preparation/develop. DB
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:Alignments.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('TrackerAlignmentRcd'),
#        tag = cms.string('Tracker10pbScenario210_mc')
        tag = cms.string('TrackerCRAFTScenario22X_mc')
    ), 
        cms.PSet(
            record = cms.string('TrackerAlignmentErrorRcd'),
#        tag = cms.string('Tracker10pbScenarioErrors210_mc')
            tag = cms.string('TrackerCRAFTScenarioErrors22X_mc')
        ))
)
process.PoolDBOutputService.DBParameters.messageLevel = 2

process.p1 = cms.Path(process.prod)



