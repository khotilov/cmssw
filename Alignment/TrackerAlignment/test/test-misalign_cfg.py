# The following comments couldn't be translated into the new config version:

#            { string record = "TrackerAlignmentRcd" string tag = "TrackerIdealGeometry210" },
#            { string record = "TrackerAlignmentErrorRcd" string tag = "TrackerIdealGeometryErrors210" }

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
# -- Load default module/services configurations -- //
# Message logger service
process.load("FWCore.MessageService.MessageLogger_cfi")

#replace MessageLogger.debugModules = { "*" }
# service = Tracer {}
# Ideal geometry producer
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

# Database output service
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# Misalignment example scenario producer
#  process.load("Alignment.TrackerAlignment.ExampleScenario_cff")
#  process.load("Alignment.TrackerAlignment.NoMovementsScenario_cff")
#  process.load("Alignment.TrackerAlignment.Tracker1000pbScenario_cff")
#  process.load("Alignment.TrackerAlignment.Tracker100pbScenario_cff")
process.load("Alignment.TrackerAlignment.Tracker10pbScenario_cff")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBSetup,
    # Writing to oracle needs the following shell variable setting (in zsh):
    # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
    # connect = cms.string('oracle://cms_orcoff_prep/CMS_COND_ALIGNMENT'),  # preparation/develop. DB
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:Alignments.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('TrackerAlignmentRcd'),
        tag = cms.string('Tracker10pbScenario210_mc')
    ), 
        cms.PSet(
            record = cms.string('TrackerAlignmentErrorRcd'),
            tag = cms.string('Tracker10pbScenarioErrors210_mc')
        ))
)

process.prod = cms.EDAnalyzer("TestAnalyzer",
    fileName = cms.untracked.string('misaligned.root')
)

process.p1 = cms.Path(process.prod)
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )
)
process.CondDBSetup.DBParameters.messageLevel = 2
process.MisalignedTracker.saveToDbase = True


