import FWCore.ParameterSet.Config as cms

from DQM.DTMonitorModule.dtChamberEfficiencyTask_cfi import *
from DQM.DTMonitorModule.dtSegmentTask_cfi import *
dtSegmentAnalysisMonitor.detailedAnalysis = True
from DQM.DTMonitorModule.dtResolutionTask_cfi import *

dqmInfoDT = cms.EDFilter("DQMEventInfo",
                         subSystemFolder = cms.untracked.string('DT')
                         )


dtDataIntegrityUnpacker = cms.EDProducer("DTUnpackingModule",
    dqmOnly = cms.untracked.bool(True),
    dataType = cms.string('DDU'),
    useStandardFEDid = cms.untracked.bool(True),
    fedbyType = cms.untracked.bool(True),
    readOutParameters = cms.PSet(
        debug = cms.untracked.bool(False),
        rosParameters = cms.PSet(
            writeSC = cms.untracked.bool(True),
            readingDDU = cms.untracked.bool(True),
            performDataIntegrityMonitor = cms.untracked.bool(True),
            readDDUIDfromDDU = cms.untracked.bool(True),
            debug = cms.untracked.bool(False),
            localDAQ = cms.untracked.bool(False)
        ),
        localDAQ = cms.untracked.bool(False),
        performDataIntegrityMonitor = cms.untracked.bool(True)
    )
)

from DQM.DTMonitorModule.dtDataIntegrityTask_cfi import *
DTDataIntegrityTask.hltMode = True

dtSources = cms.Sequence(dtDataIntegrityUnpacker *
                         dtSegmentAnalysisMonitor*
                         dtResolutionAnalysisMonitor*
                         dqmInfoDT)
