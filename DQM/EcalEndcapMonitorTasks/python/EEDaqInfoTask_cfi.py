import FWCore.ParameterSet.Config as cms

ecalEndcapDaqInfoTask = cms.EDAnalyzer("EEDaqInfoTask",
    prefixME = cms.untracked.string('Ecal'),
    enableCleanup = cms.untracked.bool(False),
    mergeRuns = cms.untracked.bool(False),
)
