import FWCore.ParameterSet.Config as cms

dataIntegrityTest = cms.EDAnalyzer("DTDataIntegrityTest",
    outputFile = cms.untracked.string('DataIntegrityTest.root'),
    folderRoot = cms.untracked.string(''),
    nTimeBin = cms.untracked.int32(9),
    doTimeHisto = cms.untracked.bool(False),
    writeHisto = cms.untracked.bool(False),
    debug = cms.untracked.bool(False),
    diagnosticPrescale = cms.untracked.int32(1)
)


