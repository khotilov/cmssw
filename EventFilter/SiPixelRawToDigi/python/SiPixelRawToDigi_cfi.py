import FWCore.ParameterSet.Config as cms

siPixelDigis = cms.EDProducer("SiPixelRawToDigi",
    Timing = cms.untracked.bool(False),
    IncludeErrors = cms.bool(True),
    InputLabel = cms.InputTag("siPixelRawData"),
    CheckPixelOrder = cms.bool(False),
    UseCablingTree = cms.untracked.bool(True)
)



