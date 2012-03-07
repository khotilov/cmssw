import FWCore.ParameterSet.Config as cms

siPixelDigis = cms.EDProducer("SiPixelRawToDigi",
    Timing = cms.untracked.bool(False),
    IncludeErrors = cms.bool(True),
    InputLabel = cms.InputTag("siPixelRawData"),
    CheckPixelOrder = cms.bool(False),
    UseQualityInfo = cms.bool(False),
    UseCablingTree = cms.untracked.bool(True),
## ErrorList: list of error codes used by tracking to invalidate modules
    ErrorList = cms.vint32(29),
## UserErrorList: list of error codes used by Pixel experts for investigation
    UserErrorList = cms.vint32(40)
)



