import FWCore.ParameterSet.Config as cms

siPixelRawData = cms.EDFilter("SiPixelDigiToRaw",
    InputLabel = cms.untracked.InputTag("simSiPixelDigis")
)



