import FWCore.ParameterSet.Config as cms

dedxUnbinned = cms.EDProducer("DeDxEstimatorProducer",
    tracks                     = cms.InputTag("generalTracks"),
    trajectoryTrackAssociation = cms.InputTag("generalTracks"),

    estimator      = cms.string('unbinnedFit'),

    UseStrip       = cms.bool(True),
    UsePixel       = cms.bool(True),
    MeVperADCStrip = cms.double(3.61e-06*250),
    MeVperADCPixel = cms.double(3.61e-06)
)



