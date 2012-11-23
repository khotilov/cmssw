import FWCore.ParameterSet.Config as cms

simCastorDigis = cms.EDProducer("CastorDigiProducer",
    doNoise = cms.bool(True),
    doTimeSlew = cms.bool(True),
    castor = cms.PSet(
        readoutFrameSize = cms.int32(6),
        binOfMaximum = cms.int32(5),
        samplingFactor = cms.double(16.75), ## pe/GeV

        doPhotoStatistics = cms.bool(True),
        photoelectronsToAnalog = cms.double(4.24),
        simHitToPhotoelectrons = cms.double(1000.0),
        syncPhase = cms.bool(True),
        timePhase = cms.double(-4.0)
    )
)


