import FWCore.ParameterSet.Config as cms

L3TrackCombiner = cms.EDProducer("L3TrackCombiner",
    labels = cms.VInputTag()
)
