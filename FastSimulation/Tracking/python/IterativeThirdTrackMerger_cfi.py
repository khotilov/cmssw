import FWCore.ParameterSet.Config as cms

iterativeThirdTrackMerging = cms.EDFilter("FastTrackMerger",
    TrackProducers = cms.VInputTag(cms.InputTag("iterativeThirdTrackCandidatesWithPairs"),
                                   cms.InputTag("iterativeThirdTracksWithPairs")),
    trackAlgo = cms.untracked.uint32(7),
    MinNumberOfTrajHits = cms.untracked.uint32(4),
    MaxLostTrajHits = cms.untracked.uint32(0)                                          
)


