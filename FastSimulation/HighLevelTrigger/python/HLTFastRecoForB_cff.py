import FWCore.ParameterSet.Config as cms


from FastSimulation.Tracking.GlobalPixelTracking_cff import *

#############################################
# Reconstruct tracks with pixel seeds
#############################################

# Take all pixel tracks for b tagging track reco (pTMin>1GeV, nHits>=8) 
hltBLifetimeRegionalCtfWithMaterialTracks = cms.EDFilter("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("globalPixelWithMaterialTracks")),
    ptMin = cms.untracked.double(1.0),
    minHits = cms.untracked.uint32(8)
)

hltBLifetimeRegionalCtfWithMaterialTracksRelaxed = cms.EDFilter("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("globalPixelWithMaterialTracks")),
    ptMin = cms.untracked.double(1.0),
    minHits = cms.untracked.uint32(8)
)

hltBLifetimeRegionalCkfTrackCandidates = cms.Sequence(globalPixelTracking)
hltBLifetimeRegionalCkfTrackCandidatesRelaxed = cms.Sequence(globalPixelTracking)


#############################################
# Reconstruct muons for MumuK
#############################################
import FWCore.ParameterSet.Config as cms

# Take all pixel-seeded tracks for b tagging track reco (pTMin>1GeV, nHits>=8) 
hltCtfWithMaterialTracksMumuk = cms.EDFilter("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("globalPixelWithMaterialTracks")),
    ptMin = cms.untracked.double(3.0),
    minHits = cms.untracked.uint32(5)
)

# produce ChargedCandidates from tracks
hltMumukAllConeTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hltCtfWithMaterialTracksMumuk"),
    particleType = cms.string('mu-')
)

hltCkfTrackCandidatesMumuk = cms.Sequence(cms.SequencePlaceholder("HLTL3muonrecoSequence"))


#############################################
# Reconstruct muons for JPsiToMumu
#############################################

# Take all pixel-seeded tracks for b tagging track reco (pTMin>1GeV, nHits>=8) 
hltCtfWithMaterialTracksMumu = cms.EDFilter("FastTrackMerger",
    SaveTracksOnly = cms.untracked.bool(True),
    TrackProducers = cms.VInputTag(cms.InputTag("hltL3Muons")),
    ptMin = cms.untracked.double(3.0),
    minHits = cms.untracked.uint32(5)
)

# produce ChargedCandidates from tracks
hltMuTracks = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("hltCtfWithMaterialTracksMumu"),
    particleType = cms.string('mu-')
)

hltCkfTrackCandidatesMumu = cms.Sequence(cms.SequencePlaceholder("HLTL3muonrecoNocandSequence"))


