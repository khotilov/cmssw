import FWCore.ParameterSet.Config as cms

# Global Pixel seeding
from FastSimulation.Tracking.GlobalPixelSeedProducer_cff import *
from FastSimulation.Tracking.GlobalPixelSeedProducerForElectrons_cff import *

# TrackCandidates
import FastSimulation.Tracking.TrackCandidateProducer_cfi

globalPixelTrackCandidates = FastSimulation.Tracking.TrackCandidateProducer_cfi.trackCandidateProducer.clone()
globalPixelTrackCandidates.SeedProducer = cms.InputTag("globalPixelSeeds","GlobalPixel")

globalPixelTrackCandidatesForElectrons = FastSimulation.Tracking.TrackCandidateProducer_cfi.trackCandidateProducer.clone()
globalPixelTrackCandidatesForElectrons.SeedProducer = cms.InputTag("globalPixelSeedsForElectrons","GlobalPixel")

# reco::Tracks (possibly with invalid hits)
import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cfi

globalPixelWithMaterialTracks = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cfi.ctfWithMaterialTracks.clone()
globalPixelWithMaterialTracks.src = 'globalPixelTrackCandidates'
globalPixelWithMaterialTracks.TTRHBuilder = 'WithoutRefit'
globalPixelWithMaterialTracks.Fitter = 'KFFittingSmoother'
globalPixelWithMaterialTracks.Propagator = 'PropagatorWithMaterial'
globalPixelWithMaterialTracks.TrajectoryInEvent = cms.bool(True)

globalPixelWithMaterialTracksForElectrons = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cfi.ctfWithMaterialTracks.clone()
globalPixelWithMaterialTracksForElectrons.src = 'globalPixelTrackCandidatesForElectrons'
globalPixelWithMaterialTracksForElectrons.TTRHBuilder = 'WithoutRefit'
globalPixelWithMaterialTracksForElectrons.Fitter = 'KFFittingSmoother'
globalPixelWithMaterialTracksForElectrons.Propagator = 'PropagatorWithMaterial'
globalPixelWithMaterialTracksForElectrons.TrajectoryInEvent = cms.bool(True)

# The sequences
globalPixelTracking = cms.Sequence(globalPixelSeeds*
                                   globalPixelTrackCandidates*
                                   globalPixelWithMaterialTracks)


globalPixelTrackingForElectrons = cms.Sequence(globalPixelSeedsForElectrons*
                                               globalPixelTrackCandidatesForElectrons*
                                               globalPixelWithMaterialTracksForElectrons)
