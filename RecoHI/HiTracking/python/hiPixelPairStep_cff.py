import FWCore.ParameterSet.Config as cms


#################################
# Filter on quality tracks
hiThirdStepFilter = cms.EDProducer("QualityFilter",
                                  TrackQuality = cms.string('highPurity'),
                                  recTracks = cms.InputTag("hiMixedTripletSelectedTracks")
                                  )

# NEW CLUSTERS (remove previously used clusters)
hiPixelPairClusters = cms.EDProducer("TrackClusterRemover",
                                clusterLessSolution= cms.bool(True),
                                oldClusterRemovalInfo = cms.InputTag("hiMixedTripletClusters"),
                                trajectories = cms.InputTag("hiThirdStepFilter"),
                                TrackQuality = cms.string('highPurity'),
                                pixelClusters = cms.InputTag("siPixelClusters"),
                                stripClusters = cms.InputTag("siStripClusters"),
                                Common = cms.PSet(
    maxChi2 = cms.double(9.0),
    ),
                                Strip = cms.PSet(
    maxChi2 = cms.double(9.0),
    #Yen-Jie's mod to preserve merged clusters
    maxSize = cms.uint32(2)   
    )
                                )


# SEEDING LAYERS
import RecoTracker.TkSeedingLayers.PixelLayerPairs_cfi
hiPixelPairSeedLayers = RecoTracker.TkSeedingLayers.PixelLayerPairs_cfi.pixellayerpairs.clone(
            ComponentName = 'hiPixelPairSeedLayers',
                        )

# SEEDS
import RecoTracker.TkSeedGenerator.GlobalSeedsFromPairsWithVertices_cff
hiPixelPairSeeds = RecoTracker.TkSeedGenerator.GlobalSeedsFromPairsWithVertices_cff.globalSeedsFromPairsWithVertices.clone()
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.VertexCollection=cms.InputTag("hiSelectedVertex")
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.ptMin = 4.0
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.originRadius = 0.005
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.nSigmaZ = 4.0
# only used for pixel tracking? -Matt
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.sigmaZVertex = 4.0
# Using a fixed error to determine the tracking region seems like a bad idea to me -Matt
hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.useFixedError = cms.bool(False)
# From pp settings, Ed's config used 0.005
#hiPixelPairSeeds.RegionFactoryPSet.RegionPSet.fixedError = 0.03
hiPixelPairSeeds.OrderedHitsFactoryPSet.SeedingLayers = cms.string('hiPixelPairSeedLayers')
hiPixelPairSeeds.OrderedHitsFactoryPSet.maxElement = 5000000
hiPixelPairSeeds.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000000
hiPixelPairSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 50000000

# QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
hiPixelPairTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone(
    ComponentName = 'hiPixelPairTrajectoryFilter',
    filterPset = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.filterPset.clone(
    minimumNumberOfHits = 6,
    minPt = 1.0
    )
    )

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
hiPixelPairChi2Est = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone(
        ComponentName = cms.string('hiPixelPairChi2Est'),
            nSigma = cms.double(3.0),
            MaxChi2 = cms.double(9.0)
        )

# TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi
hiPixelPairTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi.GroupedCkfTrajectoryBuilder.clone(
        ComponentName = 'hiPixelPairTrajectoryBuilder',
            MeasurementTrackerName = '',
            trajectoryFilterName = 'hiPixelPairTrajectoryFilter',
            clustersToSkip = cms.InputTag('hiPixelPairClusters'),
            maxCand = 3,
            estimator = cms.string('hiPixelPairChi2Est')
            )

# MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
hiPixelPairTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
        src = cms.InputTag('hiPixelPairSeeds'),
            TrajectoryBuilder = 'hiPixelPairTrajectoryBuilder'
            )


# TRACK FITTING
import RecoTracker.TrackProducer.TrackProducer_cfi
hiPixelPairGlobalPrimTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'hiPixelPairTrackCandidates',
    AlgorithmName = cms.string('iter3')
    )


#################################
# HI track selection
from RecoHI.HiTracking.HISelectedTracks_cfi import *
hiPixelPairSelectedTracks = hiSelectedTracks.clone(
    src = "hiPixelPairGlobalPrimTracks",
    min_nhits = cms.uint32(14)
    )


# Final sequence

hiPixelPairStep = cms.Sequence(hiThirdStepFilter*
                               hiPixelPairClusters*
                               hiPixelPairSeeds*
                               hiPixelPairTrackCandidates*
                               hiPixelPairGlobalPrimTracks*
                               hiPixelPairSelectedTracks)
