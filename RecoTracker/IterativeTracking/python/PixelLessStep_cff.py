import FWCore.ParameterSet.Config as cms

##########################################################################
# Large impact parameter tracking using TIB/TID/TEC stereo layer seeding #
##########################################################################


pixelLessStepClusters = cms.EDProducer("TrackClusterRemover",
    clusterLessSolution = cms.bool(True),
    oldClusterRemovalInfo = cms.InputTag("mixedTripletStepClusters"),
    trajectories = cms.InputTag("mixedTripletStepTracks"),
    overrideTrkQuals = cms.InputTag('mixedTripletStep'),
    TrackQuality = cms.string('highPurity'),
    pixelClusters = cms.InputTag("siPixelClusters"),
    stripClusters = cms.InputTag("siStripClusters"),
    Common = cms.PSet(
        maxChi2 = cms.double(9.0)
    )
)

# SEEDING LAYERS
pixelLessStepSeedLayers = cms.ESProducer("SeedingLayersESProducer",
    ComponentName = cms.string('pixelLessStepSeedLayers'),
    layerList = cms.vstring('TIB1+TIB2',
        'TID1_pos+TID2_pos','TID2_pos+TID3_pos',
        'TEC1_pos+TEC2_pos','TEC2_pos+TEC3_pos','TEC3_pos+TEC4_pos','TEC3_pos+TEC5_pos','TEC4_pos+TEC5_pos',
        'TID1_neg+TID2_neg','TID2_neg+TID3_neg',
        'TEC1_neg+TEC2_neg','TEC2_neg+TEC3_neg','TEC3_neg+TEC4_neg','TEC3_neg+TEC5_neg','TEC4_neg+TEC5_neg'),
    TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('pixelLessStepClusters')
    ),
    TID = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'),
        minRing = cms.int32(1),
        maxRing = cms.int32(2)
    ),
    TEC = cms.PSet(
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        skipClusters = cms.InputTag('pixelLessStepClusters'),
        useRingSlector = cms.bool(True),
        TTRHBuilder = cms.string('WithTrackAngle'),
        minRing = cms.int32(1),
        maxRing = cms.int32(2)
    )
)

# SEEDS
import RecoTracker.TkSeedGenerator.GlobalMixedSeeds_cff
pixelLessStepSeeds = RecoTracker.TkSeedGenerator.GlobalMixedSeeds_cff.globalMixedSeeds.clone()
pixelLessStepSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'pixelLessStepSeedLayers'
pixelLessStepSeeds.RegionFactoryPSet.RegionPSet.ptMin = 0.7
pixelLessStepSeeds.RegionFactoryPSet.RegionPSet.originHalfLength = 10.0
pixelLessStepSeeds.RegionFactoryPSet.RegionPSet.originRadius = 2.0

pixelLessStepSeeds.SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
        FilterAtHelixStage = cms.bool(True),
        FilterPixelHits = cms.bool(False),
        FilterStripHits = cms.bool(True),
        ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter')
    )

# QUALITY CUTS DURING TRACK BUILDING
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
pixelLessStepTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone(
    ComponentName = 'pixelLessStepTrajectoryFilter',
    filterPset = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.filterPset.clone(
    maxLostHits = 0,
    minimumNumberOfHits = 4,
    minPt = 0.1
    )
    )

import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
pixelLessStepChi2Est = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone(
    ComponentName = cms.string('pixelLessStepChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(16.0)
)

# TRACK BUILDING
import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi
pixelLessStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducer_cfi.GroupedCkfTrajectoryBuilder.clone(
    ComponentName = 'pixelLessStepTrajectoryBuilder',
    MeasurementTrackerName = '',
    clustersToSkip = cms.InputTag('pixelLessStepClusters'),
    trajectoryFilterName = 'pixelLessStepTrajectoryFilter',
    minNrOfHitsForRebuild = 4,
    maxCand = 2,
    alwaysUseInvalidHits = False,
    estimator = cms.string('pixelLessStepChi2Est'),
    maxPtForLooperReconstruction = cms.double(0.9) 
    )

# MAKING OF TRACK CANDIDATES
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
pixelLessStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('pixelLessStepSeeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    #onlyPixelHitsForSeedCleaner = cms.bool(True),

    TrajectoryBuilder = 'pixelLessStepTrajectoryBuilder'
)

# TRACK FITTING
import RecoTracker.TrackProducer.TrackProducer_cfi
pixelLessStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'pixelLessStepTrackCandidates',
    AlgorithmName = cms.string('iter5')
    )

import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
pixelLessStepSelector = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
    src='pixelLessStepTracks',
    trackSelectors= cms.VPSet(
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
            name = 'pixelLessStepLoose',
            chi2n_par = 0.5,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 4,
            maxNumberLostLayers = 1,
            minNumber3DLayers = 3,
            d0_par1 = ( 1.3, 4.0 ),
            dz_par1 = ( 1.3, 4.0 ),
            d0_par2 = ( 1.3, 4.0 ),
            dz_par2 = ( 1.3, 4.0 )
            ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
            name = 'pixelLessStepTight',
            preFilterName = 'pixelLessStepLoose',
            chi2n_par = 0.35,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 4,
            maxNumberLostLayers = 0,
            minNumber3DLayers = 3,
            d0_par1 = ( 1.1, 4.0 ),
            dz_par1 = ( 1.1, 4.0 ),
            d0_par2 = ( 1.1, 4.0 ),
            dz_par2 = ( 1.1, 4.0 )
            ),
        RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
            name = 'pixelLessStep',
            preFilterName = 'pixelLessStepTight',
            chi2n_par = 0.2,
            res_par = ( 0.003, 0.001 ),
            minNumberLayers = 4,
            maxNumberLostLayers = 0,
            minNumber3DLayers = 3,
            d0_par1 = ( 0.9, 4.0 ),
            dz_par1 = ( 0.9, 4.0 ),
            d0_par2 = ( 0.9, 4.0 ),
            dz_par2 = ( 0.9, 4.0 )
            ),
        ) #end of vpset
    ) #end of clone


PixelLessStep = cms.Sequence(pixelLessStepClusters*
                             pixelLessStepSeeds*
                             pixelLessStepTrackCandidates*
                             pixelLessStepTracks*
                             pixelLessStepSelector)
