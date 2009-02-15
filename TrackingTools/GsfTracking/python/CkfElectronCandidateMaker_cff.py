import FWCore.ParameterSet.Config as cms

#Chi2 estimator
import TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi
ElectronChi2 = TrackingTools.KalmanUpdators.Chi2MeasurementEstimatorESProducer_cfi.Chi2MeasurementEstimator.clone()
ElectronChi2.ComponentName = 'ElectronChi2'
ElectronChi2.MaxChi2 = 2000.
ElectronChi2.nSigma = 3.

# Trajectory Filter
from TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff import *
import TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi
TrajectoryFilterForElectrons = TrackingTools.TrajectoryFiltering.TrajectoryFilterESProducer_cfi.trajectoryFilterESProducer.clone()
TrajectoryFilterForElectrons.ComponentName = 'TrajectoryFilterForElectrons'
TrajectoryFilterForElectrons.filterPset = cms.PSet(
    chargeSignificance = cms.double(-1.0),
    minPt = cms.double(2.0),
    minHitsMinPt = cms.int32(-1),
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    maxLostHits = cms.int32(1),
    maxNumberOfHits = cms.int32(-1),
    maxConsecLostHits = cms.int32(1),
    nSigmaMinPt = cms.double(5.0),
    minimumNumberOfHits = cms.int32(5)
)

# Trajectory Builder
import RecoTracker.CkfPattern.CkfTrajectoryBuilderESProducer_cfi
TrajectoryBuilderForElectrons = RecoTracker.CkfPattern.CkfTrajectoryBuilderESProducer_cfi.CkfTrajectoryBuilder.clone()
TrajectoryBuilderForElectrons.ComponentName = 'TrajectoryBuilderForElectrons'
TrajectoryBuilderForElectrons.trajectoryFilterName = 'TrajectoryFilterForElectrons'
TrajectoryBuilderForElectrons.maxCand = 5
TrajectoryBuilderForElectrons.intermediateCleaning = False
TrajectoryBuilderForElectrons.propagatorAlong = 'fwdGsfElectronPropagator'
TrajectoryBuilderForElectrons.propagatorOpposite = 'bwdGsfElectronPropagator'
TrajectoryBuilderForElectrons.estimator = 'ElectronChi2'
TrajectoryBuilderForElectrons.MeasurementTrackerName = ''
TrajectoryBuilderForElectrons.lostHitPenalty = 90.
TrajectoryBuilderForElectrons.alwaysUseInvalidHits = True
TrajectoryBuilderForElectrons.TTRHBuilder = 'WithTrackAngle'
TrajectoryBuilderForElectrons.updator = 'KFUpdator'




# CKFTrackCandidateMaker
from RecoTracker.CkfPattern.CkfTrackCandidates_cff import *
import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
electronCkfTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone()
electronCkfTrackCandidates.src = cms.InputTag('electronMergedSeeds')
electronCkfTrackCandidates.TrajectoryBuilder = 'TrajectoryBuilderForElectrons'
electronCkfTrackCandidates.TrajectoryCleaner = 'TrajectoryCleanerBySharedHits'
electronCkfTrackCandidates.NavigationSchool = 'SimpleNavigationSchool'
electronCkfTrackCandidates.RedundantSeedCleaner = 'CachingSeedCleanerBySharedInput'







# "backward" propagator for electrons
from TrackingTools.GsfTracking.bwdGsfElectronPropagator_cff import *
# "forward" propagator for electrons
from TrackingTools.GsfTracking.fwdGsfElectronPropagator_cff import *
# TrajectoryFilter









