###############################################################
# 
# Configuration blocks for the TrajectoryFactories inheriting 
# from TrajectoryFactoryBase.
# Include this file and do e.g.
# TrajectoryFactory = cms.PSet( ReferenceTrajectoryFactory)
# 
###############################################################

import FWCore.ParameterSet.Config as cms

###############################################################
#
# Common to all TrajectoryFactories
#
###############################################################
TrajectoryFactoryBase = cms.PSet(
    PropagationDirection = cms.string('alongMomentum'), ## or "oppositeToMomentum" or "anyDirection"
    MaterialEffects = cms.string('Combined'), ## or "MultipleScattering" or "EnergyLoss" or "None"
    UseProjectedHits = cms.bool(True), ## if false, projected hits are skipped
    UseInvalidHits = cms.bool(False), ## if false, invalid hits are skipped
    UseHitWithoutDet = cms.bool(True) ## if false, RecHits that are not attached to GeomDets are skipped
)

###############################################################
#
# ReferenceTrajectoryFactory
#
###############################################################
ReferenceTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase,
    ParticleMass = cms.double(0.10565836),
    TrajectoryFactoryName = cms.string('ReferenceTrajectoryFactory')
)

###############################################################
#
# BzeroReferenceTrajectoryFactory
#
###############################################################
BzeroReferenceTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase,
    ParticleMass = cms.double(0.10565836),
    TrajectoryFactoryName = cms.string('BzeroReferenceTrajectoryFactory'),
    MomentumEstimate = cms.double(2.0)
)

###############################################################
#
# DualTrajectoryFactory
#
###############################################################
DualTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase,
    ParticleMass = cms.double(0.10565836),
    TrajectoryFactoryName = cms.string('DualTrajectoryFactory')
)

###############################################################
#
# DualBzeroTrajectoryFactory
#
###############################################################
DualBzeroTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase,
    ParticleMass = cms.double(0.10565836),
    TrajectoryFactoryName = cms.string('DualBzeroTrajectoryFactory'),
    MomentumEstimate = cms.double(2.0)
)

###############################################################
#
# TwoBodyDecayReferenceTrajectoryFactory
#
###############################################################
TwoBodyDecayTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase,
    NSigmaCut = cms.double(100.0),
    BeamSpot = cms.PSet(
        VarYY = cms.double(2.25e-06),
        VarXX = cms.double(2.25e-06),
        VarXY = cms.double(0.0),
        VarYZ = cms.double(0.0),
        MeanX = cms.double(0.0),
        MeanY = cms.double(0.0),
        MeanZ = cms.double(0.0),
        VarXZ = cms.double(0.0),
        VarZZ = cms.double(28.09)
    ),
    ParticleProperties = cms.PSet(
        PrimaryMass = cms.double(91.1876),
        PrimaryWidth = cms.double(2.4952),
        SecondaryMass = cms.double(0.105658)
    ),
    ConstructTsosWithErrors = cms.bool(False),
    UseRefittedState = cms.bool(True),
    EstimatorParameters = cms.PSet(
        MaxIterationDifference = cms.untracked.double(0.01),
        RobustificationConstant = cms.untracked.double(1.0),
        MaxIterations = cms.untracked.int32(100),
        UseInvariantMass = cms.untracked.bool(True)
    ),
    TrajectoryFactoryName = cms.string('TwoBodyDecayTrajectoryFactory')
)

###############################################################
#
# CombinedTrajectoryFactory using an instance of TwoBodyDecayTrajectoryFactory
# and ReferenceTrajectoryFactory, taking the first successful.
#
###############################################################
CombinedTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase, # will not be used!
    TrajectoryFactoryName = cms.string('CombinedTrajectoryFactory'),
    # look for PSets called TwoBody and Reference:
    TrajectoryFactoryNames = cms.vstring(
        'TwoBodyDecayTrajectoryFactory,TwoBody',  # look for PSet called TwoBody
        'ReferenceTrajectoryFactory,Reference'),  # look for PSet called Reference
    useAllFactories = cms.bool(False),
    # now one PSet for each of the configured trajectories:
    TwoBody = cms.PSet(
        TwoBodyDecayTrajectoryFactory
    ),
    Reference = cms.PSet(
        ReferenceTrajectoryFactory
    )
)
###############################################################
#
# CombinedTrajectoryFactory using two instances of BzeroReferenceTrajectoryFactory,
# one propagating alongMomentum, one oppositeToMomentum.
#
###############################################################
# First a helper object, where I'd like to do:
#BwdBzeroReferenceTrajectoryFactory = BzeroReferenceTrajectoryFactory.clone(PropagationDirection = 'oppositeToMomentum')
# Since there is no clone in cms.PSet (yet?), but clone is needed for python that works by reference, 
# take solution from https://hypernews.cern.ch/HyperNews/CMS/get/swDevelopment/1890/1.html:
import copy
BwdBzeroReferenceTrajectoryFactory = copy.deepcopy(BzeroReferenceTrajectoryFactory)
BwdBzeroReferenceTrajectoryFactory.PropagationDirection = 'oppositeToMomentum'
# now the PSet
CombinedFwdBwdBzeroTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase, # will not be used!
    TrajectoryFactoryName = cms.string('CombinedTrajectoryFactory'),

    TrajectoryFactoryNames = cms.vstring(
        'BzeroReferenceTrajectoryFactory,FwdBzero',  # look for PSet called FwdBzero
        'BzeroReferenceTrajectoryFactory,BwdBzero'), # look for PSet called BwdBzero
    useAllFactories = cms.bool(True),
    
    # now one PSet for each of the configured trajectories:
    FwdBzero = cms.PSet(BzeroReferenceTrajectoryFactory),
    BwdBzero = cms.PSet(BwdBzeroReferenceTrajectoryFactory)
)

###############################################################
#
# CombinedTrajectoryFactory using three ReferenceTrajectories:
# - two instances of BzeroReferenceTrajectoryFactory,
#   one propagating alongMomentum, one oppositeToMomentum,
# - a DualBzeroTrajectory to start in the middle.
#
###############################################################
CombinedFwdBwdDualBzeroTrajectoryFactory = cms.PSet(
    TrajectoryFactoryBase, # will not be used!
    TrajectoryFactoryName = cms.string('CombinedTrajectoryFactory'),

    TrajectoryFactoryNames = cms.vstring(
        'BzeroReferenceTrajectoryFactory,FwdBzero',  # look for PSet called FwdBzero
        'BzeroReferenceTrajectoryFactory,BwdBzero',  # look for PSet called BwdBzero
    	'DualBzeroTrajectoryFactory,DualBzero'),     # look for PSet called DualBzero
    useAllFactories = cms.bool(True),

    # now one PSet for each of the configured trajectories:
    FwdBzero  = cms.PSet(BzeroReferenceTrajectoryFactory),
    BwdBzero  = cms.PSet(BwdBzeroReferenceTrajectoryFactory), # defined above for CombinedFwdBwdBzeroTrajectoryFactory
    DualBzero = cms.PSet(DualBzeroTrajectoryFactory)
)

###############################################################
#
# DualKalmanFactory
#
###############################################################
DualKalmanFactory = cms.PSet(
    TrajectoryFactoryBase,
    ParticleMass = cms.double(0.10565836),
    TrajectoryFactoryName = cms.string('DualKalmanFactory'),
    ResidualMethod = cms.int32(2) # 1: unbiased residuals, 2: pulls
)

