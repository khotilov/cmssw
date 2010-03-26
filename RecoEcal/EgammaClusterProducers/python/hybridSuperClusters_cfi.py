import FWCore.ParameterSet.Config as cms

from RecoEcal.EgammaClusterProducers.ecalRecHitFlags_cfi import *

# Producer for Hybrid BasicClusters and SuperClusters
hybridSuperClusters = cms.EDProducer("HybridClusterProducer",
    eThreshA = cms.double(0.003),
    # seed thresold for dominos
    eseed = cms.double(0.35),
    posCalc_x0 = cms.double(0.89),
    # output collections
    clustershapecollection = cms.string(''),
    shapeAssociation = cms.string('hybridShapeAssoc'),
    # if e1x3 larger than ewing use 1x5
    # e.g. always build 1x5
    ewing = cms.double(0.0),
    # clustering parameters
    #
    # threshold on seed RecHits
    HybridBarrelSeedThr = cms.double(1.0),
    dynamicPhiRoad = cms.bool(False),
    basicclusterCollection = cms.string('hybridBarrelBasicClusters'),
    posCalc_w0 = cms.double(4.2),
    # phi road parameters
    step = cms.int32(17),
    eThreshB = cms.double(0.1),
    posCalc_t0 = cms.double(7.4),
    debugLevel = cms.string('INFO'),
    dynamicEThresh = cms.bool(False),
    # domino thresholds
    ethresh = cms.double(0.1),
    superclusterCollection = cms.string(''),
    posCalc_logweight = cms.bool(True),
    ecalhitcollection = cms.string('EcalRecHitsEB'),
    # input collection
    ecalhitproducer = cms.string('ecalRecHit'),
    # recHit flags to be excluded from seeding
    RecHitFlagToBeExcluded = cms.vint32(
        ecalRecHitFlag_kFaultyHardware,
        ecalRecHitFlag_kPoorCalib,
        ecalRecHitFlag_kSaturated,
        ecalRecHitFlag_kLeadingEdgeRecovered,
        ecalRecHitFlag_kNeighboursRecovered,
        ecalRecHitFlag_kTowerRecovered,
        ecalRecHitFlag_kDead
    ),
    RecHitSeverityToBeExcluded = cms.vint32(999),
    severityRecHitThreshold = cms.double(0.08),
    severitySpikeId = cms.int32(1),
    severitySpikeThreshold = cms.double(0)
 ) 
