import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
SeedGeneratorParameters = cms.PSet(
    EtaBound = cms.double(2.0),
    ComponentName = cms.string('TSGSmart'),
    PixelPairGeneratorSet = cms.PSet(
        ComponentName = cms.string('TSGFromOrderedHits'),
        OrderedHitsFactoryPSet = cms.PSet(
            ComponentName = cms.string('StandardHitPairGenerator'),
            SeedingLayers = cms.string('PixelLayerPairs')
        ),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    PixelTripletGeneratorSet = cms.PSet(
        ComponentName = cms.string('TSGFromOrderedHits'),
        OrderedHitsFactoryPSet = cms.PSet(
            ComponentName = cms.string('StandardHitTripletGenerator'),
            SeedingLayers = cms.string('PixelLayerTriplets'),
            GeneratorPSet = cms.PSet(
                PixelTripletHLTGenerator
            )
        ),
        TTRHBuilder = cms.string('WithTrackAngle')
    ),
    MixedGeneratorSet = cms.PSet(
        ComponentName = cms.string('TSGFromOrderedHits'),
        OrderedHitsFactoryPSet = cms.PSet(
            ComponentName = cms.string('StandardHitPairGenerator'),
            SeedingLayers = cms.string('MixedLayerPairs')
        ),
        TTRHBuilder = cms.string('WithTrackAngle')
    )
)


