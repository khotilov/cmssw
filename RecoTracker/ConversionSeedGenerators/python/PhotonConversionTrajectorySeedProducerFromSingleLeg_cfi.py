import FWCore.ParameterSet.Config as cms

photonConvTrajSeedFromSingleLeg  = cms.EDProducer("PhotonConversionTrajectorySeedProducerFromSingleLeg",
                                                  TrackRefitter        = cms.InputTag('TrackRefitter',''),
                                                  primaryVerticesTag   = cms.InputTag("offlinePrimaryVertices"),
                                                  beamSpotInputTag     = cms.InputTag("offlineBeamSpot"),
                                                  newSeedCandidates    = cms.string("convSeedCandidates"),
                                                  xcheckSeedCandidates = cms.string("xcheckSeedCandidates"),
                                                  vtxMinDoF            = cms.double(4),
                                                  maxDZSigmas          = cms.double(10.),
                                                  maxNumSelVtx         = cms.uint32(2),
                                                  applyTkVtxConstraint = cms.bool(True),
                                                  
                                                  DoxcheckSeedCandidates = cms.bool(False),
                                                  OrderedHitsFactoryPSet = cms.PSet(
                                                      maxHitPairsPerTrackAndGenerator = cms.uint32(10),
                                                      maxElement = cms.uint32(10000),
                                                      SeedingLayers = cms.string('convLayerPairs')
                                                      ),
                                                  SeedComparitorPSet = cms.PSet(
                                                      ComponentName = cms.string('none')
                                                      ),
                                                  ClusterCheckPSet = cms.PSet(
                                                      PixelClusterCollectionLabel = cms.InputTag("convClusters"),
                                                      MaxNumberOfCosmicClusters = cms.uint32(50000),
                                                      doClusterCheck = cms.bool(True),
                                                      ClusterCollectionLabel = cms.InputTag("convClusters"),
                                                      MaxNumberOfPixelClusters = cms.uint32(10000)
                                                      ),
                                                  RegionFactoryPSet = cms.PSet(
                                                      RegionPSet = cms.PSet( precise = cms.bool(True),
                                                                             beamSpot = cms.InputTag("offlineBeamSpot"),
                                                                             originRadius = cms.double(3.0),
                                                                             ptMin = cms.double(0.2),
                                                                             originHalfLength = cms.double(12.0)
                                                                             ),
                                                      ComponentName = cms.string('GlobalRegionProducerFromBeamSpot')
                                                      ),
                                                  SeedCreatorPSet = cms.PSet(
                                                      ComponentName = cms.string('SeedForPhotonConversion1Leg'),
                                                      SeedMomentumForBOFF = cms.double(5.0),
                                                      propagator = cms.string('PropagatorWithMaterial')
                                                      )
                                                  )
