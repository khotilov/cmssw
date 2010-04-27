import FWCore.ParameterSet.Config as cms
 
# import the whole HLT menu
#from HLTrigger.Configuration.HLT_8E29_cff import *
#from HLTrigger.Configuration.HLT_1E31_cff import *
from HLTrigger.Configuration.HLT_FULL_cff import *

# create the jetMET HLT reco path
DoHLTJets = cms.Path(HLTBeginSequence + 
    HLTRecoJetMETSequence +
    HLTDoJet30HTRecoSequence
    #HLTDoHTRecoSequence
)
DoHLTJetsU = cms.Path(HLTBeginSequence +
    HLTRecoJetSequenceU +
    hltMet +
    HLTDoJet15UHTRecoSequence
)
# create the muon HLT reco path
DoHltMuon = cms.Path(
    HLTBeginSequence +
    HLTL2muonrecoSequence + 
    HLTL2muonisorecoSequence + 
    HLTL3muonrecoSequence + 
    HLTL3muonisorecoSequence +
    HLTEndSequence )

# create the Egamma HLT reco paths
DoHLTPhoton = cms.Path( 
    HLTBeginSequence + 
    HLTDoRegionalEgammaEcalSequence + 
    HLTL1IsolatedEcalClustersSequence + 
    HLTL1NonIsolatedEcalClustersSequence + 
    hltL1IsoRecoEcalCandidate + 
    hltL1NonIsoRecoEcalCandidate + 
    hltL1IsolatedPhotonEcalIsol + 
    hltL1NonIsolatedPhotonEcalIsol + 
    HLTDoLocalHcalWithoutHOSequence + 
    hltL1IsolatedPhotonHcalIsol + 
    hltL1NonIsolatedPhotonHcalIsol + 
    HLTDoLocalTrackerSequence + 
    HLTL1IsoEgammaRegionalRecoTrackerSequence + 
    HLTL1NonIsoEgammaRegionalRecoTrackerSequence + 
    hltL1IsoPhotonHollowTrackIsol + 
    hltL1NonIsoPhotonHollowTrackIsol )

##DoHLTElectron = cms.Path( 
##    HLTBeginSequence + 
##    HLTDoRegionalEgammaEcalSequence + 
##    HLTL1IsolatedEcalClustersSequence + 
##    HLTL1NonIsolatedEcalClustersSequence + 
##    hltL1IsoRecoEcalCandidate + 
##    hltL1NonIsoRecoEcalCandidate + 
##    HLTDoLocalHcalWithoutHOSequence + 
##    hltL1IsolatedElectronHcalIsol + 
##    hltL1NonIsolatedElectronHcalIsol + 
##    HLTDoLocalPixelSequence + 
##    HLTDoLocalStripSequence + 
##    HLTPixelMatchElectronL1IsoSequence + 
##    HLTPixelMatchElectronL1NonIsoSequence + 
##    HLTPixelMatchElectronL1IsoTrackingSequence + 
##    HLTPixelMatchElectronL1NonIsoTrackingSequence + 
##    HLTL1IsoElectronsRegionalRecoTrackerSequence + 
##    HLTL1NonIsoElectronsRegionalRecoTrackerSequence + 
##    hltL1IsoElectronTrackIsol + 
##    hltL1NonIsoElectronTrackIsol )

##DoHLTElectronStartUpWindows = cms.Path( 
##    HLTBeginSequence + 
##    HLTDoRegionalEgammaEcalSequence + 
##    HLTL1IsolatedEcalClustersSequence + 
##    HLTL1NonIsolatedEcalClustersSequence + 
##    hltL1IsoRecoEcalCandidate + 
##    hltL1NonIsoRecoEcalCandidate + 
##    HLTDoLocalHcalWithoutHOSequence + 
##    hltL1IsolatedElectronHcalIsol + 
##    hltL1NonIsolatedElectronHcalIsol + 
##    HLTDoLocalPixelSequence + 
##    HLTDoLocalStripSequence + 
##    hltL1IsoStartUpElectronPixelSeeds + 
##    hltL1NonIsoStartUpElectronPixelSeeds + 
##    HLTPixelMatchStartUpElectronL1IsoTrackingSequence + 
##    HLTPixelMatchElectronL1NonIsoStartUpTrackingSequence + 
##    HLTL1IsoStartUpElectronsRegionalRecoTrackerSequence + 
##    HLTL1NonIsoStartUpElectronsRegionalRecoTrackerSequence + 
##    hltL1IsoStartUpElectronTrackIsol + 
##    hltL1NonIsoStartupElectronTrackIsol )

## LW ele RECO L1ISO
hltCkfL1IsoLWTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1IsoLargeWindowElectronPixelSeeds" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
                                  
)

hltCtfL1IsoLW = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltCtfL1IsoLW" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltCkfL1IsoLWTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltPixelMatchElectronsL1IsoLW = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1IsoLW" ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" )
)

## SiStrip ele RECO L1ISO
hltCkfL1IsoSSTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1IsoSiStripElectronPixelSeeds" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
                                  
)

hltCtfL1IsoSS = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltCtfL1IsoSS" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltCkfL1IsoSSTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltPixelMatchElectronsL1IsoSS = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1IsoSS" ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" )
)

## LW ele RECO L1NonISO
hltCkfL1NonIsoLWTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1NonIsoLargeWindowElectronPixelSeeds" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)

hltCtfL1NonIsoLW = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltCtfL1NonIsoLW" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltCkfL1NonIsoLWTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltPixelMatchElectronsL1NonIsoLW = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1NonIsoLW" ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" )
)

## SiStrip ele RECO L1NonISO
hltCkfL1NonIsoSSTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1NonIsoSiStripElectronPixelSeeds" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)

hltCtfL1NonIsoSS = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltCtfL1NonIsoSS" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltCkfL1NonIsoSSTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltPixelMatchElectronsL1NonIsoSS = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1NonIsoSS" ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" )
)

## LW ele TrackIso L1Iso
hltL1IsoLWEleRegPSG = cms.EDProducer( "EgammaHLTRegionalPixelSeedGeneratorProducers",
    ptMin = cms.double( 1.5 ),
    vertexZ = cms.double( 0.0 ),
    originRadius = cms.double( 0.02 ),
    originHalfLength = cms.double( 0.5 ),
    deltaEtaRegion = cms.double( 0.3 ),
    deltaPhiRegion = cms.double( 0.3 ),
    candTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    candTagEle = cms.InputTag( "hltPixelMatchElectronsL1IsoLW" ),
    UseZInVertex = cms.bool( True ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitPairGenerator" ),
      SeedingLayers = cms.string( "MixedLayerPairs" )
    )
)

hltL1IsoLWEleRegioCkfTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1IsoLWEleRegPSG" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)

hltL1IsoLWEleRegioCTF = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltL1IsoLWEleRegioCTF" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL1IsoLWEleRegioCkfTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)

## SiStrip ele TrackIso L1Iso
hltL1IsoSSEleRegPSG = cms.EDProducer( "EgammaHLTRegionalPixelSeedGeneratorProducers",
    ptMin = cms.double( 1.5 ),
    vertexZ = cms.double( 0.0 ),
    originRadius = cms.double( 0.02 ),
    originHalfLength = cms.double( 0.5 ),
    deltaEtaRegion = cms.double( 0.3 ),
    deltaPhiRegion = cms.double( 0.3 ),
    candTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    candTagEle = cms.InputTag( "hltPixelMatchElectronsL1IsoSS" ),
    UseZInVertex = cms.bool( True ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitPairGenerator" ),
      SeedingLayers = cms.string( "MixedLayerPairs" )
    )
)

hltL1IsoSSEleRegioCkfTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1IsoSSEleRegPSG" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)

hltL1IsoSSEleRegioCTF = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltL1IsoSSEleRegioCTF" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL1IsoSSEleRegioCkfTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)

## LW ele TrackIso L1NonIso
hltL1NonIsoLWEleRegPSG = cms.EDProducer( "EgammaHLTRegionalPixelSeedGeneratorProducers",
    ptMin = cms.double( 1.5 ),
    vertexZ = cms.double( 0.0 ),
    originRadius = cms.double( 0.02 ),
    originHalfLength = cms.double( 0.5 ),
    deltaEtaRegion = cms.double( 0.3 ),
    deltaPhiRegion = cms.double( 0.3 ),
    candTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    candTagEle = cms.InputTag( "hltPixelMatchElectronsL1NonIsoLW" ),
    UseZInVertex = cms.bool( True ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitPairGenerator" ),
      SeedingLayers = cms.string( "MixedLayerPairs" )
    )
)
hltL1NonIsoLWEleRegioCkfTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1NonIsoLWEleRegPSG" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)
hltL1NonIsoLWEleRegioCTF = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltL1NonIsoLWEleRegioCTF" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL1NonIsoLWEleRegioCkfTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)

## SiStrip ele TrackIso L1NonIso
hltL1NonIsoSSEleRegPSG = cms.EDProducer( "EgammaHLTRegionalPixelSeedGeneratorProducers",
    ptMin = cms.double( 1.5 ),
    vertexZ = cms.double( 0.0 ),
    originRadius = cms.double( 0.02 ),
    originHalfLength = cms.double( 0.5 ),
    deltaEtaRegion = cms.double( 0.3 ),
    deltaPhiRegion = cms.double( 0.3 ),
    candTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    candTagEle = cms.InputTag( "hltPixelMatchElectronsL1NonIsoSS" ),
    UseZInVertex = cms.bool( True ),
    BSProducer = cms.InputTag( "hltOfflineBeamSpot" ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitPairGenerator" ),
      SeedingLayers = cms.string( "MixedLayerPairs" )
    )
)
hltL1NonIsoSSEleRegioCkfTC = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltL1NonIsoSSEleRegPSG" ),
    TrajectoryBuilder = cms.string( "CkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)
hltL1NonIsoSSEleRegioCTF = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltL1NonIsoSSEleRegioCTF" ),
    Fitter = cms.string( "hltKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL1NonIsoSSEleRegioCkfTC" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)

## track iso for LW
hltL1IsoLWEleTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1IsoLW" ),
    trackProducer = cms.InputTag( "hltL1IsoLWEleRegioCTF" ),
    egTrkIsoPtMin = cms.double( 1.5 ),
    egTrkIsoConeSize = cms.double( 0.2 ),
    egTrkIsoZSpan = cms.double( 0.1 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.02 )
)

hltL1NonIsoLWEleTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIsoLW" ),
    trackProducer = cms.InputTag( "hltL1NonIsoLWEleRegioCTF" ),
    egTrkIsoPtMin = cms.double( 1.5 ),
    egTrkIsoConeSize = cms.double( 0.2 ),
    egTrkIsoZSpan = cms.double( 0.1 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.02 )
)

## track iso for SiStrip
hltL1IsoSSEleTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1IsoSS" ),
    trackProducer = cms.InputTag( "hltL1IsoSSEleRegioCTF" ),
    egTrkIsoPtMin = cms.double( 1.5 ),
    egTrkIsoConeSize = cms.double( 0.2 ),
    egTrkIsoZSpan = cms.double( 0.1 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.02 )
)

hltL1NonIsoSSEleTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIsoSS" ),
    trackProducer = cms.InputTag( "hltL1NonIsoSSEleRegioCTF" ),
    egTrkIsoPtMin = cms.double( 1.5 ),
    egTrkIsoConeSize = cms.double( 0.2 ),
    egTrkIsoZSpan = cms.double( 0.1 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.02 )
)

DoLWTracking = cms.Sequence(
    hltCkfL1IsoLWTC +
    hltCtfL1IsoLW +
    hltPixelMatchElectronsL1IsoLW +
    hltCkfL1NonIsoLWTC +
    hltCtfL1NonIsoLW +
    hltPixelMatchElectronsL1NonIsoLW +
    hltL1IsoLWEleRegPSG +
    hltL1IsoLWEleRegioCkfTC +
    hltL1IsoLWEleRegioCTF +
    hltL1IsoLWEleTrackIsol +
    hltL1NonIsoLWEleRegPSG +
    hltL1NonIsoLWEleRegioCkfTC +
    hltL1NonIsoLWEleRegioCTF +
    hltL1NonIsoLWEleTrackIsol     
    )

DoSSTracking = cms.Sequence(
    hltCkfL1IsoSSTC +
    hltCtfL1IsoSS +
    hltPixelMatchElectronsL1IsoSS +
    hltCkfL1NonIsoSSTC +
    hltCtfL1NonIsoSS +
    hltPixelMatchElectronsL1NonIsoSS +
    hltL1IsoSSEleRegPSG +
    hltL1IsoSSEleRegioCkfTC +
    hltL1IsoSSEleRegioCTF +
    hltL1IsoSSEleTrackIsol +
    hltL1NonIsoSSEleRegPSG +
    hltL1NonIsoSSEleRegioCkfTC +
    hltL1NonIsoSSEleRegioCTF +
    hltL1NonIsoSSEleTrackIsol     
    )

DoHLTElectronStartUpWindows = cms.Path( 
    HLTBeginSequence + 
    HLTDoRegionalEgammaEcalSequence + 
    HLTL1IsolatedEcalClustersSequence + 
    HLTL1NonIsolatedEcalClustersSequence + 
    hltL1IsoRecoEcalCandidate + 
    hltL1NonIsoRecoEcalCandidate + 
    HLTDoLocalHcalWithoutHOSequence + 
    hltL1IsolatedElectronHcalIsol + 
    hltL1NonIsolatedElectronHcalIsol + 
    HLTDoLocalPixelSequence + 
    HLTDoLocalStripSequence + 
    hltL1IsoStartUpElectronPixelSeeds + 
    hltL1NonIsoStartUpElectronPixelSeeds + 
    HLTPixelMatchElectronL1IsoTrackingSequence + 
    HLTPixelMatchElectronL1NonIsoTrackingSequence + 
    HLTL1IsoElectronsRegionalRecoTrackerSequence + 
    HLTL1NonIsoElectronsRegionalRecoTrackerSequence + 
    hltL1IsoElectronTrackIsol + 
    hltL1NonIsoElectronTrackIsol )

DoHLTElectronLargeWindows = cms.Path( 
    HLTBeginSequence + 
    HLTDoRegionalEgammaEcalSequence + 
    HLTL1IsolatedEcalClustersSequence + 
    HLTL1NonIsolatedEcalClustersSequence + 
    hltL1IsoRecoEcalCandidate + 
    hltL1NonIsoRecoEcalCandidate + 
    HLTDoLocalHcalWithoutHOSequence + 
    hltL1IsolatedElectronHcalIsol + 
    hltL1NonIsolatedElectronHcalIsol + 
    HLTDoLocalPixelSequence + 
    HLTDoLocalStripSequence +
    hltL1IsoLargeWindowElectronPixelSeeds +
    hltL1NonIsoLargeWindowElectronPixelSeeds +
    DoLWTracking
    )

DoHLTElectronSiStrip = cms.Path( 
    HLTBeginSequence + 
    HLTDoRegionalEgammaEcalSequence + 
    HLTL1IsolatedEcalClustersSequence + 
    HLTL1NonIsolatedEcalClustersSequence + 
    hltL1IsoRecoEcalCandidate + 
    hltL1NonIsoRecoEcalCandidate + 
    HLTDoLocalHcalWithoutHOSequence + 
    hltL1IsolatedElectronHcalIsol + 
    hltL1NonIsolatedElectronHcalIsol + 
    HLTDoLocalPixelSequence + 
    HLTDoLocalStripSequence +
    hltL1IsoSiStripElectronPixelSeeds +
    hltL1NonIsoSiStripElectronPixelSeeds +
    DoSSTracking
    )

# create the tau HLT reco path
from HLTrigger.HLTanalyzers.OpenHLT_Tau_cff import *
DoHLTTau = cms.Path(HLTBeginSequence +
                    openhltTauPrescaler +
                    OpenHLTCaloTausCreatorSequence +
                    openhltL2TauJets +
                    openhltL2TauIsolationProducer +
#                    openhltL2TauRelaxingIsolationSelector +
                    HLTDoLocalPixelSequence +
                    HLTRecopixelvertexingSequence +
                    OpenHLTL25TauTrackReconstructionSequence +
                    OpenHLTL25TauTrackIsolation +
                    TauOpenHLT+
                    HLTEndSequence)


# create the b-jet HLT paths
from HLTrigger.HLTanalyzers.OpenHLT_BJet_cff import *
DoHLTBTag = cms.Path(
    HLTBeginSequence +
    HLTBCommonL2recoSequence +
    OpenHLTBLifetimeL25recoSequence +
    OpenHLTBSoftmuonL25recoSequence +
    OpenHLTBLifetimeL3recoSequence +
    OpenHLTBLifetimeL3recoSequenceStartup +
    OpenHLTBSoftmuonL3recoSequence +
    HLTEndSequence )


DoHLTAlCaPi0Eta1E31 = cms.Path(
    HLTBeginSequence +
    hltL1sAlCaEcalPi0Eta1E31 +
    #hltPreAlCaEcalPi01E31 +
    #HLTDoRegionalPi0EtaESSequence +
    #HLTDoRegionalPi0EtaEcalSequence +
    HLTDoRegionalPi0EtaSequence +
    HLTEndSequence )

DoHLTAlCaPi0Eta8E29 = cms.Path(
    HLTBeginSequence +
    hltL1sAlCaEcalPi0Eta8E29 +
    #hltPreAlCaEcalPi08E29 +
    #HLTDoRegionalPi0EtaESSequence +
    #HLTDoRegionalPi0EtaEcalSequence +
    HLTDoRegionalPi0EtaSequence +
    HLTEndSequence )


DoHLTAlCaECALPhiSym = cms.Path(
    HLTBeginSequence +
    hltEcalRawToRecHitFacility + hltESRawToRecHitFacility + hltEcalRegionalRestFEDs + hltEcalRecHitAll +
##    hltAlCaPhiSymStream +
##    HLTDoLocalHcalSequence +
    HLTEndSequence )

#hltIsolPixelTrackProd1E31.MaxVtxDXYSeed = cms.double(101)
#hltIsolPixelTrackL2Filter1E31.MaxPtNearby = cms.double(3.0)
#hltIsolPixelTrackL2Filter1E31.MinPtTrack = cms.double(3.0)

#DoHLTIsoTrack = cms.Path(
#    HLTBeginSequence +
#    hltL1sIsoTrack1E31 +
    # hltPreIsoTrack1E31 +
#    HLTL2HcalIsolTrackSequence +
#    hltIsolPixelTrackProd1E31 +
#    hltIsolPixelTrackL2Filter1E31 +
#    HLTDoLocalStripSequence +
#    hltHITPixelPairSeedGenerator1E31 +
#    hltHITPixelTripletSeedGenerator1E31 +
#    hltHITSeedCombiner1E31 +
#    hltHITCkfTrackCandidates1E31 +
#    hltHITCtfWithMaterialTracks1E31 +
#    hltHITIPTCorrector1E31 +
#    HLTEndSequence)

#DoHLTIsoTrack8E29 = cms.Path(
#    HLTBeginSequence +
#    hltL1sIsoTrack8E29 +
    # hltPreIsoTrack8E29 +
#    HLTL2HcalIsolTrackSequence +
#    hltIsolPixelTrackProd8E29 +
#    hltIsolPixelTrackL2Filter8E29 +
#    HLTDoLocalStripSequence +
#    hltHITPixelPairSeedGenerator8E29 +
#    hltHITPixelTripletSeedGenerator8E29 +
#    hltHITSeedCombiner8E29 +
#    hltHITCkfTrackCandidates8E29 +
#    hltHITCtfWithMaterialTracks8E29 +
#    hltHITIPTCorrector8E29 +
#    HLTEndSequence)


DoHLTMinBiasPixelTracks = cms.Path(
    HLTBeginSequence + 
#    hltPreMinBiasPixelSingleTrack +
    hltSiPixelDigis + 
    hltSiPixelClusters + 
    hltSiPixelRecHits +
    hltPixelTracksForMinBias +
    hltPixelCands)

# create the Onia HLT reco path
oniaTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
    ComponentName = cms.string( "oniaTrajectoryBuilder" ),
    updator = cms.string( "KFUpdator" ),
    propagatorAlong = cms.string( "PropagatorWithMaterial" ),
    propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
    estimator = cms.string( "Chi2" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    MeasurementTrackerName = cms.string( "" ),
    trajectoryFilterName = cms.string( "oniaTrajectoryFilter" ),
    maxCand = cms.int32( 1 ),
    lostHitPenalty = cms.double( 30.0 ),
    intermediateCleaning = cms.bool( True ),
    alwaysUseInvalidHits = cms.bool( False ),
    appendToDataLabel = cms.string( "" )
)
oniaTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
    ComponentName = cms.string( "oniaTrajectoryFilter" ),
    appendToDataLabel = cms.string( "" ),
    filterPset = cms.PSet(
      chargeSignificance = cms.double( -1.0 ),
      minPt = cms.double( 1.0 ),
      minHitsMinPt = cms.int32( 3 ),
      ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
      maxLostHits = cms.int32( 1 ),
      maxNumberOfHits = cms.int32( 8 ),
      maxConsecLostHits = cms.int32( 1 ),
      nSigmaMinPt = cms.double( 5.0 ),
      minimumNumberOfHits = cms.int32( 5 )
    )
)
hltPixelTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltPixelTracks" ),
    particleType = cms.string( "mu-" )
)
HLTOniaPixelTrackSequence = cms.Sequence( HLTDoLocalPixelSequence
                                                   + hltPixelTracks + hltPixelTrackCands )
hltOniaPixelTrackSelector = cms.EDProducer("QuarkoniaTrackSelector",
    muonCandidates = cms.InputTag("hltL3MuonCandidates"),
    tracks = cms.InputTag("hltPixelTracks"),
    MinMasses = cms.vdouble(2, 8),
    MaxMasses = cms.vdouble(5, 12),
    checkCharge = cms.bool(False),
    MinTrackPt = cms.double(0.),
    MinTrackP = cms.double(2.5),
    MaxTrackEta = cms.double(999.)
)
hltOniaPixelTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltOniaPixelTrackSelector" ),
    particleType = cms.string( "mu-" )
)
hltOniaSeeds = cms.EDProducer("SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag("hltOniaPixelTrackSelector"),
    TTRHBuilder = cms.string("WithTrackAngle"),
    useProtoTrackKinematics = cms.bool(False)
)
hltOniaCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltOniaSeeds" ),
    TrajectoryBuilder = cms.string( "oniaTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "TrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet(
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" ),
      numberMeasurementsForFit = cms.int32(4)
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 100000 )
)
hltOniaCtfTracks = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltOniaCtfTracks" ),
    Fitter = cms.string( "FittingSmootherRK" ),
    Propagator = cms.string( "RungeKuttaTrackerPropagator" ),
    src = cms.InputTag( "hltOniaCkfTrackCandidates" ),
    beamSpot = cms.InputTag( "hltOfflineBeamSpot" ),
    TTRHBuilder = cms.string( "WithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltOniaCtfTrackCands = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltOniaCtfTracks" ),
    particleType = cms.string( "mu-" )
)
HLTOniaTrackSequence = cms.Sequence( HLTDoLocalStripSequence + hltOniaSeeds
                                               + hltOniaCkfTrackCandidates
                                               + hltOniaCtfTracks + hltOniaCtfTrackCands )

DoHLT_Onia_1E31 = cms.Path( HLTBeginSequence
                          + HLTL2muonrecoSequence
                          + HLTL3muonrecoSequence
                          + HLTOniaPixelTrackSequence + hltPixelTrackCands
                          + hltOniaPixelTrackSelector + hltOniaPixelTrackCands
                          + HLTOniaTrackSequence
                          + HLTEndSequence )

