# /dev/CMSSW_4_2_0/GRun/V54 (CMSSW_4_2_0_HLT2)

import FWCore.ParameterSet.Config as cms
from FastSimulation.HighLevelTrigger.HLTSetup_cff import *


HLTConfigVersion = cms.PSet(
  tableName = cms.string('/dev/CMSSW_4_2_0/GRun/V54')
)

hltESSAK5CaloL2L3 = cms.ESSource( "JetCorrectionServiceChain",
  appendToDataLabel = cms.string( "" ),
  correctors = cms.vstring( 'hltESSL2RelativeCorrectionService',
    'hltESSL3AbsoluteCorrectionService' ),
  label = cms.string( "hltESSAK5CaloL2L3" )
)
hltESSBTagRecord = cms.ESSource( "EmptyESSource",
  recordName = cms.string( "JetTagComputerRecord" ),
  iovIsRunNotTime = cms.bool( True ),
  appendToDataLabel = cms.string( "" ),
  firstValid = cms.vuint32( 1 )
)
hltESSHcalSeverityLevel = cms.ESSource( "EmptyESSource",
  recordName = cms.string( "HcalSeverityLevelComputerRcd" ),
  iovIsRunNotTime = cms.bool( True ),
  appendToDataLabel = cms.string( "" ),
  firstValid = cms.vuint32( 1 )
)
hltESSL2RelativeCorrectionService = cms.ESSource( "LXXXCorrectionService",
  appendToDataLabel = cms.string( "" ),
  level = cms.string( "L2Relative" ),
  algorithm = cms.string( "AK5Calo" ),
  section = cms.string( "" ),
  era = cms.string( "" ),
  useCondDB = cms.untracked.bool( True )
)
hltESSL3AbsoluteCorrectionService = cms.ESSource( "LXXXCorrectionService",
  appendToDataLabel = cms.string( "" ),
  level = cms.string( "L3Absolute" ),
  algorithm = cms.string( "AK5Calo" ),
  section = cms.string( "" ),
  era = cms.string( "" ),
  useCondDB = cms.untracked.bool( True )
)

AnyDirectionAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  ComponentName = cms.string( "AnyDirectionAnalyticalPropagator" ),
  PropagationDirection = cms.string( "anyDirection" ),
  MaxDPhi = cms.double( 1.6 ),
  appendToDataLabel = cms.string( "" )
)
CaloTowerGeometryFromDBEP = cms.ESProducer( "CaloTowerGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
EcalBarrelGeometryFromDBEP = cms.ESProducer( "EcalBarrelGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
EcalEndcapGeometryFromDBEP = cms.ESProducer( "EcalEndcapGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
EcalPreshowerGeometryFromDBEP = cms.ESProducer( "EcalPreshowerGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
EcalUnpackerWorkerESProducer = cms.ESProducer( "EcalUnpackerWorkerESProducer",
  ComponentName = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  DCCDataUnpacker = cms.PSet( 
    tccUnpacking = cms.bool( False ),
    orderedDCCIdList = cms.vint32( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54 ),
    srpUnpacking = cms.bool( False ),
    syncCheck = cms.bool( False ),
    feIdCheck = cms.bool( True ),
    headerUnpacking = cms.bool( True ),
    orderedFedList = cms.vint32( 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654 ),
    feUnpacking = cms.bool( True ),
    forceKeepFRData = cms.bool( False ),
    memUnpacking = cms.bool( True )
  ),
  ElectronicsMapper = cms.PSet( 
    numbXtalTSamples = cms.uint32( 10 ),
    numbTriggerTSamples = cms.uint32( 1 )
  ),
  UncalibRHAlgo = cms.PSet(  Type = cms.string( "EcalUncalibRecHitWorkerWeights" ) ),
  CalibRHAlgo = cms.PSet( 
    flagsMapDBReco = cms.vint32( 0, 0, 0, 0, 4, -1, -1, -1, 4, 4, 6, 6, 6, 7, 8 ),
    Type = cms.string( "EcalRecHitWorkerSimple" ),
    killDeadChannels = cms.bool( True ),
    ChannelStatusToBeExcluded = cms.vint32( 10, 11, 12, 13, 14, 78, 142 ),
    laserCorrection = cms.bool( False )
  )
)
HcalGeometryFromDBEP = cms.ESProducer( "HcalGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
XMLFromDBSource = cms.ESProducer( "XMLIdealGeometryESProducer",
  rootDDName = cms.string( "cms:OCMS" ),
  label = cms.string( "Extended" ),
  appendToDataLabel = cms.string( "" )
)
ZdcGeometryFromDBEP = cms.ESProducer( "ZdcGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
caloDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "CaloDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
cosmicsNavigationSchoolESProducer = cms.ESProducer( "NavigationSchoolESProducer",
  ComponentName = cms.string( "CosmicNavigationSchool" ),
  appendToDataLabel = cms.string( "" )
)
ecalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "EcalDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.02 ),
  nEta = cms.int32( 300 ),
  nPhi = cms.int32( 360 ),
  includeBadChambers = cms.bool( False )
)
ecalSeverityLevel = cms.ESProducer( "EcalSeverityLevelESProducer",
  appendToDataLabel = cms.string( "" ),
  flagMask = cms.vuint32( 1, 34, 896, 4, 49152, 6232 ),
  dbstatusMask = cms.vuint32( 1, 2046, 0, 0, 0, 64512 ),
  timeThresh = cms.double( 2.0 )
)
hcalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HcalDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
hcalRecAlgos = cms.ESProducer( "HcalRecAlgoESProducer",
  SeverityLevels = cms.VPSet( 
    cms.PSet(  RecHitFlags = cms.vstring(  ),
      ChannelStatus = cms.vstring(  ),
      Level = cms.int32( 0 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring(  ),
      ChannelStatus = cms.vstring( 'HcalCellCaloTowerProb' ),
      Level = cms.int32( 1 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring( 'HSCP_R1R2',
  'HSCP_FracLeader',
  'HSCP_OuterEnergy',
  'HSCP_ExpFit',
  'ADCSaturationBit' ),
      ChannelStatus = cms.vstring(  ),
      Level = cms.int32( 5 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring( 'HBHEHpdHitMultiplicity',
  'HFDigiTime',
  'HBHEPulseShape',
  'HOBit',
  'HFInTimeWindow',
  'ZDCBit',
  'CalibrationBit',
  'TimingErrorBit' ),
      ChannelStatus = cms.vstring(  ),
      Level = cms.int32( 8 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring( 'HFLongShort',
  'HFS8S1Ratio',
  'HFPET' ),
      ChannelStatus = cms.vstring(  ),
      Level = cms.int32( 11 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring(  ),
      ChannelStatus = cms.vstring( 'HcalCellCaloTowerMask' ),
      Level = cms.int32( 12 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring(  ),
      ChannelStatus = cms.vstring( 'HcalCellHot' ),
      Level = cms.int32( 15 )
    ),
    cms.PSet(  RecHitFlags = cms.vstring(  ),
      ChannelStatus = cms.vstring( 'HcalCellOff',
        'HcalCellDead' ),
      Level = cms.int32( 20 )
    )
  ),
  RecoveredRecHitBits = cms.vstring( 'TimingAddedBit',
    'TimingSubtractedBit' ),
  appendToDataLabel = cms.string( "" ),
  DropChannelStatusBits = cms.vstring( 'HcalCellMask',
    'HcalCellOff',
    'HcalCellDead' )
)
hltESPAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  ComponentName = cms.string( "hltESPAnalyticalPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  MaxDPhi = cms.double( 1.6 ),
  appendToDataLabel = cms.string( "" )
)
hltESPChi2EstimatorForRefit = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2EstimatorForRefit" ),
  MaxChi2 = cms.double( 100000.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPChi2MeasurementEstimator = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator" ),
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPCkf3HitTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPCkf3HitTrajectoryBuilder" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPCkf3HitTrajectoryFilter" ),
  maxCand = cms.int32( 5 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( True ),
  appendToDataLabel = cms.string( "" )
)
hltESPCkf3HitTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPCkf3HitTrajectoryFilter" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.9 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( -1 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 3 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
hltESPCkfTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPCkfTrajectoryBuilder" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPCkfTrajectoryFilter" ),
  maxCand = cms.int32( 5 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( True ),
  appendToDataLabel = cms.string( "" )
)
hltESPCkfTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPCkfTrajectoryFilter" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.9 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( -1 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 5 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
hltESPDummyDetLayerGeometry = cms.ESProducer( "DetLayerGeometryESProducer",
  ComponentName = cms.string( "hltESPDummyDetLayerGeometry" ),
  appendToDataLabel = cms.string( "" )
)
hltESPESUnpackerWorker = cms.ESProducer( "ESUnpackerWorkerESProducer",
  ComponentName = cms.string( "hltESPESUnpackerWorker" ),
  appendToDataLabel = cms.string( "" ),
  DCCDataUnpacker = cms.PSet(  LookupTable = cms.FileInPath( "EventFilter/ESDigiToRaw/data/ES_lookup_table.dat" ) ),
  RHAlgo = cms.PSet( 
    ESRecoAlgo = cms.int32( 0 ),
    Type = cms.string( "ESRecHitWorker" )
  )
)
hltESPEcalRegionCablingESProducer = cms.ESProducer( "EcalRegionCablingESProducer",
  appendToDataLabel = cms.string( "" ),
  esMapping = cms.PSet(  LookupTable = cms.FileInPath( "EventFilter/ESDigiToRaw/data/ES_lookup_table.dat" ) )
)
hltESPFastSteppingHelixPropagatorAny = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
  PropagationDirection = cms.string( "anyDirection" ),
  useInTeslaFromMagField = cms.bool( False ),
  SetVBFPointer = cms.bool( False ),
  useMagVolumes = cms.bool( True ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  ApplyRadX0Correction = cms.bool( True ),
  AssumeNoMaterial = cms.bool( False ),
  NoErrorPropagation = cms.bool( False ),
  debug = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  returnTangentPlane = cms.bool( True ),
  sendLogWarning = cms.bool( False ),
  useTuningForL2Speed = cms.bool( True ),
  useEndcapShiftsInZ = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPFastSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useInTeslaFromMagField = cms.bool( False ),
  SetVBFPointer = cms.bool( False ),
  useMagVolumes = cms.bool( True ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  ApplyRadX0Correction = cms.bool( True ),
  AssumeNoMaterial = cms.bool( False ),
  NoErrorPropagation = cms.bool( False ),
  debug = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  returnTangentPlane = cms.bool( True ),
  sendLogWarning = cms.bool( False ),
  useTuningForL2Speed = cms.bool( True ),
  useEndcapShiftsInZ = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPFittingSmootherRK = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPFittingSmootherRK" ),
  Fitter = cms.string( "hltESPTrajectoryFitterRK" ),
  Smoother = cms.string( "hltESPTrajectorySmootherRK" ),
  EstimateCut = cms.double( -1.0 ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  MinNumberOfHits = cms.int32( 5 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPHIPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPHIPixelLayerPairs" ),
  layerList = cms.vstring( 'BPix1+BPix2',
    'BPix1+BPix3',
    'BPix2+BPix3',
    'BPix1+FPix1_pos',
    'BPix1+FPix1_neg',
    'BPix1+FPix2_pos',
    'BPix1+FPix2_neg',
    'BPix2+FPix1_pos',
    'BPix2+FPix1_neg',
    'BPix2+FPix2_pos',
    'BPix2+FPix2_neg',
    'FPix1_pos+FPix2_pos',
    'FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPHIPixelLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPHIPixelLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPHITTRHBuilderWithoutRefit = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPHITTRHBuilderWithoutRefit" ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "Fake" ),
  Matcher = cms.string( "Fake" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmoother" ),
  Fitter = cms.string( "hltESPKFTrajectoryFitter" ),
  Smoother = cms.string( "hltESPKFTrajectorySmoother" ),
  EstimateCut = cms.double( -1.0 ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  MinNumberOfHits = cms.int32( 5 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFFittingSmootherForL2Muon = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmootherForL2Muon" ),
  Fitter = cms.string( "hltESPKFTrajectoryFitterForL2Muon" ),
  Smoother = cms.string( "hltESPKFTrajectorySmootherForL2Muon" ),
  EstimateCut = cms.double( -1.0 ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  MinNumberOfHits = cms.int32( 5 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( False ),
  NoInvalidHitsBeginEnd = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitter" ),
  Propagator = cms.string( "PropagatorWithMaterial" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFTrajectoryFitterForL2Muon = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitterForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmoother" ),
  Propagator = cms.string( "PropagatorWithMaterial" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFTrajectorySmootherForL2Muon = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFTrajectorySmootherForMuonTrackLoader = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
  Propagator = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPKFUpdator = cms.ESProducer( "KFUpdatorESProducer",
  ComponentName = cms.string( "hltESPKFUpdator" ),
  appendToDataLabel = cms.string( "" )
)
hltESPL3MuKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
  Propagator = cms.string( "hltESPSmartPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltESPMeasurementTracker" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  StripCPE = cms.string( "StripCPEfromTrackAngle" ),
  HitMatcher = cms.string( "StandardMatcher" ),
  Regional = cms.bool( True ),
  OnDemand = cms.bool( True ),
  UsePixelModuleQualityDB = cms.bool( True ),
  DebugPixelModuleQualityDB = cms.untracked.bool( False ),
  UsePixelROCQualityDB = cms.bool( True ),
  DebugPixelROCQualityDB = cms.untracked.bool( False ),
  UseStripModuleQualityDB = cms.bool( True ),
  DebugStripModuleQualityDB = cms.untracked.bool( False ),
  UseStripAPVFiberQualityDB = cms.bool( True ),
  DebugStripAPVFiberQualityDB = cms.untracked.bool( False ),
  MaskBadAPVFibers = cms.bool( True ),
  UseStripStripQualityDB = cms.bool( True ),
  DebugStripStripQualityDB = cms.untracked.bool( False ),
  SiStripQualityLabel = cms.string( "" ),
  switchOffPixelsIfEmpty = cms.bool( True ),
  pixelClusterProducer = cms.string( "hltSiPixelClusters" ),
  skipClusters = cms.InputTag( "" ),
  stripClusterProducer = cms.string( "hltSiStripClusters" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag(  ),
  badStripCuts = cms.PSet( 
    TID = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TEC = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TIB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    )
  )
)
hltESPMixedLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPMixedLayerPairs" ),
  layerList = cms.vstring( 'BPix1+BPix2',
    'BPix1+BPix3',
    'BPix2+BPix3',
    'BPix1+FPix1_pos',
    'BPix1+FPix1_neg',
    'BPix1+FPix2_pos',
    'BPix1+FPix2_neg',
    'BPix2+FPix1_pos',
    'BPix2+FPix1_neg',
    'BPix2+FPix2_pos',
    'BPix2+FPix2_neg',
    'FPix1_pos+FPix2_pos',
    'FPix1_neg+FPix2_neg',
    'FPix2_pos+TEC1_pos',
    'FPix2_pos+TEC2_pos',
    'TEC1_pos+TEC2_pos',
    'TEC2_pos+TEC3_pos',
    'FPix2_neg+TEC1_neg',
    'FPix2_neg+TEC2_neg',
    'TEC1_neg+TEC2_neg',
    'TEC2_neg+TEC3_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet( 
    useRingSlector = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    minRing = cms.int32( 1 ),
    maxRing = cms.int32( 1 )
  )
)
hltESPMuTrackJpsiTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPMuTrackJpsiTrajectoryBuilder" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPMuTrackJpsiTrajectoryFilter" ),
  maxCand = cms.int32( 1 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPMuTrackJpsiTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPMuTrackJpsiTrajectoryFilter" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 1.0 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 8 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 5 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
hltESPMuonCkfTrajectoryBuilder = cms.ESProducer( "MuonCkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPMuonCkfTrajectoryBuilder" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  propagatorProximity = cms.string( "SteppingHelixPropagatorAny" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPMuonCkfTrajectoryFilter" ),
  useSeedLayer = cms.bool( False ),
  rescaleErrorIfFail = cms.double( 1.0 ),
  deltaEta = cms.double( 0.1 ),
  deltaPhi = cms.double( 0.1 ),
  appendToDataLabel = cms.string( "" ),
  maxCand = cms.int32( 5 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( False ),
  alwaysUseInvalidHits = cms.bool( True )
)
hltESPMuonCkfTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPMuonCkfTrajectoryFilter" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.9 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( -1 ),
    maxConsecLostHits = cms.int32( 1 ),
    chargeSignificance = cms.double( -1.0 ),
    nSigmaMinPt = cms.double( 5.0 ),
    minimumNumberOfHits = cms.int32( 5 )
  )
)
hltESPMuonTransientTrackingRecHitBuilder = cms.ESProducer( "MuonTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
  appendToDataLabel = cms.string( "" )
)
hltESPPixelCPEGeneric = cms.ESProducer( "PixelCPEGenericESProducer",
  ComponentName = cms.string( "hltESPPixelCPEGeneric" ),
  eff_charge_cut_lowX = cms.double( 0.0 ),
  eff_charge_cut_lowY = cms.double( 0.0 ),
  eff_charge_cut_highX = cms.double( 1.0 ),
  eff_charge_cut_highY = cms.double( 1.0 ),
  size_cutX = cms.double( 3.0 ),
  size_cutY = cms.double( 3.0 ),
  EdgeClusterErrorX = cms.double( 50.0 ),
  EdgeClusterErrorY = cms.double( 85.0 ),
  inflate_errors = cms.bool( False ),
  inflate_all_errors_no_trk_angle = cms.bool( False ),
  UseErrorsFromTemplates = cms.bool( True ),
  TruncatePixelCharge = cms.bool( True ),
  IrradiationBiasCorrection = cms.bool( False ),
  DoCosmics = cms.bool( False ),
  LoadTemplatesFromDB = cms.bool( True ),
  appendToDataLabel = cms.string( "" ),
  TanLorentzAnglePerTesla = cms.double( 0.106 ),
  PixelErrorParametrization = cms.string( "NOTcmsim" ),
  Alpha2Order = cms.bool( True ),
  ClusterProbComputationFlag = cms.int32( 0 )
)
hltESPPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerPairs" ),
  layerList = cms.vstring( 'BPix1+BPix2',
    'BPix1+BPix3',
    'BPix2+BPix3',
    'BPix1+FPix1_pos',
    'BPix1+FPix1_neg',
    'BPix1+FPix2_pos',
    'BPix1+FPix2_neg',
    'BPix2+FPix1_pos',
    'BPix2+FPix1_neg',
    'BPix2+FPix2_pos',
    'BPix2+FPix2_neg',
    'FPix1_pos+FPix2_pos',
    'FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPPixelLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPPixelLayerTripletsHITHB = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTripletsHITHB" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPPixelLayerTripletsHITHE = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTripletsHITHE" ),
  layerList = cms.vstring( 'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0060 ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  FPix = cms.PSet( 
    hitErrorRZ = cms.double( 0.0036 ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    useErrorsFromParam = cms.bool( True )
  ),
  TEC = cms.PSet(  )
)
hltESPPromptTrackCountingESProducer = cms.ESProducer( "PromptTrackCountingESProducer",
  appendToDataLabel = cms.string( "" ),
  impactParameterType = cms.int32( 0 ),
  maximumDistanceToJetAxis = cms.double( 999999.0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 999999.0 ),
  maxImpactParameterSig = cms.double( 999999.0 ),
  trackQualityClass = cms.string( "any" ),
  nthTrack = cms.int32( -1 ),
  maxImpactParameter = cms.double( 0.03 ),
  deltaRmin = cms.double( 0.0 )
)
hltESPRungeKuttaTrackerPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Mass = cms.double( 0.105 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( True ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPSiStripRegionConnectivity = cms.ESProducer( "SiStripRegionConnectivity",
  EtaDivisions = cms.untracked.uint32( 20 ),
  PhiDivisions = cms.untracked.uint32( 20 ),
  EtaMax = cms.untracked.double( 2.5 )
)
hltESPSmartPropagator = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
  appendToDataLabel = cms.string( "" )
)
hltESPSmartPropagatorAny = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAny" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  appendToDataLabel = cms.string( "" )
)
hltESPSmartPropagatorAnyOpposite = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterialOpposite" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  appendToDataLabel = cms.string( "" )
)
hltESPSmartPropagatorOpposite = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterialOpposite" ),
  MuonPropagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
  appendToDataLabel = cms.string( "" )
)
hltESPSoftLeptonByDistance = cms.ESProducer( "LeptonTaggerByDistanceESProducer",
  appendToDataLabel = cms.string( "" ),
  distance = cms.double( 0.5 )
)
hltESPSoftLeptonByPt = cms.ESProducer( "LeptonTaggerByPtESProducer",
  appendToDataLabel = cms.string( "" ),
  ipSign = cms.string( "any" )
)
hltESPSteppingHelixPropagatorAlong = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  useInTeslaFromMagField = cms.bool( False ),
  SetVBFPointer = cms.bool( False ),
  useMagVolumes = cms.bool( True ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  ApplyRadX0Correction = cms.bool( True ),
  AssumeNoMaterial = cms.bool( False ),
  NoErrorPropagation = cms.bool( False ),
  debug = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  returnTangentPlane = cms.bool( True ),
  sendLogWarning = cms.bool( False ),
  useTuningForL2Speed = cms.bool( False ),
  useEndcapShiftsInZ = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  useInTeslaFromMagField = cms.bool( False ),
  SetVBFPointer = cms.bool( False ),
  useMagVolumes = cms.bool( True ),
  VBFName = cms.string( "VolumeBasedMagneticField" ),
  ApplyRadX0Correction = cms.bool( True ),
  AssumeNoMaterial = cms.bool( False ),
  NoErrorPropagation = cms.bool( False ),
  debug = cms.bool( False ),
  useMatVolumes = cms.bool( True ),
  useIsYokeFlag = cms.bool( True ),
  returnTangentPlane = cms.bool( True ),
  sendLogWarning = cms.bool( False ),
  useTuningForL2Speed = cms.bool( False ),
  useEndcapShiftsInZ = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  appendToDataLabel = cms.string( "" )
)
hltESPStraightLinePropagator = cms.ESProducer( "StraightLinePropagatorESProducer",
  ComponentName = cms.string( "hltESPStraightLinePropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  appendToDataLabel = cms.string( "" )
)
hltESPTTRHBWithTrackAngle = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBWithTrackAngle" ),
  StripCPE = cms.string( "StripCPEfromTrackAngle" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPTTRHBuilderPixelOnly = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderPixelOnly" ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPTrackCounting3D2nd = cms.ESProducer( "TrackCountingESProducer",
  appendToDataLabel = cms.string( "" ),
  nthTrack = cms.int32( 2 ),
  impactParameterType = cms.int32( 0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 5.0 ),
  maximumDistanceToJetAxis = cms.double( 0.07 ),
  trackQualityClass = cms.string( "any" )
)
hltESPTrajectoryBuilderL3 = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPTrajectoryBuilderL3" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPTrajectoryFilterL3" ),
  maxCand = cms.int32( 5 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  appendToDataLabel = cms.string( "" ),
  fractionShared = cms.double( 0.5 ),
  allowSharedFirstHit = cms.bool( False )
)
hltESPTrajectoryFilterL3 = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFilterL3" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.9 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 7 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 5 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
hltESPTrajectoryFitterRK = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFitterRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPTrajectorySmootherRK = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPTrajectorySmootherRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
hltESPbJetRegionalTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPbJetRegionalTrajectoryBuilder" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPbJetRegionalTrajectoryFilter" ),
  maxCand = cms.int32( 1 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
hltESPbJetRegionalTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPbJetRegionalTrajectoryFilter" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 1.0 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 8 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 5 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
hoDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HODetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 30 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
muonDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "MuonDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.125 ),
  nEta = cms.int32( 48 ),
  nPhi = cms.int32( 48 ),
  includeBadChambers = cms.bool( False )
)
preshowerDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "PreshowerDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.1 ),
  nEta = cms.int32( 60 ),
  nPhi = cms.int32( 30 ),
  includeBadChambers = cms.bool( False )
)
siPixelTemplateDBObjectESProducer = cms.ESProducer( "SiPixelTemplateDBObjectESProducer",
  appendToDataLabel = cms.string( "" )
)
sistripconn = cms.ESProducer( "SiStripConnectivity" )

DQMStore = cms.Service( "DQMStore",
)
DTDataIntegrityTask = cms.Service( "DTDataIntegrityTask",
  getSCInfo = cms.untracked.bool( True ),
  fedIntegrityFolder = cms.untracked.string( "DT/FEDIntegrity_EvF" ),
  processingMode = cms.untracked.string( "HLT" )
)

hltGetRaw = cms.EDAnalyzer( "HLTGetRaw",
    RawDataCollection = cms.InputTag( "rawDataCollector" )
)
hltBoolFalse = cms.EDFilter( "HLTBool",
    result = cms.bool( False )
)
hltL1sL1BscMinBiasORBptxPlusANDMinus = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BscMinBiasOR_BptxPlusANDMinus" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreActivityEcalSC7 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHybridSuperClustersActivity = cms.EDProducer( "HybridClusterProducer",
    debugLevel = cms.string( "ERROR" ),
    basicclusterCollection = cms.string( "hybridBarrelBasicClusters" ),
    superclusterCollection = cms.string( "" ),
    ecalhitproducer = cms.string( "hltEcalRecHitAll" ),
    ecalhitcollection = cms.string( "EcalRecHitsEB" ),
    HybridBarrelSeedThr = cms.double( 1.0 ),
    step = cms.int32( 17 ),
    ethresh = cms.double( 0.1 ),
    eseed = cms.double( 0.35 ),
    ewing = cms.double( 0.0 ),
    dynamicEThresh = cms.bool( False ),
    eThreshA = cms.double( 0.0030 ),
    eThreshB = cms.double( 0.1 ),
    excludeFlagged = cms.bool( False ),
    dynamicPhiRoad = cms.bool( False ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    RecHitSeverityToBeExcluded = cms.vint32( 4 ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    ),
    bremRecoveryPset = cms.PSet(  ),
    severitySpikeId = cms.int32( 2 ),
    severitySpikeThreshold = cms.double( 0.99 ),
    severityRecHitThreshold = cms.double( 4.0 )
)
hltCorrectedHybridSuperClustersActivity = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    rawSuperClusterProducer = cms.InputTag( "hltHybridSuperClustersActivity" ),
    superClusterAlgo = cms.string( "Hybrid" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 5.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 1.1 ),
      fEtEtaVec = cms.vdouble( 0.0, 1.00121, -0.63672, 0.0, 0.0, 0.0, 0.5655, 6.457, 0.5081, 8.0, 1.023, -0.00181 ),
      brLinearHighThr = cms.double( 8.0 ),
      fBremVec = cms.vdouble( -0.04382, 0.1169, 0.9267, -9.413E-4, 1.419 )
    ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet(  )
)
hltMulti5x5BasicClustersActivity = cms.EDProducer( "Multi5x5ClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    barrelHitProducer = cms.string( "hltEcalRecHitAll" ),
    endcapHitProducer = cms.string( "hltEcalRecHitAll" ),
    barrelHitCollection = cms.string( "EcalRecHitsEB" ),
    endcapHitCollection = cms.string( "EcalRecHitsEE" ),
    doEndcap = cms.bool( True ),
    doBarrel = cms.bool( False ),
    barrelClusterCollection = cms.string( "multi5x5BarrelBasicClusters" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    IslandBarrelSeedThr = cms.double( 0.5 ),
    IslandEndcapSeedThr = cms.double( 0.18 ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    )
)
hltMulti5x5SuperClustersActivity = cms.EDProducer( "Multi5x5SuperClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    endcapClusterProducer = cms.string( "hltMulti5x5BasicClustersActivity" ),
    barrelClusterProducer = cms.string( "hltMulti5x5BasicClustersActivity" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    barrelClusterCollection = cms.string( "multi5x5BarrelBasicClusters" ),
    endcapSuperclusterCollection = cms.string( "multi5x5EndcapSuperClusters" ),
    barrelSuperclusterCollection = cms.string( "multi5x5BarrelSuperClusters" ),
    doBarrel = cms.bool( False ),
    doEndcaps = cms.bool( True ),
    barrelEtaSearchRoad = cms.double( 0.06 ),
    barrelPhiSearchRoad = cms.double( 0.8 ),
    endcapEtaSearchRoad = cms.double( 0.14 ),
    endcapPhiSearchRoad = cms.double( 0.6 ),
    seedTransverseEnergyThreshold = cms.double( 1.0 ),
    dynamicPhiRoad = cms.bool( False ),
    bremRecoveryPset = cms.PSet( 
      barrel = cms.PSet( 
        cryVec = cms.vint32( 16, 13, 11, 10, 9, 8, 7, 6, 5, 4, 3 ),
        cryMin = cms.int32( 2 ),
        etVec = cms.vdouble( 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 45.0, 55.0, 135.0, 195.0, 225.0 )
      ),
      endcap = cms.PSet( 
        a = cms.double( 47.85 ),
        c = cms.double( 0.1201 ),
        b = cms.double( 108.8 )
      )
    )
)
hltMulti5x5SuperClustersWithPreshowerActivity = cms.EDProducer( "PreshowerClusterProducer",
    preshRecHitProducer = cms.InputTag( 'hltESRecHitAll','EcalRecHitsES' ),
    endcapSClusterProducer = cms.InputTag( 'hltMulti5x5SuperClustersActivity','multi5x5EndcapSuperClusters' ),
    preshClusterCollectionX = cms.string( "preshowerXClusters" ),
    preshClusterCollectionY = cms.string( "preshowerYClusters" ),
    preshNclust = cms.int32( 4 ),
    etThresh = cms.double( 0.0 ),
    assocSClusterCollection = cms.string( "" ),
    preshStripEnergyCut = cms.double( 0.0 ),
    preshSeededNstrip = cms.int32( 15 ),
    preshClusterEnergyCut = cms.double( 0.0 ),
    debugLevel = cms.string( "ERROR" )
)
hltCorrectedMulti5x5SuperClustersWithPreshowerActivity = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' ),
    rawSuperClusterProducer = cms.InputTag( "hltMulti5x5SuperClustersWithPreshowerActivity" ),
    superClusterAlgo = cms.string( "Multi5x5" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 5.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet(  ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 0.9 ),
      fEtEtaVec = cms.vdouble( 1.0, -0.4386, -32.38, 0.6372, 15.67, -0.0928, -2.462, 1.138, 20.93 ),
      brLinearHighThr = cms.double( 6.0 ),
      fBremVec = cms.vdouble( -0.05228, 0.08738, 0.9508, 0.002677, 1.221 )
    )
)
hltRecoEcalSuperClusterActivityCandidate = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( "hltCorrectedHybridSuperClustersActivity" ),
    scIslandEndcapProducer = cms.InputTag( "hltCorrectedMulti5x5SuperClustersWithPreshowerActivity" ),
    recoEcalCandidateCollection = cms.string( "" )
)
hltEcalActivitySuperClusterWrapper = cms.EDFilter( "HLTEgammaTriggerFilterObjectWrapper",
    candIsolatedTag = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    candNonIsolatedTag = cms.InputTag( "none" ),
    doIsolated = cms.bool( True )
)
hltEgammaSelectEcalSuperClustersActivityFilterSC7 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 7.0 ),
    etcutEE = cms.double( 7.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "none" )
)
hltL1sL1SingleJet16 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet16" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1SingleJet16 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sL1SingleJet36 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet36" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1SingleJet36 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTowerMakerForAll = cms.EDProducer( "CaloTowersCreator",
    EBThreshold = cms.double( 0.07 ),
    EEThreshold = cms.double( 0.3 ),
    UseEtEBTreshold = cms.bool( False ),
    UseEtEETreshold = cms.bool( False ),
    UseSymEBTreshold = cms.bool( False ),
    UseSymEETreshold = cms.bool( False ),
    HcalThreshold = cms.double( -1000.0 ),
    HBThreshold = cms.double( 0.7 ),
    HESThreshold = cms.double( 0.8 ),
    HEDThreshold = cms.double( 0.8 ),
    HOThreshold0 = cms.double( 3.5 ),
    HOThresholdPlus1 = cms.double( 3.5 ),
    HOThresholdMinus1 = cms.double( 3.5 ),
    HOThresholdPlus2 = cms.double( 3.5 ),
    HOThresholdMinus2 = cms.double( 3.5 ),
    HF1Threshold = cms.double( 0.5 ),
    HF2Threshold = cms.double( 0.85 ),
    EBWeight = cms.double( 1.0 ),
    EEWeight = cms.double( 1.0 ),
    HBWeight = cms.double( 1.0 ),
    HESWeight = cms.double( 1.0 ),
    HEDWeight = cms.double( 1.0 ),
    HOWeight = cms.double( 1.0E-99 ),
    HF1Weight = cms.double( 1.0 ),
    HF2Weight = cms.double( 1.0 ),
    EcutTower = cms.double( -1000.0 ),
    EBSumThreshold = cms.double( 0.2 ),
    EESumThreshold = cms.double( 0.45 ),
    UseHO = cms.bool( False ),
    MomConstrMethod = cms.int32( 1 ),
    MomHBDepth = cms.double( 0.2 ),
    MomHEDepth = cms.double( 0.4 ),
    MomEBDepth = cms.double( 0.3 ),
    MomEEDepth = cms.double( 0.0 ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    hoInput = cms.InputTag( "hltHoreco" ),
    hfInput = cms.InputTag( "hltHfreco" ),
    AllowMissingInputs = cms.bool( False ),
    HcalAcceptSeverityLevel = cms.uint32( 9 ),
    EcalAcceptSeverityLevel = cms.uint32( 3 ),
    UseHcalRecoveredHits = cms.bool( False ),
    UseEcalRecoveredHits = cms.bool( False ),
    UseRejectedHitsOnly = cms.bool( False ),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    EcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    UseRejectedRecoveredHcalHits = cms.bool( False ),
    UseRejectedRecoveredEcalHits = cms.bool( False ),
    EBGrid = cms.vdouble(  ),
    EBWeights = cms.vdouble(  ),
    EEGrid = cms.vdouble(  ),
    EEWeights = cms.vdouble(  ),
    HBGrid = cms.vdouble(  ),
    HBWeights = cms.vdouble(  ),
    HESGrid = cms.vdouble(  ),
    HESWeights = cms.vdouble(  ),
    HEDGrid = cms.vdouble(  ),
    HEDWeights = cms.vdouble(  ),
    HOGrid = cms.vdouble(  ),
    HOWeights = cms.vdouble(  ),
    HF1Grid = cms.vdouble(  ),
    HF1Weights = cms.vdouble(  ),
    HF2Grid = cms.vdouble(  ),
    HF2Weights = cms.vdouble(  ),
    ecalInputs = cms.VInputTag( 'hltEcalRecHitAll:EcalRecHitsEB','hltEcalRecHitAll:EcalRecHitsEE' )
)
hltAntiKT5CaloJets = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.5 ),
    src = cms.InputTag( "hltTowerMakerForAll" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltAntiKT5L2L3CorrCaloJets = cms.EDProducer( "CaloJetCorrectionProducer",
    src = cms.InputTag( "hltAntiKT5CaloJets" ),
    verbose = cms.untracked.bool( False ),
    alias = cms.untracked.string( "JetCorJetAntiKT5" ),
    correctors = cms.vstring( 'hltESSAK5CaloL2L3' )
)
hltJetIDPassedCorrJets = cms.EDProducer( "HLTJetIDProducer",
    jetsInput = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    min_EMF = cms.double( 1.0E-6 ),
    max_EMF = cms.double( 999.0 ),
    min_N90 = cms.int32( 2 )
)
hltSingleJet30 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet60 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTowerMakerForJets = cms.EDProducer( "CaloTowersCreator",
    EBThreshold = cms.double( 0.07 ),
    EEThreshold = cms.double( 0.3 ),
    UseEtEBTreshold = cms.bool( False ),
    UseEtEETreshold = cms.bool( False ),
    UseSymEBTreshold = cms.bool( False ),
    UseSymEETreshold = cms.bool( False ),
    HcalThreshold = cms.double( -1000.0 ),
    HBThreshold = cms.double( 0.7 ),
    HESThreshold = cms.double( 0.8 ),
    HEDThreshold = cms.double( 0.8 ),
    HOThreshold0 = cms.double( 3.5 ),
    HOThresholdPlus1 = cms.double( 3.5 ),
    HOThresholdMinus1 = cms.double( 3.5 ),
    HOThresholdPlus2 = cms.double( 3.5 ),
    HOThresholdMinus2 = cms.double( 3.5 ),
    HF1Threshold = cms.double( 0.5 ),
    HF2Threshold = cms.double( 0.85 ),
    EBWeight = cms.double( 1.0 ),
    EEWeight = cms.double( 1.0 ),
    HBWeight = cms.double( 1.0 ),
    HESWeight = cms.double( 1.0 ),
    HEDWeight = cms.double( 1.0 ),
    HOWeight = cms.double( 1.0E-99 ),
    HF1Weight = cms.double( 1.0 ),
    HF2Weight = cms.double( 1.0 ),
    EcutTower = cms.double( -1000.0 ),
    EBSumThreshold = cms.double( 0.2 ),
    EESumThreshold = cms.double( 0.45 ),
    UseHO = cms.bool( False ),
    MomConstrMethod = cms.int32( 1 ),
    MomHBDepth = cms.double( 0.2 ),
    MomHEDepth = cms.double( 0.4 ),
    MomEBDepth = cms.double( 0.3 ),
    MomEEDepth = cms.double( 0.0 ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    hoInput = cms.InputTag( "hltHoreco" ),
    hfInput = cms.InputTag( "hltHfreco" ),
    AllowMissingInputs = cms.bool( False ),
    HcalAcceptSeverityLevel = cms.uint32( 9 ),
    EcalAcceptSeverityLevel = cms.uint32( 3 ),
    UseHcalRecoveredHits = cms.bool( False ),
    UseEcalRecoveredHits = cms.bool( False ),
    UseRejectedHitsOnly = cms.bool( False ),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    EcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    UseRejectedRecoveredHcalHits = cms.bool( False ),
    UseRejectedRecoveredEcalHits = cms.bool( False ),
    EBGrid = cms.vdouble(  ),
    EBWeights = cms.vdouble(  ),
    EEGrid = cms.vdouble(  ),
    EEWeights = cms.vdouble(  ),
    HBGrid = cms.vdouble(  ),
    HBWeights = cms.vdouble(  ),
    HESGrid = cms.vdouble(  ),
    HESWeights = cms.vdouble(  ),
    HEDGrid = cms.vdouble(  ),
    HEDWeights = cms.vdouble(  ),
    HOGrid = cms.vdouble(  ),
    HOWeights = cms.vdouble(  ),
    HF1Grid = cms.vdouble(  ),
    HF1Weights = cms.vdouble(  ),
    HF2Grid = cms.vdouble(  ),
    HF2Weights = cms.vdouble(  ),
    ecalInputs = cms.VInputTag( 'hltEcalRegionalJetsRecHit:EcalRecHitsEB','hltEcalRegionalJetsRecHit:EcalRecHitsEE' )
)
hltAntiKT5CaloJetsRegional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.5 ),
    src = cms.InputTag( "hltTowerMakerForJets" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltAntiKT5L2L3CorrCaloJetsRegional = cms.EDProducer( "CaloJetCorrectionProducer",
    src = cms.InputTag( "hltAntiKT5CaloJetsRegional" ),
    verbose = cms.untracked.bool( False ),
    alias = cms.untracked.string( "JetCorJetAntiKT5" ),
    correctors = cms.vstring( 'hltESSAK5CaloL2L3' )
)
hltL1MatchedJetsRegional = cms.EDProducer( "HLTJetL1MatchProducer",
    jetsInput = cms.InputTag( "hltAntiKT5L2L3CorrCaloJetsRegional" ),
    L1TauJets = cms.InputTag( 'l1extraParticles','Tau' ),
    L1CenJets = cms.InputTag( 'l1extraParticles','Central' ),
    L1ForJets = cms.InputTag( 'l1extraParticles','Forward' ),
    DeltaR = cms.double( 0.5 )
)
hltJetIDPassedJetsRegional = cms.EDProducer( "HLTJetIDProducer",
    jetsInput = cms.InputTag( "hltL1MatchedJetsRegional" ),
    min_EMF = cms.double( 1.0E-6 ),
    max_EMF = cms.double( 999.0 ),
    min_N90 = cms.int32( 2 )
)
hltSingleJet60Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 60.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1SingleJet52 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet52" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreJet80 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet80Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 80.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1SingleJet68 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet68" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreJet110 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet110Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 110.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1SingleJet92 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet92" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreJet150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet150Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 150.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet190 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet190Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 190.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet240 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet240Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 240.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet300Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 300.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet370 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet370Regional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 370.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreJet370NoJetID = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleJet370RegionalNoJetID = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltL1MatchedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 370.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreDiJetAve30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve30 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 30.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve60 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve60 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 60.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve80 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve80 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 80.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve110 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve110 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 110.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve150 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 150.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve190 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve190 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 190.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve240 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve240 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 240.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve300 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 300.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltPreDiJetAve370 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJetAve370 = cms.EDFilter( "HLTDiJetAveFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( False ),
    minPtAve = cms.double( 370.0 ),
    minPtJet3 = cms.double( 99999.0 ),
    minDphi = cms.double( -1.0 )
)
hltL1sL1DoubleForJet32EtaOpp = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleForJet32_EtaOpp" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleJet30ForwardBackward = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleJet30ForwardBackward = cms.EDFilter( "HLTForwardBackwardJetsFilter",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    minPt = cms.double( 30.0 ),
    minEta = cms.double( 3.0 ),
    maxEta = cms.double( 5.1 )
)
hltPreDoubleJet60ForwardBackward = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleJet60ForwardBackward = cms.EDFilter( "HLTForwardBackwardJetsFilter",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    minPt = cms.double( 60.0 ),
    minEta = cms.double( 3.0 ),
    maxEta = cms.double( 5.1 )
)
hltPreDoubleJet70ForwardBackward = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleJet70ForwardBackward = cms.EDFilter( "HLTForwardBackwardJetsFilter",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    minPt = cms.double( 70.0 ),
    minEta = cms.double( 3.0 ),
    maxEta = cms.double( 5.1 )
)
hltL1sL1DoubleForJet44EtaOpp = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleForJet44_EtaOpp" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleJet80ForwardBackward = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleJet80ForwardBackward = cms.EDFilter( "HLTForwardBackwardJetsFilter",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    minPt = cms.double( 80.0 ),
    minEta = cms.double( 3.0 ),
    maxEta = cms.double( 5.1 )
)
hltPreDiJet130PT130 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDijet130PT130 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 3 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 130.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 130.0, 130.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreDiJet160PT160 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDijet160PT160 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 3 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 160.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 160.0, 160.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltL1sL1ETM30 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ETM30" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreCenJet80MET65 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltCenJet80CentralRegional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 80.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltMet = cms.EDProducer( "METProducer",
    src = cms.InputTag( "hltTowerMakerForAll" ),
    InputType = cms.string( "CandidateCollection" ),
    METType = cms.string( "CaloMET" ),
    alias = cms.string( "RawCaloMET" ),
    globalThreshold = cms.double( 0.3 ),
    noHF = cms.bool( True ),
    calculateSignificance = cms.bool( False ),
    onlyFiducialParticles = cms.bool( False ),
    jets = cms.InputTag( "unused" ),
    rf_type = cms.int32( 0 ),
    correctShowerTracks = cms.bool( False ),
    HO_EtResPar = cms.vdouble( 0.0, 1.3, 0.0050 ),
    HF_EtResPar = cms.vdouble( 0.0, 1.82, 0.09 ),
    HB_PhiResPar = cms.vdouble( 0.02511 ),
    HE_PhiResPar = cms.vdouble( 0.02511 ),
    EE_EtResPar = cms.vdouble( 0.2, 0.03, 0.0050 ),
    EB_PhiResPar = cms.vdouble( 0.00502 ),
    EE_PhiResPar = cms.vdouble( 0.02511 ),
    HB_EtResPar = cms.vdouble( 0.0, 1.22, 0.05 ),
    EB_EtResPar = cms.vdouble( 0.2, 0.03, 0.0050 ),
    HF_PhiResPar = cms.vdouble( 0.05022 ),
    HE_EtResPar = cms.vdouble( 0.0, 1.3, 0.05 ),
    HO_PhiResPar = cms.vdouble( 0.02511 )
)
hltMET65 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 65.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltPreCenJet80MET80HF = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltJetIDMPassedJetsRegional = cms.EDProducer( "HLTJetIDProducer",
    jetsInput = cms.InputTag( "hltL1MatchedJetsRegional" ),
    min_EMF = cms.double( 0.01 ),
    max_EMF = cms.double( 999999.0 ),
    min_N90 = cms.int32( 2 )
)
hltCenJet80MCentralRegional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDMPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 80.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltMET80 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 80.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltMetWithHF = cms.EDProducer( "METProducer",
    src = cms.InputTag( "hltTowerMakerForAll" ),
    InputType = cms.string( "CandidateCollection" ),
    METType = cms.string( "CaloMET" ),
    alias = cms.string( "RawCaloMET" ),
    globalThreshold = cms.double( 0.3 ),
    noHF = cms.bool( False ),
    calculateSignificance = cms.bool( False ),
    onlyFiducialParticles = cms.bool( False ),
    jets = cms.InputTag( "unused" ),
    rf_type = cms.int32( 0 ),
    correctShowerTracks = cms.bool( False ),
    HO_EtResPar = cms.vdouble( 0.0, 1.3, 0.0050 ),
    HF_EtResPar = cms.vdouble( 0.0, 1.82, 0.09 ),
    HB_PhiResPar = cms.vdouble( 0.02511 ),
    HE_PhiResPar = cms.vdouble( 0.02511 ),
    EE_EtResPar = cms.vdouble( 0.2, 0.03, 0.0050 ),
    EB_PhiResPar = cms.vdouble( 0.00502 ),
    EE_PhiResPar = cms.vdouble( 0.02511 ),
    HB_EtResPar = cms.vdouble( 0.0, 1.22, 0.05 ),
    EB_EtResPar = cms.vdouble( 0.2, 0.03, 0.0050 ),
    HF_PhiResPar = cms.vdouble( 0.05022 ),
    HE_EtResPar = cms.vdouble( 0.0, 1.3, 0.05 ),
    HO_PhiResPar = cms.vdouble( 0.02511 )
)
hltMETWithHF80 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMetWithHF" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 80.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltPreCenJet80MET100 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMET100 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 100.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltPreCenJet80MET160 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMET160 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 160.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1ETM20 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ETM20" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDiJet60MET45 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJet60 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 60.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 2 )
)
hltMet45 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 45.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltPre2CenJet20MET80 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hlt2CenJet20CentralRegional = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedJetsRegional" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltPre2CenJet20BtagIPMET65 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBJetHbb = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltGetJetsfromBJetHbb = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBJetHbb" )
)
hltSelectorJetsHbb = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfromBJetHbb" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 6 )
)
hltBLifetimeL25JetsHbb = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfromBJetHbb" ),
    filter = cms.bool( False ),
    etMin = cms.double( 20.0 )
)
hltBLifetimeL25AssociatorHbb = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsHbb" ),
    tracks = cms.InputTag( "hltPixelTracks" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL25TagInfosHbb = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL25AssociatorHbb" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL25BJetTagsHbb = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL25TagInfosHbb' )
)
hltBLifetimeL25FilterHbb = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL25BJetTagsHbb" ),
    MinTag = cms.double( 0.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBLifetimeL3AssociatorHbb = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsHbb" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracksHbb" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfosHbb = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3AssociatorHbb" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTagsHbb = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfosHbb' )
)
hltBLifetimeL3FilterHbb = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTagsHbb" ),
    MinTag = cms.double( 2.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltL1sL1QuadJet20Central = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_QuadJet20_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreQuadJet40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet40Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 4 )
)
hltPreQuadJet40IsoPFTau40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet40IsoPFTau40 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 4 )
)
hltAntiKT5CaloJetsEt5 = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltAntiKT5CaloJets" ),
    filter = cms.bool( False ),
    etMin = cms.double( 5.0 )
)
hltParticleFlowRecHitECAL = cms.EDProducer( "PFRecHitProducerECAL",
    ecalRecHitsEB = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    ecalRecHitsEE = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' ),
    crossBarrelEndcapBorder = cms.bool( False ),
    timing_Cleaning = cms.bool( False ),
    topological_Cleaning = cms.bool( True ),
    thresh_Cleaning_EB = cms.double( 2.0 ),
    thresh_Cleaning_EE = cms.double( 1.0E9 ),
    verbose = cms.untracked.bool( False ),
    thresh_Barrel = cms.double( 0.08 ),
    thresh_Endcap = cms.double( 0.3 )
)
hltParticleFlowRecHitHCAL = cms.EDProducer( "PFRecHitProducerHCAL",
    hcalRecHitsHBHE = cms.InputTag( "hltHbhereco" ),
    hcalRecHitsHF = cms.InputTag( "hltHfreco" ),
    caloTowers = cms.InputTag( "hltTowerMakerForAll" ),
    thresh_HF = cms.double( 0.4 ),
    navigation_HF = cms.bool( True ),
    weight_HFem = cms.double( 1.0 ),
    weight_HFhad = cms.double( 1.0 ),
    HCAL_Calib = cms.bool( True ),
    HF_Calib = cms.bool( True ),
    HCAL_Calib_29 = cms.double( 1.35 ),
    HF_Calib_29 = cms.double( 1.07 ),
    ShortFibre_Cut = cms.double( 60.0 ),
    LongFibre_Fraction = cms.double( 0.05 ),
    LongFibre_Cut = cms.double( 120.0 ),
    ShortFibre_Fraction = cms.double( 0.01 ),
    ApplyLongShortDPG = cms.bool( False ),
    LongShortFibre_Cut = cms.double( 1.0E9 ),
    MinShortTiming_Cut = cms.double( -5.0 ),
    MaxShortTiming_Cut = cms.double( 5.0 ),
    MinLongTiming_Cut = cms.double( -5.0 ),
    MaxLongTiming_Cut = cms.double( 5.0 ),
    ApplyTimeDPG = cms.bool( False ),
    ApplyPulseDPG = cms.bool( False ),
    ECAL_Compensate = cms.bool( False ),
    ECAL_Threshold = cms.double( 10.0 ),
    ECAL_Compensation = cms.double( 0.5 ),
    ECAL_Dead_Code = cms.uint32( 10 ),
    EM_Depth = cms.double( 22.0 ),
    HAD_Depth = cms.double( 47.0 ),
    verbose = cms.untracked.bool( False ),
    thresh_Barrel = cms.double( 0.4 ),
    thresh_Endcap = cms.double( 0.4 )
)
hltParticleFlowRecHitPS = cms.EDProducer( "PFRecHitProducerPS",
    ecalRecHitsES = cms.InputTag( 'hltESRecHitAll','EcalRecHitsES' ),
    verbose = cms.untracked.bool( False ),
    thresh_Barrel = cms.double( 7.0E-6 ),
    thresh_Endcap = cms.double( 7.0E-6 )
)
hltParticleFlowClusterECALForTaus = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 0.08 ),
    thresh_Seed_Barrel = cms.double( 0.23 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 4.0 ),
    thresh_Endcap = cms.double( 0.4 ),
    thresh_Seed_Endcap = cms.double( 0.6 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 100000.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 4 ),
    posCalcNCrystal = cms.int32( 5 ),
    showerSigma = cms.double( 1.5 ),
    useCornerCells = cms.bool( True ),
    cleanRBXandHPDs = cms.bool( False ),
    depthCor_Mode = cms.int32( 1 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( "hltParticleFlowRecHitECAL" ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.04, -0.024 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.04, -0.025 )
)
hltParticleFlowClusterHCAL = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 0.8 ),
    thresh_Seed_Barrel = cms.double( 0.8 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 1.0E7 ),
    thresh_Endcap = cms.double( 0.8 ),
    thresh_Seed_Endcap = cms.double( 1.1 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 1000000.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 4 ),
    posCalcNCrystal = cms.int32( 5 ),
    showerSigma = cms.double( 10.0 ),
    useCornerCells = cms.bool( True ),
    cleanRBXandHPDs = cms.bool( True ),
    depthCor_Mode = cms.int32( 2 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( "hltParticleFlowRecHitHCAL" ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.032, -0.045 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.032, -0.045 )
)
hltParticleFlowClusterHFEM = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 0.8 ),
    thresh_Seed_Barrel = cms.double( 1.4 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 80.0 ),
    thresh_Endcap = cms.double( 0.8 ),
    thresh_Seed_Endcap = cms.double( 1.4 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 80.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 0 ),
    posCalcNCrystal = cms.int32( 5 ),
    showerSigma = cms.double( 10.0 ),
    useCornerCells = cms.bool( False ),
    cleanRBXandHPDs = cms.bool( False ),
    depthCor_Mode = cms.int32( 1 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( 'hltParticleFlowRecHitHCAL','HFEM' ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.11, -0.19 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.11, -0.19 )
)
hltParticleFlowClusterHFHAD = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 0.8 ),
    thresh_Seed_Barrel = cms.double( 1.4 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 120.0 ),
    thresh_Endcap = cms.double( 0.8 ),
    thresh_Seed_Endcap = cms.double( 1.4 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 120.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 0 ),
    posCalcNCrystal = cms.int32( 5 ),
    showerSigma = cms.double( 10.0 ),
    useCornerCells = cms.bool( False ),
    cleanRBXandHPDs = cms.bool( False ),
    depthCor_Mode = cms.int32( 2 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( 'hltParticleFlowRecHitHCAL','HFHAD' ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.045, -0.08 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.045, -0.08 )
)
hltParticleFlowClusterPS = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 6.0E-5 ),
    thresh_Seed_Barrel = cms.double( 1.2E-4 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 100000.0 ),
    thresh_Endcap = cms.double( 6.0E-5 ),
    thresh_Seed_Endcap = cms.double( 1.2E-4 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 100000.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 4 ),
    posCalcNCrystal = cms.int32( -1 ),
    showerSigma = cms.double( 0.2 ),
    useCornerCells = cms.bool( False ),
    cleanRBXandHPDs = cms.bool( False ),
    depthCor_Mode = cms.int32( 1 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( "hltParticleFlowRecHitPS" ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.0, 0.0 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.0, 0.0 )
)
hltLightPFTracksForTaus = cms.EDProducer( "LightPFTrackProducer",
    UseQuality = cms.bool( False ),
    TrackQuality = cms.string( "none" ),
    TkColList = cms.VInputTag( 'hltPFJetCtfWithMaterialTracks' )
)
hltParticleFlowBlockForTaus = cms.EDProducer( "PFBlockProducer",
    RecTracks = cms.InputTag( "hltLightPFTracksForTaus" ),
    GsfRecTracks = cms.InputTag( "pfTrackElec" ),
    ConvBremGsfRecTracks = cms.InputTag( 'pfTrackElec','Secondary' ),
    RecMuons = cms.InputTag( "muons" ),
    PFNuclear = cms.InputTag( "pfDisplacedTrackerVertex" ),
    PFConversions = cms.InputTag( "pfConversions" ),
    PFV0 = cms.InputTag( "pfV0" ),
    PFClustersECAL = cms.InputTag( "hltParticleFlowClusterECALForTaus" ),
    PFClustersHCAL = cms.InputTag( "hltParticleFlowClusterHCAL" ),
    PFClustersHFEM = cms.InputTag( "hltParticleFlowClusterHFEM" ),
    PFClustersHFHAD = cms.InputTag( "hltParticleFlowClusterHFHAD" ),
    PFClustersPS = cms.InputTag( "hltParticleFlowClusterPS" ),
    useEGPhotons = cms.bool( False ),
    EGPhotons = cms.InputTag( "" ),
    usePFatHLT = cms.bool( True ),
    useNuclear = cms.bool( False ),
    useConversions = cms.bool( False ),
    useConvBremGsfTracks = cms.bool( False ),
    useConvBremPFRecTracks = cms.bool( False ),
    useV0 = cms.bool( False ),
    useIterTracking = cms.bool( False ),
    nuclearInteractionsPurity = cms.uint32( 1 ),
    pf_DPtoverPt_Cut = cms.vdouble( -1.0, -1.0, -1.0, -1.0 ),
    pf_NHit_Cut = cms.vuint32( 0, 0, 0, 0 ),
    PhotonSelectionCuts = cms.vdouble(  ),
    useRecMuons = cms.bool( False ),
    useGsfRecTracks = cms.bool( False )
)
hltParticleFlowForTaus = cms.EDProducer( "PFProducer",
    calibHF_use = cms.bool( False ),
    blocks = cms.InputTag( "hltParticleFlowBlockForTaus" ),
    muons = cms.InputTag( "" ),
    postMuonCleaning = cms.bool( False ),
    usePFElectrons = cms.bool( False ),
    usePFPhotons = cms.bool( False ),
    useEGammaElectrons = cms.bool( False ),
    egammaElectrons = cms.InputTag( "" ),
    pf_electron_output_col = cms.string( "electrons" ),
    usePFSCEleCalib = cms.bool( True ),
    useEGammaSupercluster = cms.bool( False ),
    sumEtEcalIsoForEgammaSC_barrel = cms.double( 1.0 ),
    sumEtEcalIsoForEgammaSC_endcap = cms.double( 2.0 ),
    coneEcalIsoForEgammaSC = cms.double( 0.3 ),
    sumPtTrackIsoForEgammaSC_barrel = cms.double( 4.0 ),
    sumPtTrackIsoForEgammaSC_endcap = cms.double( 4.0 ),
    coneTrackIsoForEgammaSC = cms.double( 0.3 ),
    nTrackIsoForEgammaSC = cms.uint32( 2 ),
    pf_nsigma_ECAL = cms.double( 0.0 ),
    pf_nsigma_HCAL = cms.double( 1.0 ),
    pf_electron_mvaCut = cms.double( -0.1 ),
    pf_electronID_mvaWeightFile = cms.string( "RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_PfElectrons23Jan_IntToFloat.txt" ),
    pf_electronID_crackCorrection = cms.bool( False ),
    pf_convID_mvaWeightFile = cms.string( "RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_pfConversionFeb2311.txt" ),
    pf_conv_mvaCut = cms.double( 0.0 ),
    rejectTracks_Bad = cms.bool( False ),
    rejectTracks_Step45 = cms.bool( False ),
    usePFNuclearInteractions = cms.bool( False ),
    usePFConversions = cms.bool( False ),
    usePFDecays = cms.bool( False ),
    dptRel_DispVtx = cms.double( 10.0 ),
    useCalibrationsFromDB = cms.bool( True ),
    algoType = cms.uint32( 0 ),
    nsigma_TRACK = cms.double( 1.0 ),
    pt_Error = cms.double( 1.0 ),
    usePFMuonMomAssign = cms.bool( False ),
    postHFCleaning = cms.bool( False ),
    minHFCleaningPt = cms.double( 5.0 ),
    minSignificance = cms.double( 2.5 ),
    maxSignificance = cms.double( 2.5 ),
    minSignificanceReduction = cms.double( 1.4 ),
    maxDeltaPhiPt = cms.double( 7.0 ),
    minDeltaMet = cms.double( 0.4 ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    useVerticesForNeutral = cms.bool( True ),
    calibHF_eta_step = cms.vdouble( 0.0, 2.9, 3.0, 3.2, 4.2, 4.4, 4.6, 4.8, 5.2, 5.2 ),
    calibHF_a_EMonly = cms.vdouble( 0.96945, 0.96701, 0.76309, 0.82268, 0.87583 ),
    calibHF_b_HADonly = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348 ),
    calibHF_a_EMHAD = cms.vdouble( 1.42215, 1.00496, 0.68961, 0.81656, 0.98504 ),
    calibHF_b_EMHAD = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348 ),
    calibPFSCEle_Fbrem_barrel = cms.vdouble( 0.6, 6.0, -0.0255975, 0.0576727, 0.975442, -5.46394E-4, 1.26147, 25.0, -0.02025, 0.04537, 0.9728, -8.962E-4, 1.172 ),
    calibPFSCEle_Fbrem_endcap = cms.vdouble( 0.9, 6.5, -0.0692932, 0.101776, 0.995338, -0.00236548, 0.874998, 1.653, -0.0750184, 0.147, 0.923165, 4.74665E-4, 1.10782 ),
    calibPFSCEle_barrel = cms.vdouble( 1.004, -1.536, 22.88, -1.467, 0.3555, 0.6227, 14.65, 2051.0, 25.0, 0.9932, -0.5444, 0.0, 0.5438, 0.7109, 7.645, 0.2904, 0.0 ),
    calibPFSCEle_endcap = cms.vdouble( 1.153, -16.5975, 5.668, -0.1772, 16.22, 7.326, 0.0483, -4.068, 9.406 ),
    muon_HCAL = cms.vdouble( 3.0, 3.0 ),
    muon_ECAL = cms.vdouble( 0.5, 0.5 ),
    factors_45 = cms.vdouble( 10.0, 100.0 ),
    cleanedHF = cms.VInputTag( 'hltParticleFlowRecHitHCAL:Cleaned','hltParticleFlowClusterHFHAD:Cleaned','hltParticleFlowClusterHFEM:Cleaned' ),
    iCfgCandConnector = cms.PSet( 
      bCalibSecondary = cms.bool( False ),
      bCalibPrimary = cms.bool( False ),
      bCorrect = cms.bool( False ),
      nuclCalibFactors = cms.vdouble( 0.88, 0.28, 0.04 )
    ),
    pf_clusterRecovery = cms.bool( False ),
    ecalHcalEcalEndcap = cms.vdouble( 0.46, 3.0, 1.1, 0.4, -0.02, 1.4 ),
    ecalHcalEcalBarrel = cms.vdouble( 0.67, 3.0, 1.15, 0.9, -0.06, 1.4 ),
    ecalHcalHcalBarrel = cms.vdouble( 0.46, 3.0, 1.15, 0.3, -0.02, 1.4 ),
    ecalHcalHcalEndcap = cms.vdouble( 0.46, 3.0, 1.1, 0.3, -0.02, 1.4 )
)
hltAntiKT5PFJetsForTaus = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 0 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.5 ),
    src = cms.InputTag( "hltParticleFlowForTaus" ),
    srcPVs = cms.InputTag( "hltPixelVertices" ),
    jetType = cms.string( "PFJet" ),
    jetPtMin = cms.double( 10.0 ),
    inputEtMin = cms.double( 0.0 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltAntiKT5ConvPFJetsForTaus = cms.EDProducer( "PFJetToCaloProducer",
    Source = cms.InputTag( "hltAntiKT5PFJetsForTaus" )
)
hltPFTauJetTracksAssociator = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltAntiKT5PFJetsForTaus" ),
    tracks = cms.InputTag( "hltPFJetCtfWithMaterialTracks" ),
    coneSize = cms.double( 0.5 )
)
hltPFTauTagInfo = cms.EDProducer( "PFRecoTauTagInfoProducer",
    PFCandidateProducer = cms.InputTag( "hltParticleFlowForTaus" ),
    PFJetTracksAssociatorProducer = cms.InputTag( "hltPFTauJetTracksAssociator" ),
    PVProducer = cms.InputTag( "hltPixelVertices" ),
    smearedPVsigmaX = cms.double( 0.0015 ),
    smearedPVsigmaY = cms.double( 0.0015 ),
    smearedPVsigmaZ = cms.double( 0.0050 ),
    ChargedHadrCand_AssociationCone = cms.double( 0.8 ),
    ChargedHadrCand_tkminPt = cms.double( 0.0 ),
    ChargedHadrCand_tkminPixelHitsn = cms.int32( 0 ),
    ChargedHadrCand_tkminTrackerHitsn = cms.int32( 0 ),
    ChargedHadrCand_tkmaxipt = cms.double( 0.2 ),
    ChargedHadrCand_tkmaxChi2 = cms.double( 100.0 ),
    NeutrHadrCand_HcalclusMinEt = cms.double( 0.5 ),
    GammaCand_EcalclusMinEt = cms.double( 0.5 ),
    tkminPt = cms.double( 0.0 ),
    tkminPixelHitsn = cms.int32( 0 ),
    tkminTrackerHitsn = cms.int32( 3 ),
    tkmaxipt = cms.double( 0.2 ),
    tkmaxChi2 = cms.double( 100.0 ),
    UsePVconstraint = cms.bool( True ),
    ChargedHadrCand_tkPVmaxDZ = cms.double( 0.4 ),
    tkPVmaxDZ = cms.double( 0.4 )
)
hltPFTausTightIso = cms.EDProducer( "PFRecoTauProducer",
    PFTauTagInfoProducer = cms.InputTag( "hltPFTauTagInfo" ),
    ElectronPreIDProducer = cms.InputTag( "elecpreid" ),
    PVProducer = cms.InputTag( "hltPixelVertices" ),
    Algorithm = cms.string( "ConeBased" ),
    smearedPVsigmaX = cms.double( 0.0015 ),
    smearedPVsigmaY = cms.double( 0.0015 ),
    smearedPVsigmaZ = cms.double( 0.0050 ),
    JetPtMin = cms.double( 0.0 ),
    LeadPFCand_minPt = cms.double( 0.0 ),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool( True ),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double( 0.4 ),
    LeadTrack_minPt = cms.double( 0.0 ),
    UseTrackLeadTrackDZconstraint = cms.bool( True ),
    TrackLeadTrack_maxDZ = cms.double( 0.4 ),
    MatchingConeMetric = cms.string( "DR" ),
    MatchingConeSizeFormula = cms.string( "0.2" ),
    MatchingConeSize_min = cms.double( 0.0 ),
    MatchingConeSize_max = cms.double( 0.6 ),
    TrackerSignalConeMetric = cms.string( "DR" ),
    TrackerSignalConeSizeFormula = cms.string( "0.15" ),
    TrackerSignalConeSize_min = cms.double( 0.0 ),
    TrackerSignalConeSize_max = cms.double( 0.2 ),
    TrackerIsolConeMetric = cms.string( "DR" ),
    TrackerIsolConeSizeFormula = cms.string( "0.5" ),
    TrackerIsolConeSize_min = cms.double( 0.0 ),
    TrackerIsolConeSize_max = cms.double( 0.5 ),
    ECALSignalConeMetric = cms.string( "DR" ),
    ECALSignalConeSizeFormula = cms.string( "0.15" ),
    ECALSignalConeSize_min = cms.double( 0.0 ),
    ECALSignalConeSize_max = cms.double( 0.6 ),
    ECALIsolConeMetric = cms.string( "DR" ),
    ECALIsolConeSizeFormula = cms.string( "0.5" ),
    ECALIsolConeSize_min = cms.double( 0.0 ),
    ECALIsolConeSize_max = cms.double( 0.5 ),
    HCALSignalConeMetric = cms.string( "DR" ),
    HCALSignalConeSizeFormula = cms.string( "0.2" ),
    HCALSignalConeSize_min = cms.double( 0.0 ),
    HCALSignalConeSize_max = cms.double( 0.5 ),
    HCALIsolConeMetric = cms.string( "DR" ),
    HCALIsolConeSizeFormula = cms.string( "0.5" ),
    HCALIsolConeSize_min = cms.double( 0.0 ),
    HCALIsolConeSize_max = cms.double( 0.5 ),
    Rphi = cms.double( 0.2 ),
    MaxEtInEllipse = cms.double( 2.0 ),
    AddEllipseGammas = cms.bool( False ),
    AreaMetric_recoElements_maxabsEta = cms.double( 2.5 ),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32( 0 ),
    Track_IsolAnnulus_minNhits = cms.uint32( 0 ),
    ElecPreIDLeadTkMatch_maxDR = cms.double( 0.015 ),
    EcalStripSumE_minClusEnergy = cms.double( 0.0 ),
    EcalStripSumE_deltaEta = cms.double( 0.0 ),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double( 0.0 ),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double( 0.0 ),
    maximumForElectrionPreIDOutput = cms.double( 0.0 ),
    DataType = cms.string( "AOD" ),
    emMergingAlgorithm = cms.string( "None" ),
    candOverlapCriterion = cms.string( "None" ),
    doOneProng = cms.bool( True ),
    doOneProngStrip = cms.bool( True ),
    doOneProngTwoStrips = cms.bool( True ),
    doThreeProng = cms.bool( True ),
    tauPtThreshold = cms.double( 0.0 ),
    leadPionThreshold = cms.double( 1.0 ),
    stripPtThreshold = cms.double( 0.5 ),
    chargeHadrIsolationConeSize = cms.double( 0.5 ),
    gammaIsolationConeSize = cms.double( 0.5 ),
    neutrHadrIsolationConeSize = cms.double( 0.5 ),
    useIsolationAnnulus = cms.bool( False ),
    matchingCone = cms.double( 0.2 ),
    coneMetric = cms.string( "DR" ),
    coneSizeFormula = cms.string( "2.8/ET" ),
    minimumSignalCone = cms.double( 0.0 ),
    maximumSignalCone = cms.double( 1.8 ),
    oneProngStripMassWindow = cms.vdouble( 0.0, 0.0 ),
    oneProngTwoStripsMassWindow = cms.vdouble( 0.0, 0.0 ),
    oneProngTwoStripsPi0MassWindow = cms.vdouble( 0.0, 0.0 ),
    threeProngMassWindow = cms.vdouble( 0.0, 0.0 )
)
hltPFTauTightIsoTrackFindingDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByLeadingObjectPtCut",
    PFTauProducer = cms.InputTag( "hltPFTausTightIso" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
    UseOnlyChargedHadrons = cms.bool( True ),
    MinPtLeadingObject = cms.double( 0.0 )
)
hltPFTauTightIsoTrackPt5Discriminator = cms.EDProducer( "PFRecoTauDiscriminationByLeadingObjectPtCut",
    PFTauProducer = cms.InputTag( "hltPFTausTightIso" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
    UseOnlyChargedHadrons = cms.bool( True ),
    MinPtLeadingObject = cms.double( 5.0 )
)
hltPFTauTightIsoIsolationDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByTrackerIsolation = cms.bool( True ),
    ApplyDiscriminationByECALIsolation = cms.bool( True ),
    applyOccupancyCut = cms.bool( True ),
    maximumOccupancy = cms.uint32( 0 ),
    applySumPtCut = cms.bool( False ),
    maximumSumPtCut = cms.double( 6.0 ),
    applyRelativeSumPtCut = cms.bool( False ),
    relativeSumPtCut = cms.double( 0.0 ),
    PVProducer = cms.InputTag( "hltPixelVertices" ),
    customOuterCone = cms.double( -1.0 ),
    qualityCuts = cms.PSet( 
      isolationQualityCuts = cms.PSet( 
        minTrackHits = cms.uint32( 3 ),
        minTrackPt = cms.double( 1.0 ),
        maxTrackChi2 = cms.double( 100.0 ),
        minTrackPixelHits = cms.uint32( 0 ),
        minGammaEt = cms.double( 1.5 ),
        useTracksInsteadOfPFHadrons = cms.bool( False ),
        maxDeltaZ = cms.double( 0.2 ),
        maxTransverseImpactParameter = cms.double( 0.05 )
      ),
      signalQualityCuts = cms.PSet( 
        maxDeltaZ = cms.double( 0.5 ),
        minTrackPt = cms.double( 0.0 ),
        maxTrackChi2 = cms.double( 1000.0 ),
        useTracksInsteadOfPFHadrons = cms.bool( False ),
        minGammaEt = cms.double( 0.5 ),
        minTrackPixelHits = cms.uint32( 0 ),
        minTrackHits = cms.uint32( 3 ),
        maxTransverseImpactParameter = cms.double( 0.2 )
      )
    ),
    PFTauProducer = cms.InputTag( "hltPFTausTightIso" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) )
)
hltSelectedPFTausTightIsoTrackFinding = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackFindingDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltSelectedPFTausTightIsoTrackPt5 = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackPt5Discriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltSelectedPFTausTightIsoTrackFindingIsolation = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackFindingDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      ),
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoIsolationDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltSelectedPFTausTightIsoTrackPt5Isolation = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackPt5Discriminator" ),
        selectionCut = cms.double( 0.5 )
      ),
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoIsolationDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltConvPFTausTightIsoTrackFinding = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTightIsoTrackFinding" )
)
hltConvPFTausTightIsoTrackFindingIsolation = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTightIsoTrackFindingIsolation" )
)
hltConvPFTausTightIso = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltPFTausTightIso" )
)
hltConvPFTausTightIsoTrackPt5 = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTightIsoTrackPt5" )
)
hltConvPFTausTightIsoTrackPt5Isolation = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTightIsoTrackPt5Isolation" )
)
hltPFTau5Track = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIsoTrackFinding" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 5.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau5Track5 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIsoTrackPt5" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 5.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltFilterPFTauTrack5TightIsoL1QuadJet20Central = cms.EDProducer( "L1HLTJetsMatching",
    JetSrc = cms.InputTag( "hltConvPFTausTightIsoTrackPt5Isolation" ),
    L1TauTrigger = cms.InputTag( "hltL1sL1QuadJet20Central" ),
    EtMin = cms.double( 0.0 )
)
hltFilterPFTauTrack5TightIsoL1QuadJet20CentralPFTau40 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltFilterPFTauTrack5TightIsoL1QuadJet20Central" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPreQuadJet50 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet50Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 50.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 4 )
)
hltSelector4Jets = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBLifetimeL25Jets = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelector4Jets" ),
    filter = cms.bool( False ),
    etMin = cms.double( 25.0 )
)
hltBLifetimeL25Associator = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25Jets" ),
    tracks = cms.InputTag( "hltPixelTracks" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL25TagInfos = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL25Associator" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL25BJetTags = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL25TagInfos' )
)
hltBLifetimeL25Filter = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL25BJetTags" ),
    MinTag = cms.double( 0.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBLifetimeL3Associator = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25Jets" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracks" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfos = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3Associator" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTags = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfos' )
)
hltBLifetimeL3Filter = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTags" ),
    MinTag = cms.double( 2.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltPreQuadJet50Jet40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPentaJet40Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 5 )
)
hltPreQuadJet60 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet60 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 60.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 4 )
)
hltPreQuadJet70 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet70 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 70.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 4 )
)
hltPreExclDiJet60HFOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltExclDiJet60HFOR = cms.EDFilter( "HLTExclDiJetFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minPtJet = cms.double( 60.0 ),
    minHFe = cms.double( 70.0 ),
    HF_OR = cms.bool( True )
)
hltL1sL1SingleJet36FwdVeto = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet36_FwdVeto" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreExclDiJet60HFAND = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltExclDiJet60HFAND = cms.EDFilter( "HLTExclDiJetFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minPtJet = cms.double( 60.0 ),
    minHFe = cms.double( 70.0 ),
    HF_OR = cms.bool( False )
)
hltL1sL1HTT50 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HTT50" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreHT150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT150 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 150.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltL1sL1HTT75 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HTT75" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreHLTHT150AlphaT0p6 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT150AlphaT0p6 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 150.0 ),
    minAlphaT = cms.double( 0.6 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT200 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 200.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltL1sL1HTT100 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HTT100" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreHLTHT200AlphaT0p53 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT200AlphaT0p53 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 200.0 ),
    minAlphaT = cms.double( 0.53 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHLTHT200AlphaT0p6 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT200AlphaT0p6 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 200.0 ),
    minAlphaT = cms.double( 0.6 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT250 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT250 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 250.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHLTHT250AlphaT0p53 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT250AlphaT0p53 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 250.0 ),
    minAlphaT = cms.double( 0.53 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHLTHT250AlphaT0p54 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT250AlphaT0p54 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 250.0 ),
    minAlphaT = cms.double( 0.54 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT250MHT60 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMHT60 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 60.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 30.0, 30.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT250MHT70 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMHT70 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 70.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 30.0, 30.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT300 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 300.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHT300MHT75 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMHT75 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 75.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 30.0, 30.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT300PFMHT55 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltParticleFlowClusterECAL = cms.EDProducer( "PFClusterProducer",
    thresh_Barrel = cms.double( 0.8 ),
    thresh_Seed_Barrel = cms.double( 0.8 ),
    thresh_Pt_Barrel = cms.double( 0.0 ),
    thresh_Pt_Seed_Barrel = cms.double( 0.0 ),
    thresh_Clean_Barrel = cms.double( 100000.0 ),
    thresh_Endcap = cms.double( 0.8 ),
    thresh_Seed_Endcap = cms.double( 1.1 ),
    thresh_Pt_Endcap = cms.double( 0.0 ),
    thresh_Pt_Seed_Endcap = cms.double( 0.0 ),
    thresh_Clean_Endcap = cms.double( 100000.0 ),
    thresh_DoubleSpike_Barrel = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Barrel = cms.double( -1.0 ),
    thresh_DoubleSpike_Endcap = cms.double( 1.0E9 ),
    minS6S2_DoubleSpike_Endcap = cms.double( -1.0 ),
    nNeighbours = cms.int32( 4 ),
    posCalcNCrystal = cms.int32( 5 ),
    showerSigma = cms.double( 10.0 ),
    useCornerCells = cms.bool( True ),
    cleanRBXandHPDs = cms.bool( False ),
    depthCor_Mode = cms.int32( 1 ),
    depthCor_A = cms.double( 0.89 ),
    depthCor_B = cms.double( 7.4 ),
    depthCor_A_preshower = cms.double( 0.89 ),
    depthCor_B_preshower = cms.double( 4.0 ),
    PFRecHits = cms.InputTag( "hltParticleFlowRecHitECAL" ),
    minS4S1_Clean_Barrel = cms.vdouble( 0.04, -0.024 ),
    minS4S1_Clean_Endcap = cms.vdouble( 0.04, -0.025 )
)
hltLightPFTracks = cms.EDProducer( "LightPFTrackProducer",
    UseQuality = cms.bool( False ),
    TrackQuality = cms.string( "none" ),
    TkColList = cms.VInputTag( 'hltPFlowTrackSelectionHighPurity' )
)
hltParticleFlowBlock = cms.EDProducer( "PFBlockProducer",
    RecTracks = cms.InputTag( "hltLightPFTracks" ),
    GsfRecTracks = cms.InputTag( "pfTrackElec" ),
    ConvBremGsfRecTracks = cms.InputTag( 'pfTrackElec','Secondary' ),
    RecMuons = cms.InputTag( "muons" ),
    PFNuclear = cms.InputTag( "pfDisplacedTrackerVertex" ),
    PFConversions = cms.InputTag( "pfConversions" ),
    PFV0 = cms.InputTag( "pfV0" ),
    PFClustersECAL = cms.InputTag( "hltParticleFlowClusterECAL" ),
    PFClustersHCAL = cms.InputTag( "hltParticleFlowClusterHCAL" ),
    PFClustersHFEM = cms.InputTag( "hltParticleFlowClusterHFEM" ),
    PFClustersHFHAD = cms.InputTag( "hltParticleFlowClusterHFHAD" ),
    PFClustersPS = cms.InputTag( "hltParticleFlowClusterPS" ),
    useEGPhotons = cms.bool( False ),
    EGPhotons = cms.InputTag( "" ),
    usePFatHLT = cms.bool( True ),
    useNuclear = cms.bool( False ),
    useConversions = cms.bool( False ),
    useConvBremGsfTracks = cms.bool( False ),
    useConvBremPFRecTracks = cms.bool( False ),
    useV0 = cms.bool( False ),
    useIterTracking = cms.bool( False ),
    nuclearInteractionsPurity = cms.uint32( 1 ),
    pf_DPtoverPt_Cut = cms.vdouble( -1.0, -1.0, -1.0, -1.0 ),
    pf_NHit_Cut = cms.vuint32( 0, 0, 0, 0 ),
    PhotonSelectionCuts = cms.vdouble(  ),
    useRecMuons = cms.bool( False ),
    useGsfRecTracks = cms.bool( False )
)
hltParticleFlow = cms.EDProducer( "PFProducer",
    calibHF_use = cms.bool( False ),
    blocks = cms.InputTag( "hltParticleFlowBlock" ),
    muons = cms.InputTag( "" ),
    postMuonCleaning = cms.bool( False ),
    usePFElectrons = cms.bool( False ),
    usePFPhotons = cms.bool( False ),
    useEGammaElectrons = cms.bool( False ),
    egammaElectrons = cms.InputTag( "" ),
    pf_electron_output_col = cms.string( "electrons" ),
    usePFSCEleCalib = cms.bool( True ),
    useEGammaSupercluster = cms.bool( False ),
    sumEtEcalIsoForEgammaSC_barrel = cms.double( 1.0 ),
    sumEtEcalIsoForEgammaSC_endcap = cms.double( 2.0 ),
    coneEcalIsoForEgammaSC = cms.double( 0.3 ),
    sumPtTrackIsoForEgammaSC_barrel = cms.double( 4.0 ),
    sumPtTrackIsoForEgammaSC_endcap = cms.double( 4.0 ),
    coneTrackIsoForEgammaSC = cms.double( 0.3 ),
    nTrackIsoForEgammaSC = cms.uint32( 2 ),
    pf_nsigma_ECAL = cms.double( 0.0 ),
    pf_nsigma_HCAL = cms.double( 1.0 ),
    pf_electron_mvaCut = cms.double( -0.1 ),
    pf_electronID_mvaWeightFile = cms.string( "RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_PfElectrons23Jan_IntToFloat.txt" ),
    pf_electronID_crackCorrection = cms.bool( False ),
    pf_convID_mvaWeightFile = cms.string( "RecoParticleFlow/PFProducer/data/MVAnalysis_BDT.weights_pfConversionFeb2311.txt" ),
    pf_conv_mvaCut = cms.double( 0.0 ),
    rejectTracks_Bad = cms.bool( False ),
    rejectTracks_Step45 = cms.bool( False ),
    usePFNuclearInteractions = cms.bool( False ),
    usePFConversions = cms.bool( False ),
    usePFDecays = cms.bool( False ),
    dptRel_DispVtx = cms.double( 10.0 ),
    useCalibrationsFromDB = cms.bool( True ),
    algoType = cms.uint32( 0 ),
    nsigma_TRACK = cms.double( 1.0 ),
    pt_Error = cms.double( 1.0 ),
    usePFMuonMomAssign = cms.bool( False ),
    postHFCleaning = cms.bool( False ),
    minHFCleaningPt = cms.double( 5.0 ),
    minSignificance = cms.double( 2.5 ),
    maxSignificance = cms.double( 2.5 ),
    minSignificanceReduction = cms.double( 1.4 ),
    maxDeltaPhiPt = cms.double( 7.0 ),
    minDeltaMet = cms.double( 0.4 ),
    vertexCollection = cms.InputTag( "hltPixelVertices" ),
    useVerticesForNeutral = cms.bool( True ),
    calibHF_eta_step = cms.vdouble( 0.0, 2.9, 3.0, 3.2, 4.2, 4.4, 4.6, 4.8, 5.2, 5.2 ),
    calibHF_a_EMonly = cms.vdouble( 0.96945, 0.96701, 0.76309, 0.82268, 0.87583 ),
    calibHF_b_HADonly = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348 ),
    calibHF_a_EMHAD = cms.vdouble( 1.42215, 1.00496, 0.68961, 0.81656, 0.98504 ),
    calibHF_b_EMHAD = cms.vdouble( 1.27541, 0.85361, 0.86333, 0.89091, 0.94348 ),
    calibPFSCEle_Fbrem_barrel = cms.vdouble( 0.6, 6.0, -0.0255975, 0.0576727, 0.975442, -5.46394E-4, 1.26147, 25.0, -0.02025, 0.04537, 0.9728, -8.962E-4, 1.172 ),
    calibPFSCEle_Fbrem_endcap = cms.vdouble( 0.9, 6.5, -0.0692932, 0.101776, 0.995338, -0.00236548, 0.874998, 1.653, -0.0750184, 0.147, 0.923165, 4.74665E-4, 1.10782 ),
    calibPFSCEle_barrel = cms.vdouble( 1.004, -1.536, 22.88, -1.467, 0.3555, 0.6227, 14.65, 2051.0, 25.0, 0.9932, -0.5444, 0.0, 0.5438, 0.7109, 7.645, 0.2904, 0.0 ),
    calibPFSCEle_endcap = cms.vdouble( 1.153, -16.5975, 5.668, -0.1772, 16.22, 7.326, 0.0483, -4.068, 9.406 ),
    muon_HCAL = cms.vdouble( 3.0, 3.0 ),
    muon_ECAL = cms.vdouble( 0.5, 0.5 ),
    factors_45 = cms.vdouble( 10.0, 100.0 ),
    cleanedHF = cms.VInputTag( 'hltParticleFlowRecHitHCAL:Cleaned','hltParticleFlowClusterHFHAD:Cleaned','hltParticleFlowClusterHFEM:Cleaned' ),
    iCfgCandConnector = cms.PSet( 
      bCalibSecondary = cms.bool( False ),
      bCalibPrimary = cms.bool( False ),
      bCorrect = cms.bool( False ),
      nuclCalibFactors = cms.vdouble( 0.88, 0.28, 0.04 )
    ),
    pf_clusterRecovery = cms.bool( False ),
    ecalHcalEcalEndcap = cms.vdouble( 0.46, 3.0, 1.1, 0.4, -0.02, 1.4 ),
    ecalHcalEcalBarrel = cms.vdouble( 0.67, 3.0, 1.15, 0.9, -0.06, 1.4 ),
    ecalHcalHcalBarrel = cms.vdouble( 0.46, 3.0, 1.15, 0.3, -0.02, 1.4 ),
    ecalHcalHcalEndcap = cms.vdouble( 0.46, 3.0, 1.1, 0.3, -0.02, 1.4 )
)
hltAntiKT5PFJets = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 0 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "AntiKt" ),
    rParam = cms.double( 0.5 ),
    src = cms.InputTag( "hltParticleFlow" ),
    srcPVs = cms.InputTag( "hltPixelVertices" ),
    jetType = cms.string( "PFJet" ),
    jetPtMin = cms.double( 15.0 ),
    inputEtMin = cms.double( 0.0 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltAntiKT5ConvPFJets = cms.EDProducer( "PFJetToCaloProducer",
    Source = cms.InputTag( "hltAntiKT5PFJets" )
)
hltPFMHT55Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 55.0 ),
    minNJet = cms.int32( 3 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreHT300BTagIP = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBJetRA2b = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
hltGetJetsfromBJetRA2b = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBJetRA2b" )
)
hltSelectorJetsRA2b = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfromBJetRA2b" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 6 )
)
hltBLifetimeL25JetsRA2b = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelectorJetsRA2b" ),
    filter = cms.bool( False ),
    etMin = cms.double( 30.0 )
)
hltBLifetimeL3AssociatorRA2b = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsRA2b" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracksRA2b" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfosRA2b = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3AssociatorRA2b" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTagsRA2b = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfosRA2b' )
)
hltBLifetimeL3FilterRA2b = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTagsRA2b" ),
    MinTag = cms.double( 4.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltPreHT300BTagIPPFMHT55 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreHT300BTagIPPFMHT75 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPFMHT75Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 75.0 ),
    minNJet = cms.int32( 3 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreHLTHT300AlphaT0p52 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT300AlphaT0p52 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 300.0 ),
    minAlphaT = cms.double( 0.52 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHLTHT300AlphaT0p53 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT300AlphaT0p53 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 300.0 ),
    minAlphaT = cms.double( 0.53 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT350 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT350 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 350.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHLTHT350AlphaT0p51 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT350AlphaT0p51 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 350.0 ),
    minAlphaT = cms.double( 0.51 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHLTHT350AlphaT0p53 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT350AlphaT0p53 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 350.0 ),
    minAlphaT = cms.double( 0.53 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT400 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT400 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 400.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHLTHT400AlphaT0p51 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT400AlphaT0p51 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 5 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 400.0 ),
    minAlphaT = cms.double( 0.51 ),
    minPtJet = cms.vdouble( 40.0, 40.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPreHT450 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT450 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 450.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHT500 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT500 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 500.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPreHT550 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHT550 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 0.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 4 ),
    usePt = cms.bool( False ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 550.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 40.0 ),
    etaJet = cms.vdouble( 3.0 )
)
hltPrePFMHT150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPFMHT150Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 150.0 ),
    minNJet = cms.int32( 3 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreMET100 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreMET120 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMET120 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 120.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltPreMET200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMET200 = cms.EDFilter( "HLT1CaloMET",
    inputTag = cms.InputTag( "hltMet" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 200.0 ),
    MaxEta = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1DoubleJet36Central = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleJet36_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreR014MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJet56NoID = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 56.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltRHemisphere = cms.EDFilter( "HLTRHemisphere",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minJetPt = cms.double( 40.0 ),
    maxEta = cms.double( 3.0 ),
    maxNJ = cms.int32( 7 ),
    acceptNJ = cms.bool( True )
)
hltR014MR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.14 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR014MR150BTag = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBJetRAzr = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
hltGetJetsfromBJetRAzr = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBJetRAzr" )
)
hltSelectorJetsRAzr = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfromBJetRAzr" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 6 )
)
hltBLifetimeL25JetsRAzr = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelectorJetsRAzr" ),
    filter = cms.bool( False ),
    etMin = cms.double( 30.0 )
)
hltBLifetimeL3AssociatorRAzr = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsRAzr" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracksRAzr" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfosRAzr = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3AssociatorRAzr" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTagsRAzr = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfosRAzr' )
)
hltBLifetimeL3FilterRAzr = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTagsRAzr" ),
    MinTag = cms.double( 4.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltPreR014MR450BTag = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR014MR450 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.14 ),
    minMR = cms.double( 450.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR014MR600 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR014MR600 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.14 ),
    minMR = cms.double( 600.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR020MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR020MR350BTag = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR350 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 350.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR020MR500 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR500 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 500.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR020MR550 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR550 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 550.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR025MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR025MR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.25 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR025MR250BTag = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR025MR250 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.25 ),
    minMR = cms.double( 250.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR025MR400 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR025MR400 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.25 ),
    minMR = cms.double( 400.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR025MR450 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR025MR450 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.25 ),
    minMR = cms.double( 450.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR033MR300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR033MR300 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.33 ),
    minMR = cms.double( 300.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR033MR350 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR033MR350 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.33 ),
    minMR = cms.double( 350.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR038MR200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR038MR200 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.38 ),
    minMR = cms.double( 200.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreR038MR250 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR038MR250 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.38 ),
    minMR = cms.double( 250.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltL1sL1SingleMuOpen = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMuOpen" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1SingleMuOpen = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenL1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpen" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltPreL1SingleMuOpenDT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenL1FilteredDT = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpen" ),
    MaxEta = cms.double( 1.25 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltL1sL1SingleMu10 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu10" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1Mu10 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu10L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu10" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1sL1SingleMu20 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu20" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1Mu20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu20L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu20" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1sL1DoubleMu0 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleMu0" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1DoubleMu0 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiMuonL1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu0" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltPreL2Mu10 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDt1DRecHits = cms.EDProducer( "DTRecHitProducer",
    dtDigiLabel = cms.InputTag( "simMuonDTDigis" ),
    recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
    recAlgoConfig = cms.PSet( 
      minTime = cms.double( -3.0 ),
      debug = cms.untracked.bool( False ),
      tTrigModeConfig = cms.PSet( 
        vPropWire = cms.double( 24.4 ),
        doTOFCorrection = cms.bool( True ),
        tofCorrType = cms.int32( 0 ),
        wirePropCorrType = cms.int32( 0 ),
        tTrigLabel = cms.string( "" ),
        doWirePropCorrection = cms.bool( True ),
        doT0Correction = cms.bool( True ),
        debug = cms.untracked.bool( False )
      ),
      maxTime = cms.double( 420.0 ),
      tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
      stepTwoFromDigi = cms.bool( False ),
      doVdriftCorr = cms.bool( False )
    )
)
hltDt4DSegments = cms.EDProducer( "DTRecSegment4DProducer",
    debug = cms.untracked.bool( False ),
    recHits1DLabel = cms.InputTag( "hltDt1DRecHits" ),
    recHits2DLabel = cms.InputTag( "dt2DSegments" ),
    Reco4DAlgoName = cms.string( "DTCombinatorialPatternReco4D" ),
    Reco4DAlgoConfig = cms.PSet( 
      segmCleanerMode = cms.int32( 2 ),
      Reco2DAlgoName = cms.string( "DTCombinatorialPatternReco" ),
      recAlgoConfig = cms.PSet( 
        minTime = cms.double( -3.0 ),
        debug = cms.untracked.bool( False ),
        tTrigModeConfig = cms.PSet( 
          vPropWire = cms.double( 24.4 ),
          doTOFCorrection = cms.bool( True ),
          tofCorrType = cms.int32( 0 ),
          wirePropCorrType = cms.int32( 0 ),
          tTrigLabel = cms.string( "" ),
          doWirePropCorrection = cms.bool( True ),
          doT0Correction = cms.bool( True ),
          debug = cms.untracked.bool( False )
        ),
        maxTime = cms.double( 420.0 ),
        tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
        stepTwoFromDigi = cms.bool( False ),
        doVdriftCorr = cms.bool( False )
      ),
      nSharedHitsMax = cms.int32( 2 ),
      hit_afterT0_resolution = cms.double( 0.03 ),
      Reco2DAlgoConfig = cms.PSet( 
        segmCleanerMode = cms.int32( 2 ),
        recAlgoConfig = cms.PSet( 
          minTime = cms.double( -3.0 ),
          debug = cms.untracked.bool( False ),
          tTrigModeConfig = cms.PSet( 
            vPropWire = cms.double( 24.4 ),
            doTOFCorrection = cms.bool( True ),
            tofCorrType = cms.int32( 0 ),
            wirePropCorrType = cms.int32( 0 ),
            tTrigLabel = cms.string( "" ),
            doWirePropCorrection = cms.bool( True ),
            doT0Correction = cms.bool( True ),
            debug = cms.untracked.bool( False )
          ),
          maxTime = cms.double( 420.0 ),
          tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
          stepTwoFromDigi = cms.bool( False ),
          doVdriftCorr = cms.bool( False )
        ),
        nSharedHitsMax = cms.int32( 2 ),
        AlphaMaxPhi = cms.double( 1.0 ),
        hit_afterT0_resolution = cms.double( 0.03 ),
        MaxAllowedHits = cms.uint32( 50 ),
        performT0_vdriftSegCorrection = cms.bool( False ),
        AlphaMaxTheta = cms.double( 0.9 ),
        debug = cms.untracked.bool( False ),
        recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
        nUnSharedHitsMin = cms.int32( 2 ),
        performT0SegCorrection = cms.bool( False )
      ),
      performT0_vdriftSegCorrection = cms.bool( False ),
      debug = cms.untracked.bool( False ),
      recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
      nUnSharedHitsMin = cms.int32( 2 ),
      AllDTRecHits = cms.bool( True ),
      performT0SegCorrection = cms.bool( False )
    )
)
hltCsc2DRecHits = cms.EDProducer( "CSCRecHitDProducer",
    CSCUseCalibrations = cms.bool( True ),
    CSCUseStaticPedestals = cms.bool( False ),
    CSCUseTimingCorrections = cms.bool( True ),
    stripDigiTag = cms.InputTag( 'simMuonCSCDigis','MuonCSCStripDigi' ),
    wireDigiTag = cms.InputTag( 'simMuonCSCDigis','MuonCSCWireDigi' ),
    CSCstripWireDeltaTime = cms.int32( 8 ),
    CSCNoOfTimeBinsForDynamicPedestal = cms.int32( 2 ),
    CSCStripPeakThreshold = cms.double( 10.0 ),
    CSCStripClusterChargeCut = cms.double( 25.0 ),
    CSCWireClusterDeltaT = cms.int32( 1 ),
    CSCStripxtalksOffset = cms.double( 0.03 ),
    NoiseLevel_ME1a = cms.double( 7.0 ),
    XTasymmetry_ME1a = cms.double( 0.0 ),
    ConstSyst_ME1a = cms.double( 0.022 ),
    NoiseLevel_ME1b = cms.double( 8.0 ),
    XTasymmetry_ME1b = cms.double( 0.0 ),
    ConstSyst_ME1b = cms.double( 0.0070 ),
    NoiseLevel_ME12 = cms.double( 9.0 ),
    XTasymmetry_ME12 = cms.double( 0.0 ),
    ConstSyst_ME12 = cms.double( 0.0 ),
    NoiseLevel_ME13 = cms.double( 8.0 ),
    XTasymmetry_ME13 = cms.double( 0.0 ),
    ConstSyst_ME13 = cms.double( 0.0 ),
    NoiseLevel_ME21 = cms.double( 9.0 ),
    XTasymmetry_ME21 = cms.double( 0.0 ),
    ConstSyst_ME21 = cms.double( 0.0 ),
    NoiseLevel_ME22 = cms.double( 9.0 ),
    XTasymmetry_ME22 = cms.double( 0.0 ),
    ConstSyst_ME22 = cms.double( 0.0 ),
    NoiseLevel_ME31 = cms.double( 9.0 ),
    XTasymmetry_ME31 = cms.double( 0.0 ),
    ConstSyst_ME31 = cms.double( 0.0 ),
    NoiseLevel_ME32 = cms.double( 9.0 ),
    XTasymmetry_ME32 = cms.double( 0.0 ),
    ConstSyst_ME32 = cms.double( 0.0 ),
    NoiseLevel_ME41 = cms.double( 9.0 ),
    XTasymmetry_ME41 = cms.double( 0.0 ),
    ConstSyst_ME41 = cms.double( 0.0 ),
    readBadChannels = cms.bool( True ),
    readBadChambers = cms.bool( True ),
    UseAverageTime = cms.bool( False ),
    UseParabolaFit = cms.bool( False ),
    UseFivePoleFit = cms.bool( True )
)
hltCscSegments = cms.EDProducer( "CSCSegmentProducer",
    inputObjects = cms.InputTag( "hltCsc2DRecHits" ),
    algo_type = cms.int32( 1 ),
    algo_psets = cms.VPSet( 
      cms.PSet(  chamber_types = cms.vstring( 'ME1/a',
  'ME1/b',
  'ME1/2',
  'ME1/3',
  'ME2/1',
  'ME2/2',
  'ME3/1',
  'ME3/2',
  'ME4/1',
  'ME4/2' ),
        algo_name = cms.string( "CSCSegAlgoST" ),
        parameters_per_chamber_type = cms.vint32( 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
        algo_psets = cms.VPSet( 
          cms.PSet(  maxRatioResidualPrune = cms.double( 3.0 ),
            yweightPenalty = cms.double( 1.5 ),
            maxRecHitsInCluster = cms.int32( 20 ),
            dPhiFineMax = cms.double( 0.025 ),
            preClusteringUseChaining = cms.bool( True ),
            ForceCovariance = cms.bool( False ),
            hitDropLimit6Hits = cms.double( 0.3333 ),
            NormChi2Cut2D = cms.double( 20.0 ),
            BPMinImprovement = cms.double( 10000.0 ),
            Covariance = cms.double( 0.0 ),
            tanPhiMax = cms.double( 0.5 ),
            SeedBig = cms.double( 0.0015 ),
            onlyBestSegment = cms.bool( False ),
            dRPhiFineMax = cms.double( 8.0 ),
            SeedSmall = cms.double( 2.0E-4 ),
            curvePenalty = cms.double( 2.0 ),
            dXclusBoxMax = cms.double( 4.0 ),
            BrutePruning = cms.bool( True ),
            curvePenaltyThreshold = cms.double( 0.85 ),
            CorrectTheErrors = cms.bool( True ),
            hitDropLimit4Hits = cms.double( 0.6 ),
            useShowering = cms.bool( False ),
            CSCDebug = cms.untracked.bool( False ),
            tanThetaMax = cms.double( 1.2 ),
            NormChi2Cut3D = cms.double( 10.0 ),
            minHitsPerSegment = cms.int32( 3 ),
            ForceCovarianceAll = cms.bool( False ),
            yweightPenaltyThreshold = cms.double( 1.0 ),
            prePrunLimit = cms.double( 3.17 ),
            hitDropLimit5Hits = cms.double( 0.8 ),
            preClustering = cms.bool( True ),
            prePrun = cms.bool( True ),
            maxDPhi = cms.double( 999.0 ),
            maxDTheta = cms.double( 999.0 ),
            Pruning = cms.bool( True ),
            dYclusBoxMax = cms.double( 8.0 )
          ),
          cms.PSet(  maxRatioResidualPrune = cms.double( 3.0 ),
            yweightPenalty = cms.double( 1.5 ),
            maxRecHitsInCluster = cms.int32( 24 ),
            dPhiFineMax = cms.double( 0.025 ),
            preClusteringUseChaining = cms.bool( True ),
            ForceCovariance = cms.bool( False ),
            hitDropLimit6Hits = cms.double( 0.3333 ),
            NormChi2Cut2D = cms.double( 20.0 ),
            BPMinImprovement = cms.double( 10000.0 ),
            Covariance = cms.double( 0.0 ),
            tanPhiMax = cms.double( 0.5 ),
            SeedBig = cms.double( 0.0015 ),
            onlyBestSegment = cms.bool( False ),
            dRPhiFineMax = cms.double( 8.0 ),
            SeedSmall = cms.double( 2.0E-4 ),
            curvePenalty = cms.double( 2.0 ),
            dXclusBoxMax = cms.double( 4.0 ),
            BrutePruning = cms.bool( True ),
            curvePenaltyThreshold = cms.double( 0.85 ),
            CorrectTheErrors = cms.bool( True ),
            hitDropLimit4Hits = cms.double( 0.6 ),
            useShowering = cms.bool( False ),
            CSCDebug = cms.untracked.bool( False ),
            tanThetaMax = cms.double( 1.2 ),
            NormChi2Cut3D = cms.double( 10.0 ),
            minHitsPerSegment = cms.int32( 3 ),
            ForceCovarianceAll = cms.bool( False ),
            yweightPenaltyThreshold = cms.double( 1.0 ),
            prePrunLimit = cms.double( 3.17 ),
            hitDropLimit5Hits = cms.double( 0.8 ),
            preClustering = cms.bool( True ),
            prePrun = cms.bool( True ),
            maxDPhi = cms.double( 999.0 ),
            maxDTheta = cms.double( 999.0 ),
            Pruning = cms.bool( True ),
            dYclusBoxMax = cms.double( 8.0 )
          )
        )
      )
    )
)
hltRpcRecHits = cms.EDProducer( "RPCRecHitProducer",
    rpcDigiLabel = cms.InputTag( "simMuonRPCDigis" ),
    recAlgo = cms.string( "RPCRecHitStandardAlgo" ),
    maskSource = cms.string( "File" ),
    maskvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCMaskVec.dat" ),
    deadSource = cms.string( "File" ),
    deadvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCDeadVec.dat" ),
    recAlgoConfig = cms.PSet(  )
)
hltL2MuonSeeds = cms.EDProducer( "L2MuonSeedGenerator",
    InputObjects = cms.InputTag( "l1extraParticles" ),
    GMTReadoutCollection = cms.InputTag( "gmtDigis" ),
    Propagator = cms.string( "SteppingHelixPropagatorAny" ),
    L1MinPt = cms.double( 0.0 ),
    L1MaxEta = cms.double( 2.5 ),
    L1MinQuality = cms.uint32( 1 ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'SteppingHelixPropagatorAny' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    )
)
hltL2Muons = cms.EDProducer( "L2MuonProducer",
    InputObjects = cms.InputTag( "hltL2MuonSeeds" ),
    L2TrajBuilderParameters = cms.PSet( 
      DoRefit = cms.bool( False ),
      SeedPropagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
      FilterParameters = cms.PSet( 
        NumberOfSigma = cms.double( 3.0 ),
        FitDirection = cms.string( "insideOut" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        MaxChi2 = cms.double( 1000.0 ),
        MuonTrajectoryUpdatorParameters = cms.PSet( 
          MaxChi2 = cms.double( 25.0 ),
          RescaleErrorFactor = cms.double( 100.0 ),
          Granularity = cms.int32( 0 ),
          ExcludeRPCFromFit = cms.bool( False ),
          UseInvalidHits = cms.bool( True ),
          RescaleError = cms.bool( False )
        ),
        EnableRPCMeasurement = cms.bool( True ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        EnableDTMeasurement = cms.bool( True ),
        RPCRecSegmentLabel = cms.InputTag( "hltRpcRecHits" ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        EnableCSCMeasurement = cms.bool( True )
      ),
      NavigationType = cms.string( "Standard" ),
      SeedTransformerParameters = cms.PSet( 
        Fitter = cms.string( "hltESPKFFittingSmootherForL2Muon" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        NMinRecHits = cms.uint32( 2 ),
        UseSubRecHits = cms.bool( False ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        RescaleError = cms.double( 100.0 )
      ),
      DoBackwardFilter = cms.bool( True ),
      SeedPosition = cms.string( "in" ),
      BWFilterParameters = cms.PSet( 
        NumberOfSigma = cms.double( 3.0 ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        FitDirection = cms.string( "outsideIn" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        MaxChi2 = cms.double( 100.0 ),
        MuonTrajectoryUpdatorParameters = cms.PSet( 
          MaxChi2 = cms.double( 25.0 ),
          RescaleErrorFactor = cms.double( 100.0 ),
          Granularity = cms.int32( 2 ),
          ExcludeRPCFromFit = cms.bool( False ),
          UseInvalidHits = cms.bool( True ),
          RescaleError = cms.bool( False )
        ),
        EnableRPCMeasurement = cms.bool( True ),
        BWSeedType = cms.string( "fromGenerator" ),
        EnableDTMeasurement = cms.bool( True ),
        RPCRecSegmentLabel = cms.InputTag( "hltRpcRecHits" ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
        EnableCSCMeasurement = cms.bool( True )
      ),
      DoSeedRefit = cms.bool( False )
    ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPFastSteppingHelixPropagatorAny',
        'hltESPFastSteppingHelixPropagatorOpposite' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    TrackLoaderParameters = cms.PSet( 
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
      DoSmoothing = cms.bool( False ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPosition = cms.vdouble( 0.0, 0.0, 0.0 ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 )
      ),
      VertexConstraint = cms.bool( True )
    )
)
hltL2MuonCandidates = cms.EDProducer( "L2MuonCandidateProducer",
    InputObjects = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
)
hltL2Mu10L2Filtered10 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu10L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 10.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1sL1SingleMu12 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu12" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL2Mu20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu12L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu12" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2Mu20L2Filtered20 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu12L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 20.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreL2DoubleMu0 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiMuonL2PreFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDiMuonL1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMuOpenL1Filtered = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpen" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltSingleMu3L2Filtered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuOpenL1Filtered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL3TrajSeedOIState = cms.EDProducer( "TSGFromL2Muon",
    PtCut = cms.double( 1.0 ),
    PCut = cms.double( 2.5 ),
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPSteppingHelixPropagatorOpposite',
        'hltESPSteppingHelixPropagatorAlong' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    MuonTrackingRegionBuilder = cms.PSet(  ),
    TkSeedGenerator = cms.PSet( 
      propagatorCompatibleName = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
      option = cms.uint32( 3 ),
      maxChi2 = cms.double( 40.0 ),
      errorMatrixPset = cms.PSet( 
        atIP = cms.bool( True ),
        action = cms.string( "use" ),
        errorMatrixValuesPSet = cms.PSet( 
          pf3_V12 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V13 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V11 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V14 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V15 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V34 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          yAxis = cms.vdouble( 0.0, 1.0, 1.4, 10.0 ),
          pf3_V33 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V45 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V44 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          xAxis = cms.vdouble( 0.0, 13.0, 30.0, 70.0, 1000.0 ),
          pf3_V23 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V22 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V55 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          zAxis = cms.vdouble( -3.14159, 3.14159 ),
          pf3_V35 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V25 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V24 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          )
        )
      ),
      propagatorName = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
      manySeeds = cms.bool( False ),
      copyMuonRecHit = cms.bool( False ),
      ComponentName = cms.string( "TSGForRoadSearch" )
    ),
    TrackerSeedCleaner = cms.PSet(  ),
    TSGFromMixedPairs = cms.PSet(  ),
    TSGFromPixelTriplets = cms.PSet(  ),
    TSGFromPixelPairs = cms.PSet(  ),
    TSGForRoadSearchOI = cms.PSet(  ),
    TSGForRoadSearchIOpxl = cms.PSet(  ),
    TSGFromPropagation = cms.PSet(  ),
    TSGFromCombinedHits = cms.PSet(  )
)
hltL3TkTracksFromL2OIState = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "" ),
    Fitter = cms.string( "hltESPKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL3TrackCandidateFromL2OIState" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltL3MuonsOIState = cms.EDProducer( "L3MuonProducer",
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    L3TrajBuilderParameters = cms.PSet( 
      ScaleTECyFactor = cms.double( -1.0 ),
      GlbRefitterParameters = cms.PSet( 
        TrackerSkipSection = cms.int32( -1 ),
        DoPredictionsOnly = cms.bool( False ),
        PropDirForCosmics = cms.bool( False ),
        HitThreshold = cms.int32( 1 ),
        MuonHitsOption = cms.int32( 1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        TrackerSkipSystem = cms.int32( -1 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      MuonTrackingRegionBuilder = cms.PSet( 
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        OnDemand = cms.double( -1.0 ),
        Rescale_Dz = cms.double( 3.0 ),
        Eta_min = cms.double( 0.05 ),
        Rescale_phi = cms.double( 3.0 ),
        Eta_fixed = cms.double( 0.2 ),
        DeltaZ_Region = cms.double( 15.9 ),
        MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Phi_fixed = cms.double( 0.2 ),
        DeltaR = cms.double( 0.2 ),
        EscapePt = cms.double( 1.5 ),
        UseFixedRegion = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        Rescale_eta = cms.double( 3.0 ),
        Phi_min = cms.double( 0.05 ),
        UseVertex = cms.bool( False ),
        beamSpot = cms.InputTag( "offlineBeamSpot" )
      ),
      RefitRPCHits = cms.bool( True ),
      PCut = cms.double( 2.5 ),
      TrackTransformer = cms.PSet( 
        DoPredictionsOnly = cms.bool( False ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" )
      ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        MinP = cms.double( 2.5 ),
        MinPt = cms.double( 1.0 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        LocChi2Cut = cms.double( 0.0010 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_3 = cms.double( 7.0 ),
        Quality_2 = cms.double( 15.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        Quality_1 = cms.double( 20.0 )
      ),
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      tkTrajLabel = cms.InputTag( "hltL3TkTracksFromL2OIState" )
    ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    TrackLoaderParameters = cms.PSet( 
      PutTkTrackIntoEvent = cms.untracked.bool( True ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      SmoothTkTrack = cms.untracked.bool( False ),
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 )
      ),
      VertexConstraint = cms.bool( False ),
      DoSmoothing = cms.bool( True )
    )
)
hltL3TkTracksFromL2OIHit = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "" ),
    Fitter = cms.string( "hltESPKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL3TrackCandidateFromL2OIHit" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltL3MuonsOIHit = cms.EDProducer( "L3MuonProducer",
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    L3TrajBuilderParameters = cms.PSet( 
      ScaleTECyFactor = cms.double( -1.0 ),
      GlbRefitterParameters = cms.PSet( 
        TrackerSkipSection = cms.int32( -1 ),
        DoPredictionsOnly = cms.bool( False ),
        PropDirForCosmics = cms.bool( False ),
        HitThreshold = cms.int32( 1 ),
        MuonHitsOption = cms.int32( 1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        TrackerSkipSystem = cms.int32( -1 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      MuonTrackingRegionBuilder = cms.PSet( 
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        OnDemand = cms.double( -1.0 ),
        Rescale_Dz = cms.double( 3.0 ),
        Eta_min = cms.double( 0.05 ),
        Rescale_phi = cms.double( 3.0 ),
        Eta_fixed = cms.double( 0.2 ),
        DeltaZ_Region = cms.double( 15.9 ),
        MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Phi_fixed = cms.double( 0.2 ),
        DeltaR = cms.double( 0.2 ),
        EscapePt = cms.double( 1.5 ),
        UseFixedRegion = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        Rescale_eta = cms.double( 3.0 ),
        Phi_min = cms.double( 0.05 ),
        UseVertex = cms.bool( False ),
        beamSpot = cms.InputTag( "offlineBeamSpot" )
      ),
      RefitRPCHits = cms.bool( True ),
      PCut = cms.double( 2.5 ),
      TrackTransformer = cms.PSet( 
        DoPredictionsOnly = cms.bool( False ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" )
      ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        MinP = cms.double( 2.5 ),
        MinPt = cms.double( 1.0 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        LocChi2Cut = cms.double( 0.0010 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_3 = cms.double( 7.0 ),
        Quality_2 = cms.double( 15.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        Quality_1 = cms.double( 20.0 )
      ),
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      tkTrajLabel = cms.InputTag( "hltL3TkTracksFromL2OIHit" )
    ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    TrackLoaderParameters = cms.PSet( 
      PutTkTrackIntoEvent = cms.untracked.bool( True ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      SmoothTkTrack = cms.untracked.bool( False ),
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 )
      ),
      VertexConstraint = cms.bool( False ),
      DoSmoothing = cms.bool( True )
    )
)
hltL3TkFromL2OICombination = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3MuonsOIState','hltL3MuonsOIHit' )
)
hltL3TkTracksFromL2IOHit = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "" ),
    Fitter = cms.string( "hltESPKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL3TrackCandidateFromL2IOHit" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltL3MuonsIOHit = cms.EDProducer( "L3MuonProducer",
    MuonCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    L3TrajBuilderParameters = cms.PSet( 
      ScaleTECyFactor = cms.double( -1.0 ),
      GlbRefitterParameters = cms.PSet( 
        TrackerSkipSection = cms.int32( -1 ),
        DoPredictionsOnly = cms.bool( False ),
        PropDirForCosmics = cms.bool( False ),
        HitThreshold = cms.int32( 1 ),
        MuonHitsOption = cms.int32( 1 ),
        Chi2CutRPC = cms.double( 1.0 ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        DTRecSegmentLabel = cms.InputTag( "hltDt4DSegments" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        CSCRecSegmentLabel = cms.InputTag( "hltCscSegments" ),
        Chi2CutCSC = cms.double( 150.0 ),
        Chi2CutDT = cms.double( 10.0 ),
        RefitRPCHits = cms.bool( True ),
        SkipStation = cms.int32( -1 ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" ),
        TrackerSkipSystem = cms.int32( -1 )
      ),
      ScaleTECxFactor = cms.double( -1.0 ),
      TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
      MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
      MuonTrackingRegionBuilder = cms.PSet( 
        EtaR_UpperLimit_Par1 = cms.double( 0.25 ),
        EtaR_UpperLimit_Par2 = cms.double( 0.15 ),
        OnDemand = cms.double( -1.0 ),
        Rescale_Dz = cms.double( 3.0 ),
        Eta_min = cms.double( 0.05 ),
        Rescale_phi = cms.double( 3.0 ),
        Eta_fixed = cms.double( 0.2 ),
        DeltaZ_Region = cms.double( 15.9 ),
        MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        PhiR_UpperLimit_Par2 = cms.double( 0.2 ),
        vertexCollection = cms.InputTag( "pixelVertices" ),
        Phi_fixed = cms.double( 0.2 ),
        DeltaR = cms.double( 0.2 ),
        EscapePt = cms.double( 1.5 ),
        UseFixedRegion = cms.bool( False ),
        PhiR_UpperLimit_Par1 = cms.double( 0.6 ),
        Rescale_eta = cms.double( 3.0 ),
        Phi_min = cms.double( 0.05 ),
        UseVertex = cms.bool( False ),
        beamSpot = cms.InputTag( "offlineBeamSpot" )
      ),
      RefitRPCHits = cms.bool( True ),
      PCut = cms.double( 2.5 ),
      TrackTransformer = cms.PSet( 
        DoPredictionsOnly = cms.bool( False ),
        Fitter = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
        TrackerRecHitBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
        Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
        MuonRecHitBuilder = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
        RefitDirection = cms.string( "insideOut" ),
        RefitRPCHits = cms.bool( True ),
        Propagator = cms.string( "hltESPSmartPropagatorAny" )
      ),
      GlobalMuonTrackMatcher = cms.PSet( 
        Pt_threshold1 = cms.double( 0.0 ),
        DeltaDCut_3 = cms.double( 15.0 ),
        MinP = cms.double( 2.5 ),
        MinPt = cms.double( 1.0 ),
        Chi2Cut_1 = cms.double( 50.0 ),
        Pt_threshold2 = cms.double( 9.99999999E8 ),
        LocChi2Cut = cms.double( 0.0010 ),
        Eta_threshold = cms.double( 1.2 ),
        Quality_3 = cms.double( 7.0 ),
        Quality_2 = cms.double( 15.0 ),
        Chi2Cut_2 = cms.double( 50.0 ),
        Chi2Cut_3 = cms.double( 200.0 ),
        DeltaDCut_1 = cms.double( 40.0 ),
        DeltaRCut_2 = cms.double( 0.2 ),
        DeltaRCut_3 = cms.double( 1.0 ),
        DeltaDCut_2 = cms.double( 10.0 ),
        DeltaRCut_1 = cms.double( 0.1 ),
        Propagator = cms.string( "hltESPSmartPropagator" ),
        Quality_1 = cms.double( 20.0 )
      ),
      PtCut = cms.double( 1.0 ),
      TrackerPropagator = cms.string( "SteppingHelixPropagatorAny" ),
      tkTrajLabel = cms.InputTag( "hltL3TkTracksFromL2IOHit" )
    ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPSmartPropagatorAny',
        'SteppingHelixPropagatorAny',
        'hltESPSmartPropagator',
        'hltESPSteppingHelixPropagatorOpposite' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    TrackLoaderParameters = cms.PSet( 
      PutTkTrackIntoEvent = cms.untracked.bool( True ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      SmoothTkTrack = cms.untracked.bool( False ),
      MuonSeededTracksInstance = cms.untracked.string( "L2Seeded" ),
      Smoother = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        Propagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 )
      ),
      VertexConstraint = cms.bool( False ),
      DoSmoothing = cms.bool( True )
    )
)
hltL3TrajectorySeed = cms.EDProducer( "L3MuonTrajectorySeedCombiner",
    labels = cms.VInputTag( 'hltL3TrajSeedIOHit','hltL3TrajSeedOIState','hltL3TrajSeedOIHit' )
)
hltL3TrackCandidateFromL2 = cms.EDProducer( "L3TrackCandCombiner",
    labels = cms.VInputTag( 'hltL3TrackCandidateFromL2IOHit','hltL3TrackCandidateFromL2OIHit','hltL3TrackCandidateFromL2OIState' )
)
hltL3TkTracksFromL2 = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3TkTracksFromL2IOHit','hltL3TkTracksFromL2OIHit','hltL3TkTracksFromL2OIState' )
)
hltL3MuonsLinksCombination = cms.EDProducer( "L3TrackLinksCombiner",
    labels = cms.VInputTag( 'hltL3MuonsOIState','hltL3MuonsOIHit','hltL3MuonsIOHit' )
)
hltL3Muons = cms.EDProducer( "L3TrackCombiner",
    labels = cms.VInputTag( 'hltL3MuonsOIState','hltL3MuonsOIHit','hltL3MuonsIOHit' )
)
hltL3MuonCandidates = cms.EDProducer( "L3MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltL3Muons" )
)
hltSingleMu3L3Filtered3 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMu3L2Filtered0" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1sL1SingleMu3 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu3" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreMu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu3L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu3" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltSingleMu5L2Filtered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu3L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMu5L3Filtered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMu5L2Filtered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu8 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL2Mu3L2Filtered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu3L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltSingleMu8L3Filtered8 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu3L2Filtered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 8.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1sL1SingleMu7 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu7" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreMu12 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu7L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu7" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2Mu7L2Filtered7 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu7L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 7.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltSingleMu12L3Filtered12 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu7L2Filtered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMu15L3Filtered15 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu10L2Filtered10" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 15.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMu12L2Filtered12 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu12L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMu20L3Filtered20 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMu12L2Filtered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 20.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu24 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL2Mu12L2Filtered12 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu12L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltSingleMu24L3Filtered24 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu12L2Filtered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 24.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMu30L3Filtered30 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu12L2Filtered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 30.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1sL1SingleMu16 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu16" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreMu40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1SingleMu16L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu16" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2Mu16L2Filtered16 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMu16L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 16.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltSingleMu40L3Filtered40 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu16L2Filtered16" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 40.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreIsoMu12 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTowerMakerForMuons = cms.EDProducer( "CaloTowersCreator",
    EBThreshold = cms.double( 0.07 ),
    EEThreshold = cms.double( 0.3 ),
    UseEtEBTreshold = cms.bool( False ),
    UseEtEETreshold = cms.bool( False ),
    UseSymEBTreshold = cms.bool( False ),
    UseSymEETreshold = cms.bool( False ),
    HcalThreshold = cms.double( -1000.0 ),
    HBThreshold = cms.double( 0.7 ),
    HESThreshold = cms.double( 0.8 ),
    HEDThreshold = cms.double( 0.8 ),
    HOThreshold0 = cms.double( 3.5 ),
    HOThresholdPlus1 = cms.double( 3.5 ),
    HOThresholdMinus1 = cms.double( 3.5 ),
    HOThresholdPlus2 = cms.double( 3.5 ),
    HOThresholdMinus2 = cms.double( 3.5 ),
    HF1Threshold = cms.double( 0.5 ),
    HF2Threshold = cms.double( 0.85 ),
    EBWeight = cms.double( 1.0 ),
    EEWeight = cms.double( 1.0 ),
    HBWeight = cms.double( 1.0 ),
    HESWeight = cms.double( 1.0 ),
    HEDWeight = cms.double( 1.0 ),
    HOWeight = cms.double( 1.0E-99 ),
    HF1Weight = cms.double( 1.0 ),
    HF2Weight = cms.double( 1.0 ),
    EcutTower = cms.double( -1000.0 ),
    EBSumThreshold = cms.double( 0.2 ),
    EESumThreshold = cms.double( 0.45 ),
    UseHO = cms.bool( False ),
    MomConstrMethod = cms.int32( 1 ),
    MomHBDepth = cms.double( 0.2 ),
    MomHEDepth = cms.double( 0.4 ),
    MomEBDepth = cms.double( 0.3 ),
    MomEEDepth = cms.double( 0.0 ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    hoInput = cms.InputTag( "hltHoreco" ),
    hfInput = cms.InputTag( "hltHfreco" ),
    AllowMissingInputs = cms.bool( False ),
    HcalAcceptSeverityLevel = cms.uint32( 9 ),
    EcalAcceptSeverityLevel = cms.uint32( 3 ),
    UseHcalRecoveredHits = cms.bool( False ),
    UseEcalRecoveredHits = cms.bool( False ),
    UseRejectedHitsOnly = cms.bool( False ),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    EcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    UseRejectedRecoveredHcalHits = cms.bool( False ),
    UseRejectedRecoveredEcalHits = cms.bool( False ),
    EBGrid = cms.vdouble(  ),
    EBWeights = cms.vdouble(  ),
    EEGrid = cms.vdouble(  ),
    EEWeights = cms.vdouble(  ),
    HBGrid = cms.vdouble(  ),
    HBWeights = cms.vdouble(  ),
    HESGrid = cms.vdouble(  ),
    HESWeights = cms.vdouble(  ),
    HEDGrid = cms.vdouble(  ),
    HEDWeights = cms.vdouble(  ),
    HOGrid = cms.vdouble(  ),
    HOWeights = cms.vdouble(  ),
    HF1Grid = cms.vdouble(  ),
    HF1Weights = cms.vdouble(  ),
    HF2Grid = cms.vdouble(  ),
    HF2Weights = cms.vdouble(  ),
    ecalInputs = cms.VInputTag( 'hltEcalRegionalMuonsRecHit:EcalRecHitsEB','hltEcalRegionalMuonsRecHit:EcalRecHitsEE' )
)
hltL2MuonIsolations = cms.EDProducer( "L2MuonIsolationProducer",
    StandAloneCollectionLabel = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
    ExtractorPSet = cms.PSet( 
      DR_Veto_H = cms.double( 0.1 ),
      Vertex_Constraint_Z = cms.bool( False ),
      Threshold_H = cms.double( 0.5 ),
      ComponentName = cms.string( "CaloExtractor" ),
      Threshold_E = cms.double( 0.2 ),
      DR_Max = cms.double( 0.24 ),
      DR_Veto_E = cms.double( 0.07 ),
      Weight_E = cms.double( 1.5 ),
      Vertex_Constraint_XY = cms.bool( False ),
      DepositLabel = cms.untracked.string( "EcalPlusHcal" ),
      CaloTowerCollectionLabel = cms.InputTag( "hltTowerMakerForMuons" ),
      Weight_H = cms.double( 1.0 )
    ),
    IsolatorPSet = cms.PSet( 
      ConeSizes = cms.vdouble( 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24 ),
      ComponentName = cms.string( "SimpleCutsIsolator" ),
      EtaBounds = cms.vdouble( 0.0435, 0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 1.4355, 1.5225, 1.6095, 1.6965, 1.785, 1.88, 1.9865, 2.1075, 2.247, 2.411 ),
      Thresholds = cms.vdouble( 4.0, 3.7, 4.0, 3.5, 3.4, 3.4, 3.2, 3.4, 3.1, 2.9, 2.9, 2.7, 3.1, 3.0, 2.4, 2.1, 2.0, 2.3, 2.2, 2.4, 2.5, 2.5, 2.6, 2.9, 3.1, 2.9 )
    )
)
hltSingleMuIsoL2IsoFiltered7 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu7L2Filtered7" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    DepTag = cms.VInputTag( 'hltL2MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltSingleMuIsoL3PreFiltered12 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL2IsoFiltered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltL3MuonIsolations = cms.EDProducer( "L3MuonIsolationProducer",
    inputMuonCollection = cms.InputTag( "hltL3Muons" ),
    OutputMuIsoDeposits = cms.bool( True ),
    TrackPt_Min = cms.double( -1.0 ),
    CutsPSet = cms.PSet( 
      ConeSizes = cms.vdouble( 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24 ),
      ComponentName = cms.string( "SimpleCuts" ),
      Thresholds = cms.vdouble( 1.1, 1.1, 1.1, 1.1, 1.2, 1.1, 1.2, 1.1, 1.2, 1.0, 1.1, 1.0, 1.0, 1.1, 1.0, 1.0, 1.1, 0.9, 1.1, 0.9, 1.1, 1.0, 1.0, 0.9, 0.8, 0.1 ),
      maxNTracks = cms.int32( -1 ),
      EtaBounds = cms.vdouble( 0.0435, 0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 1.4355, 1.5225, 1.6095, 1.6965, 1.785, 1.88, 1.9865, 2.1075, 2.247, 2.411 ),
      applyCutsORmaxNTracks = cms.bool( False )
    ),
    ExtractorPSet = cms.PSet( 
      Chi2Prob_Min = cms.double( -1.0 ),
      Diff_z = cms.double( 0.2 ),
      inputTrackCollection = cms.InputTag( "hltPixelTracks" ),
      ReferenceRadius = cms.double( 6.0 ),
      BeamSpotLabel = cms.InputTag( "offlineBeamSpot" ),
      ComponentName = cms.string( "PixelTrackExtractor" ),
      DR_Max = cms.double( 0.24 ),
      Diff_r = cms.double( 0.1 ),
      VetoLeadingTrack = cms.bool( True ),
      DR_VetoPt = cms.double( 0.025 ),
      DR_Veto = cms.double( 0.01 ),
      NHits_Min = cms.uint32( 0 ),
      Chi2Ndof_Max = cms.double( 1.0E64 ),
      Pt_Min = cms.double( -1.0 ),
      DepositLabel = cms.untracked.string( "PXLS" ),
      BeamlineOption = cms.string( "BeamSpotFromEvent" ),
      PropagateTracksToRadius = cms.bool( True ),
      PtVeto_Min = cms.double( 2.0 )
    )
)
hltSingleMuIsoL3IsoFiltered12 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL3PreFiltered12" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltPreIsoMu15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMuIsoL2IsoFiltered10 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu10L2Filtered10" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    DepTag = cms.VInputTag( 'hltL2MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltSingleMuIsoL3PreFiltered15 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL2IsoFiltered10" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 15.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMuIsoL3IsoFiltered15 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL3PreFiltered15" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltPreIsoMu17 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMuIsoL3PreFiltered17 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL2IsoFiltered10" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 17.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMuIsoL3IsoFiltered17 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL3PreFiltered17" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltPreIsoMu24 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMuIsoL2IsoFiltered12 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu12L2Filtered12" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    DepTag = cms.VInputTag( 'hltL2MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltSingleMuIsoL3PreFiltered24 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL2IsoFiltered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 24.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMuIsoL3IsoFiltered24 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL3PreFiltered24" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltPreIsoMu30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleMuIsoL3PreFiltered30 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL2IsoFiltered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 30.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMuIsoL3IsoFiltered30 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMuIsoL3PreFiltered30" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltL1sL1DoubleMu3 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleMu3" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL2DoubleMu23NoVertex = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1DoubleMuon3L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu3" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2MuonCandidatesNoVtx = cms.EDProducer( "L2MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltL2Muons" )
)
hltL2DoubleMu23NoVertexL2PreFiltered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidatesNoVtx" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMuon3L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 3.0 ),
    MinNhits = cms.int32( 1 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 23.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiMuonL3PreFiltered3 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDiMuonL2PreFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu6 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiMuon3L2PreFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMuon3L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDiMuonL3PreFiltered6 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDiMuon3L2PreFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 6.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu7 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiMuonL3PreFiltered7 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDiMuon3L2PreFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 7.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu0Bs = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDimuonL1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu0" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltDimuonL2PreFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDimuonL3PreFiltered2Bs = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 2.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDoubleMu2BsL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 4.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 4.8 ),
    MaxInvMass = cms.double( 6.0 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu4Excl = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1DoubleMuon3L1Filtered3 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu3" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 3.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltL2DoubleMu3L2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMuon3L1Filtered3" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 3.0 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDiMuonL3PreFiltered4 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2DoubleMu3L2Filtered" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 4.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDoubleMu4ExclL3PreFiltered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2DoubleMu3L2Filtered" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 0.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 9999.0 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 0.3 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 99999.9 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltPreDimuon0Jpsi = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltJpsiL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 0.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 2.8 ),
    MaxInvMass = cms.double( 3.35 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerJpsi0 = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJpsiL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterJpsi = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerJpsi0" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon0Upsilon = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltUpsilonL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 0.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 8.5 ),
    MaxInvMass = cms.double( 11.5 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 2.5 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerUpsilon = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltUpsilonL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterUpsilon = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerUpsilon" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDoubleMu2BarrelBs = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleMu2BarrelBsL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 1.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 3.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 2.0 ),
    MinInvMass = cms.double( 4.8 ),
    MaxInvMass = cms.double( 6.0 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 1.4 ),
    saveTags = cms.bool( True )
)
hltPreDimuon5BarrelUpsilon = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBarrelUpsilonL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 4.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 8.5 ),
    MaxInvMass = cms.double( 11.5 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 1.25 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerUpsilonBarrel = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltBarrelUpsilonL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterUpsilonBarrel = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerUpsilonBarrel" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDoubleMu2Dimuon6Bs = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleMu2Dimuon6BsL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 5.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 2.0 ),
    MinInvMass = cms.double( 4.8 ),
    MaxInvMass = cms.double( 6.0 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltPreDimuon7LowMassDisplaced = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltLowMassDisplacedL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.2 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 6.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 3.0 ),
    MinInvMass = cms.double( 1.0 ),
    MaxInvMass = cms.double( 4.8 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerLowMass = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltLowMassDisplacedL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltDisplacedmumuFilterLowMass = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 3.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.05 ),
    MinCosinePointingAngle = cms.double( 0.9 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerLowMass" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon7JpsiDisplaced = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltJpsiDisplacedL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 6.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 2.9 ),
    MaxInvMass = cms.double( 3.3 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerJpsi = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJpsiDisplacedL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltDisplacedmumuFilterJpsi = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 3.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( 0.9 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerJpsi" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon7JpsiXBarrel = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltJpsiXBarrelL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 6.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 2.95 ),
    MaxInvMass = cms.double( 3.25 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 1.25 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerJpsiXBarrel = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJpsiXBarrelL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterJpsiXBarrel = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerJpsiXBarrel" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon7PsiPrime = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPsiPrimeL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 6.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 3.35 ),
    MaxInvMass = cms.double( 4.05 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 2.5 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerPsiPrime = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltPsiPrimeL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterPsiPrime = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerPsiPrime" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon10BarrelJpsi = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBarrelJpsiL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 9.9 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 2.8 ),
    MaxInvMass = cms.double( 3.35 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 1.25 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerJpsiBarrel = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltBarrelJpsiL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterJpsiBarrel = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerJpsiBarrel" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon0JpsiMuon = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTripleMuonL1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu0" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 3 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltTripleMuonL2PreFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltTripleMuonL1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 3 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltTripleMuL3PreFiltered0 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltTripleMuonL2PreFiltered0" ),
    MinN = cms.int32( 3 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltJpsiMuonL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 0.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 2.8 ),
    MaxInvMass = cms.double( 3.35 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 999999.0 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerJpsiMuon = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltJpsiMuonL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterJpsiMuon = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerJpsiMuon" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreDimuon0UpsilonMuon = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltUpsilonMuonL3Filtered = cms.EDFilter( "HLTMuonDimuonL3Filter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDimuonL2PreFiltered0" ),
    FastAccept = cms.bool( False ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    ChargeOpt = cms.int32( -1 ),
    MinPtPair = cms.double( 0.0 ),
    MinPtMax = cms.double( 0.0 ),
    MinPtMin = cms.double( 0.0 ),
    MinInvMass = cms.double( 8.5 ),
    MaxInvMass = cms.double( 11.5 ),
    MinAcop = cms.double( -999.0 ),
    MaxAcop = cms.double( 999.0 ),
    MinPtBalance = cms.double( -1.0 ),
    MaxPtBalance = cms.double( 999999.0 ),
    NSigmaPt = cms.double( 0.0 ),
    MaxDCAMuMu = cms.double( 0.5 ),
    MaxRapidityPair = cms.double( 2.5 ),
    saveTags = cms.bool( True )
)
hltDisplacedmumuVtxProducerUpsilonMuon = cms.EDProducer( "HLTDisplacedmumuVtxProducer",
    Src = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltUpsilonMuonL3Filtered" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinPtPair = cms.double( 0.0 ),
    MinInvMass = cms.double( 0.0 ),
    MaxInvMass = cms.double( 999999.0 ),
    ChargeOpt = cms.int32( -1 )
)
hltVertexmumuFilterUpsilonMuon = cms.EDFilter( "HLTDisplacedmumuFilter",
    FastAccept = cms.bool( True ),
    MinLxySignificance = cms.double( 0.0 ),
    MaxNormalisedChi2 = cms.double( 999999.0 ),
    MinVtxProbability = cms.double( 0.0050 ),
    MinCosinePointingAngle = cms.double( -2.0 ),
    saveTags = cms.bool( True ),
    DisplacedVertexTag = cms.InputTag( "hltDisplacedmumuVtxProducerUpsilonMuon" ),
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    MuonTag = cms.InputTag( "hltL3MuonCandidates" )
)
hltPreMu13Mu8 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1DoubleMuon3L2Filtered7 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMuon3L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 7.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltDiMuonL3PreFiltered8 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltDiMuon3L2PreFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 8.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltSingleMu13L3Filtered13 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMuon3L2Filtered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 13.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreTripleMu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1DoubleMu3L1TriMuFiltered3 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu3" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 3.0 ),
    MinN = cms.int32( 3 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1DoubleMu3L2TriMuFiltered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMu3L1TriMuFiltered3" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 3 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1DoubleMu3L3TriMuFiltered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1DoubleMu3L2TriMuFiltered3" ),
    MinN = cms.int32( 3 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu5L2Mu2Jpsi = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMu5L2Mu2L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMu0" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltMu5L2Mu2L2PreFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltMu5L2Mu2L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 2.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltMu5L2Mu2L3Filtered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltMu5L2Mu2L2PreFiltered0" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltMu5L2Mu2JpsiTrackMassFiltered = cms.EDFilter( "HLTMuonTrackMassFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    TrackTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltMu5L2Mu2L3Filtered5" ),
    saveTags = cms.bool( True ),
    checkCharge = cms.bool( True ),
    MinTrackPt = cms.double( 2.0 ),
    MinTrackP = cms.double( 0.0 ),
    MaxTrackEta = cms.double( 999.0 ),
    MaxTrackDxy = cms.double( 999.0 ),
    MaxTrackDz = cms.double( 999.0 ),
    MinTrackHits = cms.int32( 2 ),
    MaxTrackNormChi2 = cms.double( 1.0E10 ),
    MaxDCAMuonTrack = cms.double( 99999.9 ),
    CutCowboys = cms.bool( False ),
    MinMasses = cms.vdouble( 1.8 ),
    MaxMasses = cms.vdouble( 4.5 )
)
hltL1sL1SingleEG12 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG12" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPrePhoton20CaloIdVLIsoL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHybridSuperClustersL1Isolated = cms.EDProducer( "EgammaHLTHybridClusterProducer",
    debugLevel = cms.string( "INFO" ),
    basicclusterCollection = cms.string( "" ),
    superclusterCollection = cms.string( "" ),
    ecalhitproducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalhitcollection = cms.string( "EcalRecHitsEB" ),
    l1TagIsolated = cms.InputTag( 'l1extraParticles','Isolated' ),
    l1TagNonIsolated = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    doIsolated = cms.bool( True ),
    l1LowerThr = cms.double( 5.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 999.0 ),
    regionEtaMargin = cms.double( 0.14 ),
    regionPhiMargin = cms.double( 0.4 ),
    HybridBarrelSeedThr = cms.double( 1.5 ),
    step = cms.int32( 17 ),
    ethresh = cms.double( 0.1 ),
    eseed = cms.double( 0.35 ),
    ewing = cms.double( 0.0 ),
    dynamicEThresh = cms.bool( False ),
    eThreshA = cms.double( 0.0030 ),
    eThreshB = cms.double( 0.1 ),
    excludeFlagged = cms.bool( False ),
    dynamicPhiRoad = cms.bool( False ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    RecHitSeverityToBeExcluded = cms.vint32( 4 ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    ),
    bremRecoveryPset = cms.PSet(  )
)
hltCorrectedHybridSuperClustersL1Isolated = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    rawSuperClusterProducer = cms.InputTag( "hltHybridSuperClustersL1Isolated" ),
    superClusterAlgo = cms.string( "Hybrid" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.03 ),
    etThresh = cms.double( 1.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 1.1 ),
      fBremVec = cms.vdouble( -0.05208, 0.1331, 0.9196, -5.735E-4, 1.343 ),
      brLinearHighThr = cms.double( 8.0 ),
      fEtEtaVec = cms.vdouble( 1.0012, -0.5714, 0.0, 0.0, 0.0, 0.5549, 12.74, 1.0448, 0.0, 0.0, 0.0, 0.0, 8.0, 1.023, -0.00181, 0.0, 0.0 )
    ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet(  )
)
hltMulti5x5BasicClustersL1Isolated = cms.EDProducer( "EgammaHLTMulti5x5ClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    doBarrel = cms.bool( False ),
    doEndcaps = cms.bool( True ),
    doIsolated = cms.bool( True ),
    barrelHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    endcapHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    barrelHitCollection = cms.string( "EcalRecHitsEB" ),
    endcapHitCollection = cms.string( "EcalRecHitsEE" ),
    barrelClusterCollection = cms.string( "notused" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    Multi5x5BarrelSeedThr = cms.double( 0.5 ),
    Multi5x5EndcapSeedThr = cms.double( 0.18 ),
    l1TagIsolated = cms.InputTag( 'l1extraParticles','Isolated' ),
    l1TagNonIsolated = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    l1LowerThr = cms.double( 5.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 999.0 ),
    regionEtaMargin = cms.double( 0.3 ),
    regionPhiMargin = cms.double( 0.4 ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    )
)
hltMulti5x5SuperClustersL1Isolated = cms.EDProducer( "Multi5x5SuperClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    endcapClusterProducer = cms.string( "hltMulti5x5BasicClustersL1Isolated" ),
    barrelClusterProducer = cms.string( "notused" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    barrelClusterCollection = cms.string( "multi5x5BarrelBasicClusters" ),
    endcapSuperclusterCollection = cms.string( "multi5x5EndcapSuperClusters" ),
    barrelSuperclusterCollection = cms.string( "multi5x5BarrelSuperClusters" ),
    doBarrel = cms.bool( False ),
    doEndcaps = cms.bool( True ),
    barrelEtaSearchRoad = cms.double( 0.06 ),
    barrelPhiSearchRoad = cms.double( 0.8 ),
    endcapEtaSearchRoad = cms.double( 0.14 ),
    endcapPhiSearchRoad = cms.double( 0.6 ),
    seedTransverseEnergyThreshold = cms.double( 1.0 ),
    dynamicPhiRoad = cms.bool( False ),
    bremRecoveryPset = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet( 
        a = cms.double( 47.85 ),
        c = cms.double( 0.1201 ),
        b = cms.double( 108.8 )
      ),
      doEndcaps = cms.bool( True ),
      doBarrel = cms.bool( False )
    )
)
hltMulti5x5EndcapSuperClustersWithPreshowerL1Isolated = cms.EDProducer( "PreshowerClusterProducer",
    preshRecHitProducer = cms.InputTag( 'hltESRegionalEgammaRecHit','EcalRecHitsES' ),
    endcapSClusterProducer = cms.InputTag( 'hltMulti5x5SuperClustersL1Isolated','multi5x5EndcapSuperClusters' ),
    preshClusterCollectionX = cms.string( "preshowerXClusters" ),
    preshClusterCollectionY = cms.string( "preshowerYClusters" ),
    preshNclust = cms.int32( 4 ),
    etThresh = cms.double( 5.0 ),
    assocSClusterCollection = cms.string( "" ),
    preshStripEnergyCut = cms.double( 0.0 ),
    preshSeededNstrip = cms.int32( 15 ),
    preshClusterEnergyCut = cms.double( 0.0 ),
    debugLevel = cms.string( "" )
)
hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1Isolated = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' ),
    rawSuperClusterProducer = cms.InputTag( "hltMulti5x5EndcapSuperClustersWithPreshowerL1Isolated" ),
    superClusterAlgo = cms.string( "Multi5x5" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 1.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet(  ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 0.6 ),
      fBremVec = cms.vdouble( -0.04163, 0.08552, 0.95048, -0.002308, 1.077 ),
      brLinearHighThr = cms.double( 6.0 ),
      fEtEtaVec = cms.vdouble( 0.9746, -6.512, 0.0, 0.0, 0.02771, 4.983, 0.0, 0.0, -0.007288, -0.9446, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 )
    )
)
hltHybridSuperClustersL1NonIsolated = cms.EDProducer( "EgammaHLTHybridClusterProducer",
    debugLevel = cms.string( "INFO" ),
    basicclusterCollection = cms.string( "" ),
    superclusterCollection = cms.string( "" ),
    ecalhitproducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalhitcollection = cms.string( "EcalRecHitsEB" ),
    l1TagIsolated = cms.InputTag( 'l1extraParticles','Isolated' ),
    l1TagNonIsolated = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    doIsolated = cms.bool( False ),
    l1LowerThr = cms.double( 5.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 999.0 ),
    regionEtaMargin = cms.double( 0.14 ),
    regionPhiMargin = cms.double( 0.4 ),
    HybridBarrelSeedThr = cms.double( 1.5 ),
    step = cms.int32( 17 ),
    ethresh = cms.double( 0.1 ),
    eseed = cms.double( 0.35 ),
    ewing = cms.double( 0.0 ),
    dynamicEThresh = cms.bool( False ),
    eThreshA = cms.double( 0.0030 ),
    eThreshB = cms.double( 0.1 ),
    excludeFlagged = cms.bool( False ),
    dynamicPhiRoad = cms.bool( False ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    RecHitSeverityToBeExcluded = cms.vint32( 4 ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    ),
    bremRecoveryPset = cms.PSet(  )
)
hltCorrectedHybridSuperClustersL1NonIsolatedTemp = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    rawSuperClusterProducer = cms.InputTag( "hltHybridSuperClustersL1NonIsolated" ),
    superClusterAlgo = cms.string( "Hybrid" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.03 ),
    etThresh = cms.double( 1.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 1.1 ),
      fBremVec = cms.vdouble( -0.05208, 0.1331, 0.9196, -5.735E-4, 1.343 ),
      brLinearHighThr = cms.double( 8.0 ),
      fEtEtaVec = cms.vdouble( 1.0012, -0.5714, 0.0, 0.0, 0.0, 0.5549, 12.74, 1.0448, 0.0, 0.0, 0.0, 0.0, 8.0, 1.023, -0.00181, 0.0, 0.0 )
    ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet(  )
)
hltCorrectedHybridSuperClustersL1NonIsolated = cms.EDProducer( "EgammaHLTRemoveDuplicatedSC",
    L1NonIsoUskimmedSC = cms.InputTag( "hltCorrectedHybridSuperClustersL1NonIsolatedTemp" ),
    L1IsoSC = cms.InputTag( "hltCorrectedHybridSuperClustersL1Isolated" ),
    L1NonIsoSkimmedCollection = cms.string( "" )
)
hltMulti5x5BasicClustersL1NonIsolated = cms.EDProducer( "EgammaHLTMulti5x5ClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    doBarrel = cms.bool( False ),
    doEndcaps = cms.bool( True ),
    doIsolated = cms.bool( False ),
    barrelHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    endcapHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    barrelHitCollection = cms.string( "EcalRecHitsEB" ),
    endcapHitCollection = cms.string( "EcalRecHitsEE" ),
    barrelClusterCollection = cms.string( "notused" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    Multi5x5BarrelSeedThr = cms.double( 0.5 ),
    Multi5x5EndcapSeedThr = cms.double( 0.18 ),
    l1TagIsolated = cms.InputTag( 'l1extraParticles','Isolated' ),
    l1TagNonIsolated = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    l1LowerThr = cms.double( 5.0 ),
    l1UpperThr = cms.double( 999.0 ),
    l1LowerThrIgnoreIsolation = cms.double( 999.0 ),
    regionEtaMargin = cms.double( 0.3 ),
    regionPhiMargin = cms.double( 0.4 ),
    RecHitFlagToBeExcluded = cms.vint32(  ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    )
)
hltMulti5x5SuperClustersL1NonIsolated = cms.EDProducer( "Multi5x5SuperClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    endcapClusterProducer = cms.string( "hltMulti5x5BasicClustersL1NonIsolated" ),
    barrelClusterProducer = cms.string( "notused" ),
    endcapClusterCollection = cms.string( "multi5x5EndcapBasicClusters" ),
    barrelClusterCollection = cms.string( "multi5x5BarrelBasicClusters" ),
    endcapSuperclusterCollection = cms.string( "multi5x5EndcapSuperClusters" ),
    barrelSuperclusterCollection = cms.string( "multi5x5BarrelSuperClusters" ),
    doBarrel = cms.bool( False ),
    doEndcaps = cms.bool( True ),
    barrelEtaSearchRoad = cms.double( 0.06 ),
    barrelPhiSearchRoad = cms.double( 0.8 ),
    endcapEtaSearchRoad = cms.double( 0.14 ),
    endcapPhiSearchRoad = cms.double( 0.6 ),
    seedTransverseEnergyThreshold = cms.double( 1.0 ),
    dynamicPhiRoad = cms.bool( False ),
    bremRecoveryPset = cms.PSet( 
      barrel = cms.PSet(  ),
      endcap = cms.PSet( 
        a = cms.double( 47.85 ),
        c = cms.double( 0.1201 ),
        b = cms.double( 108.8 )
      ),
      doEndcaps = cms.bool( True ),
      doBarrel = cms.bool( False )
    )
)
hltMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated = cms.EDProducer( "PreshowerClusterProducer",
    preshRecHitProducer = cms.InputTag( 'hltESRegionalEgammaRecHit','EcalRecHitsES' ),
    endcapSClusterProducer = cms.InputTag( 'hltMulti5x5SuperClustersL1NonIsolated','multi5x5EndcapSuperClusters' ),
    preshClusterCollectionX = cms.string( "preshowerXClusters" ),
    preshClusterCollectionY = cms.string( "preshowerYClusters" ),
    preshNclust = cms.int32( 4 ),
    etThresh = cms.double( 5.0 ),
    assocSClusterCollection = cms.string( "" ),
    preshStripEnergyCut = cms.double( 0.0 ),
    preshSeededNstrip = cms.int32( 15 ),
    preshClusterEnergyCut = cms.double( 0.0 ),
    debugLevel = cms.string( "" )
)
hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolatedTemp = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' ),
    rawSuperClusterProducer = cms.InputTag( "hltMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated" ),
    superClusterAlgo = cms.string( "Multi5x5" ),
    applyEnergyCorrection = cms.bool( True ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 1.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet(  ),
    isl_fCorrPset = cms.PSet(  ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 0.6 ),
      fBremVec = cms.vdouble( -0.04163, 0.08552, 0.95048, -0.002308, 1.077 ),
      brLinearHighThr = cms.double( 6.0 ),
      fEtEtaVec = cms.vdouble( 0.9746, -6.512, 0.0, 0.0, 0.02771, 4.983, 0.0, 0.0, -0.007288, -0.9446, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0 )
    )
)
hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated = cms.EDProducer( "EgammaHLTRemoveDuplicatedSC",
    L1NonIsoUskimmedSC = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolatedTemp" ),
    L1IsoSC = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1Isolated" ),
    L1NonIsoSkimmedCollection = cms.string( "" )
)
hltL1IsoRecoEcalCandidate = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( "hltCorrectedHybridSuperClustersL1Isolated" ),
    scIslandEndcapProducer = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1Isolated" ),
    recoEcalCandidateCollection = cms.string( "" )
)
hltL1NonIsoRecoEcalCandidate = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( "hltCorrectedHybridSuperClustersL1NonIsolated" ),
    scIslandEndcapProducer = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated" ),
    recoEcalCandidateCollection = cms.string( "" )
)
hltEGRegionalL1SingleEG12 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG12" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG20EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG12" ),
    etcutEB = cms.double( 20.0 ),
    etcutEE = cms.double( 20.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsoHLTClusterShape = cms.EDProducer( "EgammaHLTClusterShapeProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' ),
    isIeta = cms.bool( True )
)
hltL1NonIsoHLTClusterShape = cms.EDProducer( "EgammaHLTClusterShapeProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' ),
    isIeta = cms.bool( True )
)
hltEG20CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsolatedPhotonEcalIsol = cms.EDProducer( "EgammaHLTEcalRecIsolationProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    ecalBarrelRecHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalBarrelRecHitCollection = cms.InputTag( "EcalRecHitsEB" ),
    ecalEndcapRecHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalEndcapRecHitCollection = cms.InputTag( "EcalRecHitsEE" ),
    etMinBarrel = cms.double( -9999.0 ),
    eMinBarrel = cms.double( 0.08 ),
    etMinEndcap = cms.double( 0.1 ),
    eMinEndcap = cms.double( -9999.0 ),
    intRadiusBarrel = cms.double( 3.0 ),
    intRadiusEndcap = cms.double( 3.0 ),
    extRadius = cms.double( 0.3 ),
    jurassicWidth = cms.double( 3.0 ),
    useIsolEt = cms.bool( True ),
    tryBoth = cms.bool( True ),
    subtract = cms.bool( False ),
    useNumCrystals = cms.bool( True )
)
hltL1NonIsolatedPhotonEcalIsol = cms.EDProducer( "EgammaHLTEcalRecIsolationProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    ecalBarrelRecHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalBarrelRecHitCollection = cms.InputTag( "EcalRecHitsEB" ),
    ecalEndcapRecHitProducer = cms.InputTag( "hltEcalRegionalEgammaRecHit" ),
    ecalEndcapRecHitCollection = cms.InputTag( "EcalRecHitsEE" ),
    etMinBarrel = cms.double( -9999.0 ),
    eMinBarrel = cms.double( 0.08 ),
    etMinEndcap = cms.double( 0.1 ),
    eMinEndcap = cms.double( -9999.0 ),
    intRadiusBarrel = cms.double( 3.0 ),
    intRadiusEndcap = cms.double( 3.0 ),
    extRadius = cms.double( 0.3 ),
    jurassicWidth = cms.double( 3.0 ),
    useIsolEt = cms.bool( True ),
    tryBoth = cms.bool( True ),
    subtract = cms.bool( False ),
    useNumCrystals = cms.bool( True )
)
hltPhoton20CaloIdVLIsoLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG20CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.5 ),
    thrRegularEE = cms.double( 5.5 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsolatedPhotonHcalForHE = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.0 ),
    outerCone = cms.double( 0.14 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( False )
)
hltL1NonIsolatedPhotonHcalForHE = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.0 ),
    outerCone = cms.double( 0.14 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( False )
)
hltPhoton20CaloIdVLIsoLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVLIsoLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsolatedPhotonHcalIsol = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.16 ),
    outerCone = cms.double( 0.29 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( True )
)
hltL1NonIsolatedPhotonHcalIsol = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.16 ),
    outerCone = cms.double( 0.29 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( True )
)
hltPhoton20CaloIdVLIsoLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVLIsoLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsolatedPhotonHollowTrackIsol = cms.EDProducer( "EgammaHLTPhotonTrackIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    trackProducer = cms.InputTag( "hltL1IsoEgammaRegionalCTFFinalFitWithMaterial" ),
    countTracks = cms.bool( False ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.29 ),
    egTrkIsoZSpan = cms.double( 999999.0 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.06 ),
    egTrkIsoStripBarrel = cms.double( 0.03 ),
    egTrkIsoStripEndcap = cms.double( 0.03 )
)
hltL1NonIsolatedPhotonHollowTrackIsol = cms.EDProducer( "EgammaHLTPhotonTrackIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    trackProducer = cms.InputTag( "hltL1NonIsoEgammaRegionalCTFFinalFitWithMaterial" ),
    countTracks = cms.bool( False ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.29 ),
    egTrkIsoZSpan = cms.double( 999999.0 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.06 ),
    egTrkIsoStripBarrel = cms.double( 0.03 ),
    egTrkIsoStripEndcap = cms.double( 0.03 )
)
hltPhoton20CaloIdVLIsoLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVLIsoLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton20R9IdPhoton18R9Id = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton20R9IdPhoton18R9IdEgammaLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsoR9ID = cms.EDProducer( "EgammaHLTR9IDProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' )
)
hltL1NonIsoR9ID = cms.EDProducer( "EgammaHLTR9IDProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRegionalEgammaRecHit','EcalRecHitsEE' )
)
hltPhoton20R9IdPhoton18R9IdEgammaR9IDFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton20R9IdPhoton18R9IdEgammaLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoR9ID" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoR9ID" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG18EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 18.0 ),
    etcutEE = cms.double( 18.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityPhotonHcalForHE = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.0 ),
    outerCone = cms.double( 0.14 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( False )
)
hltDoubleIsoEG18HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityR9ID = cms.EDProducer( "EgammaHLTR9IDProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' )
)
hltPhoton20R9IdPhoton18R9IdEgammaR9IDDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityR9ID" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton20CaloIdVTIsoTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.0 ),
    thrRegularEE = cms.double( 5.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.0 ),
    thrRegularEE = cms.double( 3.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.0 ),
    thrRegularEE = cms.double( 3.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG8EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityPhotonClusterShape = cms.EDProducer( "EgammaHLTClusterShapeProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    ecalRechitEB = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    ecalRechitEE = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' ),
    isIeta = cms.bool( True )
)
hltEle8CaloIdLCaloIsoVLNoL1SeedClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityPhotonEcalIsol = cms.EDProducer( "EgammaHLTEcalRecIsolationProducer",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    ecalBarrelRecHitProducer = cms.InputTag( "hltEcalRecHitAll" ),
    ecalBarrelRecHitCollection = cms.InputTag( "EcalRecHitsEB" ),
    ecalEndcapRecHitProducer = cms.InputTag( "hltEcalRecHitAll" ),
    ecalEndcapRecHitCollection = cms.InputTag( "EcalRecHitsEE" ),
    etMinBarrel = cms.double( -9999.0 ),
    eMinBarrel = cms.double( 0.08 ),
    etMinEndcap = cms.double( 0.1 ),
    eMinEndcap = cms.double( -9999.0 ),
    intRadiusBarrel = cms.double( 3.0 ),
    intRadiusEndcap = cms.double( 3.0 ),
    extRadius = cms.double( 0.3 ),
    jurassicWidth = cms.double( 3.0 ),
    useIsolEt = cms.bool( True ),
    tryBoth = cms.bool( True ),
    subtract = cms.bool( False ),
    useNumCrystals = cms.bool( True )
)
hltEle8CaloIdLCaloIsoVLNoL1SeedEcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLNoL1SeedClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle8CaloIdLCaloIsoVLNoL1SeedHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLNoL1SeedEcalIsolFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityPhotonHcalIsol = cms.EDProducer( "EgammaHLTHcalIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    hbheRecHitProducer = cms.InputTag( "hltHbhereco" ),
    eMinHB = cms.double( 0.7 ),
    eMinHE = cms.double( 0.8 ),
    etMinHB = cms.double( -1.0 ),
    etMinHE = cms.double( -1.0 ),
    innerCone = cms.double( 0.16 ),
    outerCone = cms.double( 0.29 ),
    depth = cms.int32( -1 ),
    doEtSum = cms.bool( True )
)
hltEle8CaloIdLCaloIsoVLNoL1SeedHcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLNoL1SeedHEFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityStartUpElectronPixelSeeds = cms.EDProducer( "ElectronSeedProducer",
    barrelSuperClusters = cms.InputTag( "hltCorrectedHybridSuperClustersActivity" ),
    endcapSuperClusters = cms.InputTag( "hltCorrectedMulti5x5SuperClustersWithPreshowerActivity" ),
    SeedConfiguration = cms.PSet( 
      searchInTIDTEC = cms.bool( True ),
      HighPtThreshold = cms.double( 35.0 ),
      r2MinF = cms.double( -0.15 ),
      OrderedHitsFactoryPSet = cms.PSet( 
        maxElement = cms.uint32( 0 ),
        ComponentName = cms.string( "StandardHitPairGenerator" ),
        SeedingLayers = cms.string( "hltESPMixedLayerPairs" ),
        useOnDemandTracker = cms.untracked.int32( 0 )
      ),
      DeltaPhi1Low = cms.double( 0.23 ),
      DeltaPhi1High = cms.double( 0.08 ),
      ePhiMin1 = cms.double( -0.08 ),
      PhiMin2 = cms.double( -0.0040 ),
      LowPtThreshold = cms.double( 3.0 ),
      RegionPSet = cms.PSet( 
        deltaPhiRegion = cms.double( 0.4 ),
        originHalfLength = cms.double( 15.0 ),
        useZInVertex = cms.bool( True ),
        deltaEtaRegion = cms.double( 0.1 ),
        ptMin = cms.double( 1.5 ),
        originRadius = cms.double( 0.2 ),
        VertexProducer = cms.InputTag( "dummyVertices" )
      ),
      maxHOverE = cms.double( 999999.0 ),
      dynamicPhiRoad = cms.bool( False ),
      ePhiMax1 = cms.double( 0.04 ),
      DeltaPhi2 = cms.double( 0.0040 ),
      measurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
      SizeWindowENeg = cms.double( 0.675 ),
      nSigmasDeltaZ1 = cms.double( 5.0 ),
      rMaxI = cms.double( 0.2 ),
      rMinI = cms.double( -0.2 ),
      preFilteredSeeds = cms.bool( False ),
      r2MaxF = cms.double( 0.15 ),
      pPhiMin1 = cms.double( -0.04 ),
      initialSeeds = cms.InputTag( "globalPixelSeeds:GlobalPixel" ),
      pPhiMax1 = cms.double( 0.08 ),
      hbheModule = cms.string( "hbhereco" ),
      SCEtCut = cms.double( 3.0 ),
      z2MaxB = cms.double( 0.09 ),
      fromTrackerSeeds = cms.bool( True ),
      hcalRecHits = cms.InputTag( "hltHbhereco" ),
      z2MinB = cms.double( -0.09 ),
      hbheInstance = cms.string( "" ),
      PhiMax2 = cms.double( 0.0040 ),
      hOverEConeSize = cms.double( 0.0 ),
      hOverEHBMinE = cms.double( 999999.0 ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      applyHOverECut = cms.bool( False ),
      hOverEHFMinE = cms.double( 999999.0 )
    )
)
hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLNoL1SeedHcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLDoubleLegCombFilter = cms.EDFilter( "HLTEgammaDoubleLegCombFilter",
    firstLegLastFilter = cms.InputTag( "hltPhoton20CaloIdVTIsoTTrackIsoFilter" ),
    secondLegLastFilter = cms.InputTag( "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" ),
    nrRequiredFirstLeg = cms.int32( 1 ),
    nrRequiredSecondLeg = cms.int32( 1 ),
    maxMatchDR = cms.double( 0.3 )
)
hltPrePhoton26Photon18 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG26EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG12" ),
    etcutEB = cms.double( 26.0 ),
    etcutEE = cms.double( 26.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26Photon18EgammaLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG26EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG18HELastFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton26IsoVLPhoton18 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG26HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG26EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26IsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG26HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26IsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26IsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26IsoVLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26IsoVLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton26IsoVLPhoton18IsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton26IsoVLHcalIsoLastFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26IsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26IsoVLPhoton18IsoVLEgammaEcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26IsoVLPhoton18IsoVLEgammaHcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26IsoVLPhoton18IsoVLEgammaEcalIsolDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltActivityPhotonHollowTrackIsol = cms.EDProducer( "EgammaHLTPhotonTrackIsolationProducersRegional",
    recoEcalCandidateProducer = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    trackProducer = cms.InputTag( "hltEcalActivityEgammaRegionalCTFFinalFitWithMaterial" ),
    countTracks = cms.bool( False ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.29 ),
    egTrkIsoZSpan = cms.double( 999999.0 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.06 ),
    egTrkIsoStripBarrel = cms.double( 0.03 ),
    egTrkIsoStripEndcap = cms.double( 0.03 )
)
hltPhoton26IsoVLPhoton18IsoVLEgammaTrackIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26IsoVLPhoton18IsoVLEgammaHcalIsolDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton26CaloIdLIsoVLPhoton18 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG26CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG26HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG26CaloIdLIsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG26CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG26CaloIdLIsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG26CaloIdLIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG26CaloIdLIsoVLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG26CaloIdLIsoVLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton26CaloIdLIsoVLPhoton18R9Id = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR9IdFilterUnseeded = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityR9ID" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaClusterShapeFilterUnseeded = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaEcalIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaClusterShapeFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaHcalIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaEcalIsolFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaTrackIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaHcalIsolFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaDoubleLegCombFilter = cms.EDFilter( "HLTEgammaDoubleLegCombFilter",
    firstLegLastFilter = cms.InputTag( "hltR9IdFilterUnseeded" ),
    secondLegLastFilter = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaTrackIsolFilterUnseeded" ),
    nrRequiredFirstLeg = cms.int32( 1 ),
    nrRequiredSecondLeg = cms.int32( 1 ),
    maxMatchDR = cms.double( 0.01 )
)
hltPrePhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG26CaloIdLIsoVLHcalIsoLastFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG26CaloIdLIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaEcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaClusterShapeDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaHcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaEcalIsolDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaTrackIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaHcalIsolDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton26R9IdPhoton18CaloIdLIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG26R9IdFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG26HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoR9ID" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoR9ID" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaClusterShapeFilterUnseeded = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaEcalIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaClusterShapeFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaHcalIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaEcalIsolFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaTrackIsolFilterUnseeded = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaHcalIsolFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaDoubleLegCombFilter = cms.EDFilter( "HLTEgammaDoubleLegCombFilter",
    firstLegLastFilter = cms.InputTag( "hltR9IdFilterUnseeded" ),
    secondLegLastFilter = cms.InputTag( "hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaTrackIsolFilterUnseeded" ),
    nrRequiredFirstLeg = cms.int32( 1 ),
    nrRequiredSecondLeg = cms.int32( 1 ),
    maxMatchDR = cms.double( 0.01 )
)
hltPrePhoton26R9IdPhoton18R9Id = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton26R9IdPhoton18R9IdEgammaLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26R9IdPhoton18R9IdEgammaR9IDFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton26R9IdPhoton18R9IdEgammaLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoR9ID" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoR9ID" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton26R9IdPhoton18R9IdEgammaR9IDDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG18HEFilter" ),
    isoTag = cms.InputTag( "hltActivityR9ID" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( False ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.8 ),
    thrRegularEE = cms.double( 0.8 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltL1sL1SingleEG15 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG15" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPrePhoton30CaloIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1SingleEG15 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG15" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG30EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG15" ),
    etcutEB = cms.double( 30.0 ),
    etcutEE = cms.double( 30.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG30CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG30EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG30CaloIdVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG30CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton30CaloIdVLIsoL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton30CaloIdVLIsoLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG30CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.5 ),
    thrRegularEE = cms.double( 5.5 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton30CaloIdVLIsoLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton30CaloIdVLIsoLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton30CaloIdVLIsoLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton30CaloIdVLIsoLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton30CaloIdVLIsoLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton30CaloIdVLIsoLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1sL1SingleEG20 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG20" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPrePhoton32CaloIdLPhoton26CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1SingleEG20 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG20" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG32EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 32.0 ),
    etcutEE = cms.double( 32.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG32CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG32CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32CaloIdLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG26EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 26.0 ),
    etcutEE = cms.double( 26.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton32CaloIdLPhoton26CaloIdLEgammaLHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG26EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton32CaloIdLPhoton26CaloIdLEgammaClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton32CaloIdLPhoton26CaloIdLEgammaLHEDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton36IsoVLPhoton22 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG36EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 36.0 ),
    etcutEE = cms.double( 36.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG36HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG36EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton36IsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG36HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton36IsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton36IsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton36IsoVLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton36IsoVLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG22EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 22.0 ),
    etcutEE = cms.double( 22.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoubleIsoEG22HELastFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG22EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton36CaloIdLPhoton22CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG36CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG36EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG36CaloIdLClusterShapeLastFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG36CaloIdLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton36CaloIdLPhoton22CaloIdLEgammaLHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG22EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton36CaloIdLPhoton22CaloIdLEgammaClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton36CaloIdLPhoton22CaloIdLEgammaLHEDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton36CaloIdLIsoVLPhoton22 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG36CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG36CaloIdLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG36CaloIdLIsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG36CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG36CaloIdLIsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG36CaloIdLIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG36CaloIdLIsoVLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG36CaloIdLIsoVLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton40CaloIdLPhoton28CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG40EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 40.0 ),
    etcutEE = cms.double( 40.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG40CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG40EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG40CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG40CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG28EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 28.0 ),
    etcutEE = cms.double( 28.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton40CaloIdLPhoton28CaloIdLEgammaLHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG28EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPhoton40CaloIdLPhoton28CaloIdLEgammaClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton40CaloIdLPhoton28CaloIdLEgammaLHEDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPrePhoton50CaloIdVLIsoL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG50EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 50.0 ),
    etcutEE = cms.double( 50.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG50CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG50EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton50CaloIdVLIsoLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG50CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.5 ),
    thrRegularEE = cms.double( 5.5 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton50CaloIdVLIsoLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton50CaloIdVLIsoLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton50CaloIdVLIsoLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton50CaloIdVLIsoLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton50CaloIdVLIsoLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton50CaloIdVLIsoLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton70CaloIdLHT300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG70EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 70.0 ),
    etcutEE = cms.double( 70.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG70CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG70EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG70CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG70CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton70CaloIdLHT350 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhoton70CaloIdLMHT50 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltMHT50 = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 50.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 30.0, 30.0 ),
    etaJet = cms.vdouble( 3.0, 3.0 )
)
hltPrePhoton70CaloIdLMHT70 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhoton75CaloIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG75EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 75.0 ),
    etcutEE = cms.double( 75.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG75CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG75EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton75CaloIdVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG75CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton75CaloIdVLIsoL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton75CaloIdVLIsoLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG75CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.5 ),
    thrRegularEE = cms.double( 5.5 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton75CaloIdVLIsoLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton75CaloIdVLIsoLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton75CaloIdVLIsoLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton75CaloIdVLIsoLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton75CaloIdVLIsoLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton75CaloIdVLIsoLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton90CaloIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG90EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 90.0 ),
    etcutEE = cms.double( 90.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG90CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG90EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton90CaloIdVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG90CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton90CaloIdVLIsoL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPhoton90CaloIdVLIsoLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltEG90CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.5 ),
    thrRegularEE = cms.double( 5.5 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton90CaloIdVLIsoLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton90CaloIdVLIsoLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton90CaloIdVLIsoLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton90CaloIdVLIsoLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton90CaloIdVLIsoLTrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton90CaloIdVLIsoLHcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.5 ),
    thrRegularEE = cms.double( 3.5 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton125 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG125EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 125.0 ),
    etcutEE = cms.double( 125.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton125HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG125EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton200NoHE = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG200EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 200.0 ),
    etcutEE = cms.double( 200.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreDoublePhoton33 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG33EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 33.0 ),
    etcutEE = cms.double( 33.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG33HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG33EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG33EtFilterUnseededTight = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 33.0 ),
    etcutEE = cms.double( 33.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoublePhoton33EgammaLHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG33EtFilterUnseededTight" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreDoublePhoton33HEVT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG33HEVTFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG33EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleEG33EtDoubleFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 33.0 ),
    etcutEE = cms.double( 33.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoubleEG33HEVTDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG33EtDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreDoublePhoton50 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG50HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG50EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleEG50EtDoubleFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 50.0 ),
    etcutEE = cms.double( 50.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoubleEG50HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG50EtDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreDoublePhoton60 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG60EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 60.0 ),
    etcutEE = cms.double( 60.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG60HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG60EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleEG60EtDoubleFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 60.0 ),
    etcutEE = cms.double( 60.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoubleEG60HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG60EtDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltL1sL1DoubleEG2FwdVeto = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleEG2_FwdVeto" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoublePhoton5IsoVLCEP = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1DoubleEG2FwdVeto = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1DoubleEG2FwdVeto" ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltDoublePhoton5IsoVLEtPhiFilter = cms.EDFilter( "HLTEgammaDoubleEtDeltaPhiFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1DoubleEG2FwdVeto" ),
    etcut = cms.double( 5.0 ),
    minDeltaPhi = cms.double( 2.5 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoublePhoton5IsoVLEgammaHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoublePhoton5IsoVLEtPhiFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.15 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoublePhoton5IsoVLEgammaEcalIsolFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltDoublePhoton5IsoVLEgammaHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 6.0 ),
    thrRegularEE = cms.double( 6.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoublePhoton5IsoVLEgammaHcalIsolFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltDoublePhoton5IsoVLEgammaEcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoublePhoton5IsoVLEgammaTrackIsolFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltDoublePhoton5IsoVLEgammaHcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 4.0 ),
    thrRegularEE = cms.double( 4.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltTowerMakerForHcal = cms.EDProducer( "CaloTowersCreator",
    EBThreshold = cms.double( 0.07 ),
    EEThreshold = cms.double( 0.3 ),
    UseEtEBTreshold = cms.bool( False ),
    UseEtEETreshold = cms.bool( False ),
    UseSymEBTreshold = cms.bool( False ),
    UseSymEETreshold = cms.bool( False ),
    HcalThreshold = cms.double( -1000.0 ),
    HBThreshold = cms.double( 0.7 ),
    HESThreshold = cms.double( 0.8 ),
    HEDThreshold = cms.double( 0.8 ),
    HOThreshold0 = cms.double( 3.5 ),
    HOThresholdPlus1 = cms.double( 3.5 ),
    HOThresholdMinus1 = cms.double( 3.5 ),
    HOThresholdPlus2 = cms.double( 3.5 ),
    HOThresholdMinus2 = cms.double( 3.5 ),
    HF1Threshold = cms.double( 0.5 ),
    HF2Threshold = cms.double( 0.85 ),
    EBWeight = cms.double( 1.0E-99 ),
    EEWeight = cms.double( 1.0E-99 ),
    HBWeight = cms.double( 1.0 ),
    HESWeight = cms.double( 1.0 ),
    HEDWeight = cms.double( 1.0 ),
    HOWeight = cms.double( 1.0E-99 ),
    HF1Weight = cms.double( 1.0 ),
    HF2Weight = cms.double( 1.0 ),
    EcutTower = cms.double( -1000.0 ),
    EBSumThreshold = cms.double( 0.2 ),
    EESumThreshold = cms.double( 0.45 ),
    UseHO = cms.bool( False ),
    MomConstrMethod = cms.int32( 1 ),
    MomHBDepth = cms.double( 0.2 ),
    MomHEDepth = cms.double( 0.4 ),
    MomEBDepth = cms.double( 0.3 ),
    MomEEDepth = cms.double( 0.0 ),
    hbheInput = cms.InputTag( "hltHbhereco" ),
    hoInput = cms.InputTag( "hltHoreco" ),
    hfInput = cms.InputTag( "hltHfreco" ),
    AllowMissingInputs = cms.bool( True ),
    HcalAcceptSeverityLevel = cms.uint32( 9 ),
    EcalAcceptSeverityLevel = cms.uint32( 3 ),
    UseHcalRecoveredHits = cms.bool( False ),
    UseEcalRecoveredHits = cms.bool( False ),
    UseRejectedHitsOnly = cms.bool( False ),
    HcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    EcalAcceptSeverityLevelForRejectedHit = cms.uint32( 9999 ),
    UseRejectedRecoveredHcalHits = cms.bool( False ),
    UseRejectedRecoveredEcalHits = cms.bool( False ),
    EBGrid = cms.vdouble(  ),
    EBWeights = cms.vdouble(  ),
    EEGrid = cms.vdouble(  ),
    EEWeights = cms.vdouble(  ),
    HBGrid = cms.vdouble(  ),
    HBWeights = cms.vdouble(  ),
    HESGrid = cms.vdouble(  ),
    HESWeights = cms.vdouble(  ),
    HEDGrid = cms.vdouble(  ),
    HEDWeights = cms.vdouble(  ),
    HOGrid = cms.vdouble(  ),
    HOWeights = cms.vdouble(  ),
    HF1Grid = cms.vdouble(  ),
    HF1Weights = cms.vdouble(  ),
    HF2Grid = cms.vdouble(  ),
    HF2Weights = cms.vdouble(  ),
    ecalInputs = cms.VInputTag(  )
)
hltHcalTowerFilter = cms.EDFilter( "HLTHcalTowerFilter",
    inputTag = cms.InputTag( "hltTowerMakerForHcal" ),
    saveTags = cms.bool( False ),
    MinE_HB = cms.double( 1.5 ),
    MinE_HE = cms.double( 2.5 ),
    MinE_HF = cms.double( 9.0 ),
    MaxN_HB = cms.int32( 2 ),
    MaxN_HE = cms.int32( 2 ),
    MaxN_HF = cms.int32( 8 )
)
hltL1sL1SingleEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1SingleEG5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreL1SingleEG12 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreEle8 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1SingleEG5 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG5" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG8EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG5" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1IsoStartUpElectronPixelSeeds = cms.EDProducer( "ElectronSeedProducer",
    barrelSuperClusters = cms.InputTag( "hltCorrectedHybridSuperClustersL1Isolated" ),
    endcapSuperClusters = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1Isolated" ),
    SeedConfiguration = cms.PSet( 
      searchInTIDTEC = cms.bool( True ),
      HighPtThreshold = cms.double( 35.0 ),
      r2MinF = cms.double( -0.15 ),
      OrderedHitsFactoryPSet = cms.PSet( 
        maxElement = cms.uint32( 0 ),
        ComponentName = cms.string( "StandardHitPairGenerator" ),
        SeedingLayers = cms.string( "hltESPMixedLayerPairs" ),
        useOnDemandTracker = cms.untracked.int32( 0 )
      ),
      DeltaPhi1Low = cms.double( 0.23 ),
      DeltaPhi1High = cms.double( 0.08 ),
      ePhiMin1 = cms.double( -0.08 ),
      PhiMin2 = cms.double( -0.0040 ),
      LowPtThreshold = cms.double( 3.0 ),
      RegionPSet = cms.PSet( 
        deltaPhiRegion = cms.double( 0.4 ),
        originHalfLength = cms.double( 15.0 ),
        useZInVertex = cms.bool( True ),
        deltaEtaRegion = cms.double( 0.1 ),
        ptMin = cms.double( 1.5 ),
        originRadius = cms.double( 0.2 ),
        VertexProducer = cms.InputTag( "dummyVertices" )
      ),
      maxHOverE = cms.double( 999999.0 ),
      dynamicPhiRoad = cms.bool( False ),
      ePhiMax1 = cms.double( 0.04 ),
      DeltaPhi2 = cms.double( 0.0040 ),
      measurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
      SizeWindowENeg = cms.double( 0.675 ),
      nSigmasDeltaZ1 = cms.double( 5.0 ),
      rMaxI = cms.double( 0.2 ),
      PhiMax2 = cms.double( 0.0040 ),
      preFilteredSeeds = cms.bool( False ),
      r2MaxF = cms.double( 0.15 ),
      pPhiMin1 = cms.double( -0.04 ),
      initialSeeds = cms.InputTag( "globalPixelSeeds:GlobalPixel" ),
      pPhiMax1 = cms.double( 0.08 ),
      hbheModule = cms.string( "hbhereco" ),
      SCEtCut = cms.double( 3.0 ),
      z2MaxB = cms.double( 0.09 ),
      fromTrackerSeeds = cms.bool( True ),
      hcalRecHits = cms.InputTag( "hltHbhereco" ),
      z2MinB = cms.double( -0.09 ),
      hbheInstance = cms.string( "" ),
      rMinI = cms.double( -0.2 ),
      hOverEConeSize = cms.double( 0.0 ),
      hOverEHBMinE = cms.double( 999999.0 ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      applyHOverECut = cms.bool( False ),
      hOverEHFMinE = cms.double( 999999.0 )
    )
)
hltL1NonIsoStartUpElectronPixelSeeds = cms.EDProducer( "ElectronSeedProducer",
    barrelSuperClusters = cms.InputTag( "hltCorrectedHybridSuperClustersL1NonIsolated" ),
    endcapSuperClusters = cms.InputTag( "hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated" ),
    SeedConfiguration = cms.PSet( 
      searchInTIDTEC = cms.bool( True ),
      HighPtThreshold = cms.double( 35.0 ),
      r2MinF = cms.double( -0.15 ),
      OrderedHitsFactoryPSet = cms.PSet( 
        maxElement = cms.uint32( 0 ),
        ComponentName = cms.string( "StandardHitPairGenerator" ),
        SeedingLayers = cms.string( "hltESPMixedLayerPairs" ),
        useOnDemandTracker = cms.untracked.int32( 0 )
      ),
      DeltaPhi1Low = cms.double( 0.23 ),
      DeltaPhi1High = cms.double( 0.08 ),
      ePhiMin1 = cms.double( -0.08 ),
      PhiMin2 = cms.double( -0.0040 ),
      LowPtThreshold = cms.double( 3.0 ),
      RegionPSet = cms.PSet( 
        deltaPhiRegion = cms.double( 0.4 ),
        originHalfLength = cms.double( 15.0 ),
        useZInVertex = cms.bool( True ),
        deltaEtaRegion = cms.double( 0.1 ),
        ptMin = cms.double( 1.5 ),
        originRadius = cms.double( 0.2 ),
        VertexProducer = cms.InputTag( "dummyVertices" )
      ),
      maxHOverE = cms.double( 999999.0 ),
      dynamicPhiRoad = cms.bool( False ),
      ePhiMax1 = cms.double( 0.04 ),
      DeltaPhi2 = cms.double( 0.0040 ),
      measurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
      SizeWindowENeg = cms.double( 0.675 ),
      nSigmasDeltaZ1 = cms.double( 5.0 ),
      rMaxI = cms.double( 0.2 ),
      PhiMax2 = cms.double( 0.0040 ),
      preFilteredSeeds = cms.bool( False ),
      r2MaxF = cms.double( 0.15 ),
      pPhiMin1 = cms.double( -0.04 ),
      initialSeeds = cms.InputTag( "globalPixelSeeds:GlobalPixel" ),
      pPhiMax1 = cms.double( 0.08 ),
      hbheModule = cms.string( "hbhereco" ),
      SCEtCut = cms.double( 3.0 ),
      z2MaxB = cms.double( 0.09 ),
      fromTrackerSeeds = cms.bool( True ),
      hcalRecHits = cms.InputTag( "hltHbhereco" ),
      z2MinB = cms.double( -0.09 ),
      hbheInstance = cms.string( "" ),
      rMinI = cms.double( -0.2 ),
      hOverEConeSize = cms.double( 0.0 ),
      hOverEHBMinE = cms.double( 999999.0 ),
      beamSpot = cms.InputTag( "offlineBeamSpot" ),
      applyHOverECut = cms.bool( False ),
      hOverEHFMinE = cms.double( 999999.0 )
    )
)
hltEle8PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle8HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreEle8CaloIdLCaloIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle8CaloIdLCaloIsoVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLCaloIsoVLEcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLCaloIsoVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLEcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLCaloIsoVLHcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLCaloIsoVLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLHcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreEle8CaloIdLTrkIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle8CaloIdLTrkIdVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLTrkIdVLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPixelMatchElectronsL1Iso = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1IsoWithMaterialTracks" ),
    BSProducer = cms.InputTag( "offlineBeamSpot" )
)
hltPixelMatchElectronsL1NonIso = cms.EDProducer( "EgammaHLTPixelMatchElectronProducers",
    TrackProducer = cms.InputTag( "hltCtfL1NonIsoWithMaterialTracks" ),
    BSProducer = cms.InputTag( "offlineBeamSpot" )
)
hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltElectronL1IsoDetaDphi = cms.EDProducer( "EgammaHLTElectronDetaDphiProducer",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    BSProducer = cms.InputTag( "offlineBeamSpot" )
)
hltElectronL1NonIsoDetaDphi = cms.EDProducer( "EgammaHLTElectronDetaDphiProducer",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    BSProducer = cms.InputTag( "offlineBeamSpot" )
)
hltEle8CaloIdLTrkIdVLDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle8CaloIdLTrkIdVLDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG15EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG12" ),
    etcutEB = cms.double( 15.0 ),
    etcutEE = cms.double( 15.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdTCaloIsoTEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTCaloIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdTCaloIsoTEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTCaloIsoTHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTCaloIsoTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTCaloIsoTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle15CaloIdVTCaloIsoTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle15CaloIdVTCaloIsoTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1IsoElectronTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    trackProducer = cms.InputTag( "hltL1IsoEgammaRegionalCTFFinalFitWithMaterial" ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.3 ),
    egTrkIsoZSpan = cms.double( 0.15 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.03 ),
    egTrkIsoStripBarrel = cms.double( 0.03 ),
    egTrkIsoStripEndcap = cms.double( 0.03 )
)
hltL1NonIsoElectronTrackIsol = cms.EDProducer( "EgammaHLTElectronTrackIsolationProducers",
    electronProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    trackProducer = cms.InputTag( "hltL1NonIsoEgammaRegionalCTFFinalFitWithMaterial" ),
    egTrkIsoPtMin = cms.double( 1.0 ),
    egTrkIsoConeSize = cms.double( 0.3 ),
    egTrkIsoZSpan = cms.double( 0.15 ),
    egTrkIsoRSpan = cms.double( 999999.0 ),
    egTrkIsoVetoConeSize = cms.double( 0.03 ),
    egTrkIsoStripBarrel = cms.double( 0.03 ),
    egTrkIsoStripEndcap = cms.double( 0.03 )
)
hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle17CaloIdLCaloIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG17EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG12" ),
    etcutEB = cms.double( 17.0 ),
    etcutEE = cms.double( 17.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG17CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG17CaloIdLCaloIsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG17CaloIdLCaloIsoVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17CaloIdLCaloIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG17CaloIdLCaloIsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17CaloIdLCaloIsoVLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdLCaloIsoVLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG17CaloIdLCaloIsoVLHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleEG8EtFilterUnseeded = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdIsoEle8CaloIdIsoClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdIsoEle8CaloIdIsoEcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdIsoEle8CaloIdIsoClusterShapeDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdIsoEle8CaloIdIsoHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdIsoEle8CaloIdIsoEcalIsolDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdIsoEle8CaloIdIsoHcalIsolDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdIsoEle8CaloIdIsoHEDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdIsoEle8CaloIdIsoHcalIsolDoubleFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8EcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8EcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.05 ),
    thrOverPtEE = cms.double( 0.05 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter = cms.EDFilter( "HLTPMMassFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEDoubleFilter" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    lowerMassCut = cms.double( 30.0 ),
    upperMassCut = cms.double( 999999.0 ),
    nZcandcut = cms.int32( 1 ),
    isElectron1 = cms.untracked.bool( False ),
    isElectron2 = cms.untracked.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8L1MatchFilterRegional = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG12" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8L1MatchFilterRegional" ),
    etcutEB = cms.double( 17.0 ),
    etcutEE = cms.double( 17.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.05 ),
    thrOverPtEE = cms.double( 0.05 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DoubleEtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DoubleEtFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchDoubleFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEDoubleFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter = cms.EDFilter( "HLTPMMassFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchDoubleFilter" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    lowerMassCut = cms.double( 30.0 ),
    upperMassCut = cms.double( 999999.0 ),
    nZcandcut = cms.int32( 1 ),
    isElectron1 = cms.untracked.bool( False ),
    isElectron2 = cms.untracked.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreEle17CaloIdLCaloIsoVLEle15HFL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHFEMClusters = cms.EDProducer( "HFEMClusterProducer",
    hits = cms.InputTag( "hltHfreco" ),
    minTowerEnergy = cms.double( 4.0 ),
    seedThresholdET = cms.double( 5.0 ),
    maximumSL = cms.double( 0.98 ),
    maximumRenergy = cms.double( 50.0 ),
    usePMTFlag = cms.bool( True ),
    usePulseFlag = cms.bool( True ),
    forcePulseFlagMC = cms.bool( False ),
    correctionType = cms.int32( 1 )
)
hltHFRecoEcalCandidate = cms.EDProducer( "HFRecoEcalCandidateProducer",
    hfclusters = cms.InputTag( "hltHFEMClusters" ),
    Correct = cms.bool( True ),
    e9e25Cut = cms.double( 0.9 ),
    intercept2DCut = cms.double( 0.2 ),
    e1e9Cut = cms.vdouble( -1.0, 99.0 ),
    eCOREe9Cut = cms.vdouble( -1.0, 99.0 ),
    eSeLCut = cms.vdouble( -1.0, 99.9 )
)
hltHFEMFilter = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltHFRecoEcalCandidate" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG27EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG15" ),
    etcutEB = cms.double( 27.0 ),
    etcutEE = cms.double( 27.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG27CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG27EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVLCaloIsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG27CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVLCaloIsoVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVLCaloIsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVLCaloIsoVLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVLCaloIsoVLOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle27CaloIdVLCaloIsoVLTrkIdVLDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle27CaloIdVLCaloIsoVLTrkIdVLDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLTrkIdVLDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVLCaloIsoVLTrkIdVLDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.2 ),
    thrOverPtEE = cms.double( 0.2 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle27CaloIdVTCaloIsoTTrkIdTTrkIsoT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG27CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG27EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdTCaloIsoTEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG27CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVTCaloIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdTCaloIsoTEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVTCaloIsoTHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 999999.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVTCaloIsoTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle27CaloIdVTCaloIsoTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle27CaloIdVTCaloIsoTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle27CaloIdVTCaloIsoTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle27CaloIdVTCaloIsoTTrkIdTDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG32CaloIdVLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.024 ),
    thrRegularEE = cms.double( 0.04 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVLCaloIsoVLEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32CaloIdVLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVLCaloIsoVLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVLCaloIsoVLHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVLCaloIsoVLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVLCaloIsoVLOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle32CaloIdVLCaloIsoVLTrkIdVLDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdVLCaloIsoVLTrkIdVLDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLTrkIdVLDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVLCaloIsoVLTrkIdVLDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.2 ),
    thrOverPtEE = cms.double( 0.2 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle32CaloIdVTCaloIsoTTrkIdTTrkIsoT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG32CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTEcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG32CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVTCaloIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTEcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVTCaloIsoTHcalIsoFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 999999.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVTCaloIsoTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTHcalIsoFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdVTCaloIsoTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle32CaloIdVTCaloIsoTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdVTCaloIsoTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdVTCaloIsoTTrkIdTDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17L1MatchFilterRegional = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1SingleEG20" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17L1MatchFilterRegional" ),
    etcutEB = cms.double( 32.0 ),
    etcutEE = cms.double( 32.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.1 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DoubleEtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 17.0 ),
    etcutEE = cms.double( 17.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DoubleEtFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreEle45CaloIdVTTrkIdT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG45EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 45.0 ),
    etcutEE = cms.double( 45.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG45CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG45EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG45CaloIdVTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG45CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle45CaloIdVTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG45CaloIdVTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle45CaloIdVTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle45CaloIdVTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle45CaloIdVTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle45CaloIdVTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle45CaloIdVTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle45CaloIdVTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle52CaloIdVTTrkIdT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG52EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG20" ),
    etcutEB = cms.double( 52.0 ),
    etcutEE = cms.double( 52.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG52CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG52EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG52CaloIdVTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG52CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle52CaloIdVTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG52CaloIdVTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle52CaloIdVTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle52CaloIdVTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle52CaloIdVTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle52CaloIdVTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle52CaloIdVTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle52CaloIdVTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1sL1DoubleEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleEG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleEle8CaloIdLTrkIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1DoubleEG5 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1DoubleEG5" ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltDoubleEG8EtFilterL1DoubleEG5 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1DoubleEG5" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterL1DoubleEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8DetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreDoubleEle33 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle33PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG33HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleEG33HEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG33EtDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDiEle33PixelMatchDoubleFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltDoubleEG33HEDoubleFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreDoubleEle33CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG33CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG33HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle33CaloIdLPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG33CaloIdLClusterShapeFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleEG33CaloIdLClusterShapeDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG33HEDoubleFilter" ),
    isoTag = cms.InputTag( "hltActivityPhotonClusterShape" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDiEle33CaloIdLPixelMatchDoubleFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltDoubleEG33CaloIdLClusterShapeDoubleFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltL1sSingleIsoTau35Trk20MET45 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleTauJet52 OR L1_SingleJet68" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreSingleIsoTau35Trk20MET45 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltCaloTowersTau1Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Tau' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 0 )
)
hltIconeTau1Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersTau1Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersTau2Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Tau' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 1 )
)
hltIconeTau2Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersTau2Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersTau3Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Tau' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 2 )
)
hltIconeTau3Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersTau3Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersTau4Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Tau' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 3 )
)
hltIconeTau4Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersTau4Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersCentral1Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Central' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 0 )
)
hltIconeCentral1Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersCentral1Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersCentral2Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Central' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 1 )
)
hltIconeCentral2Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersCentral2Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersCentral3Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Central' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 2 )
)
hltIconeCentral3Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersCentral3Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltCaloTowersCentral4Regional = cms.EDProducer( "CaloTowerCreatorForTauHLT",
    towers = cms.InputTag( "hltTowerMakerForJets" ),
    UseTowersInCone = cms.double( 0.8 ),
    TauTrigger = cms.InputTag( 'l1extraParticles','Central' ),
    minimumEt = cms.double( 0.5 ),
    minimumE = cms.double( 0.8 ),
    TauId = cms.int32( 3 )
)
hltIconeCentral4Regional = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.2 ),
    src = cms.InputTag( "hltCaloTowersCentral4Regional" ),
    srcPVs = cms.InputTag( "offlinePrimaryVertices" ),
    jetType = cms.string( "CaloJet" ),
    jetPtMin = cms.double( 1.0 ),
    inputEtMin = cms.double( 0.3 ),
    inputEMin = cms.double( 0.0 ),
    doPVCorrection = cms.bool( False ),
    doPUOffsetCorr = cms.bool( False ),
    nSigmaPU = cms.double( 1.0 ),
    radiusPU = cms.double( 0.5 ),
    Active_Area_Repeats = cms.int32( 5 ),
    GhostArea = cms.double( 0.01 ),
    Ghost_EtaMax = cms.double( 6.0 ),
    maxBadEcalCells = cms.uint32( 9999999 ),
    maxRecoveredEcalCells = cms.uint32( 9999999 ),
    maxProblematicEcalCells = cms.uint32( 9999999 ),
    maxBadHcalCells = cms.uint32( 9999999 ),
    maxRecoveredHcalCells = cms.uint32( 9999999 ),
    maxProblematicHcalCells = cms.uint32( 9999999 ),
    doAreaFastjet = cms.bool( False ),
    doRhoFastjet = cms.bool( False ),
    subtractorName = cms.string( "" ),
    sumRecHits = cms.bool( False )
)
hltL2TauJets = cms.EDProducer( "L2TauJetsMerger",
    EtMin = cms.double( 20.0 ),
    JetSrc = cms.VInputTag( 'hltIconeTau1Regional','hltIconeTau2Regional','hltIconeTau3Regional','hltIconeTau4Regional','hltIconeCentral1Regional','hltIconeCentral2Regional','hltIconeCentral3Regional','hltIconeCentral4Regional' )
)
hltFilterL2EtCutSingleIsoPFTau35Trk20MET45 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltL2TauJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 35.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
hltPFTauTightIso35 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIso" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 35.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTauTightIso35Track = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIsoTrackFinding" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 35.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTauTightIsoTrackPt20Discriminator = cms.EDProducer( "PFRecoTauDiscriminationByLeadingObjectPtCut",
    PFTauProducer = cms.InputTag( "hltPFTausTightIso" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
    UseOnlyChargedHadrons = cms.bool( True ),
    MinPtLeadingObject = cms.double( 20.0 )
)
hltSelectedPFTauTightIsoTrackPt20 = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackPt20Discriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltConvPFTauTightIsoTrackPt20 = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTauTightIsoTrackPt20" )
)
hltFilterSingleIsoPFTau35Trk20LeadTrackPt20 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTauTightIsoTrackPt20" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 35.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltSelectedPFTauTightIsoTrackPt20Isolation = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTausTightIso" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoTrackPt20Discriminator" ),
        selectionCut = cms.double( 0.5 )
      ),
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTightIsoIsolationDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltConvPFTauTightIsoTrackPt20Isolation = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTauTightIsoTrackPt20Isolation" )
)
hltL1HLTSingleIsoPFTau35Trk20Met45JetsMatch = cms.EDProducer( "L1HLTJetsMatching",
    JetSrc = cms.InputTag( "hltConvPFTauTightIsoTrackPt20Isolation" ),
    L1TauTrigger = cms.InputTag( "hltL1sSingleIsoTau35Trk20MET45" ),
    EtMin = cms.double( 0.0 )
)
hltFilterSingleIsoPFTau35Trk20MET45LeadTrack20MET45IsolationL1HLTMatched = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltL1HLTSingleIsoPFTau35Trk20Met45JetsMatch" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 35.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltL1sDoubleIsoTau25Trk5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleTauJet28 OR L1_DoubleJet52" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleIsoTau25Trk5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltFilterL2EtCutDoublePFIsoTau25Trk5 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltL2TauJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.1 ),
    MinN = cms.int32( 2 )
)
hltDoublePFTauTightIso25Track = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIsoTrackFinding" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.1 ),
    MinN = cms.int32( 2 )
)
hltDoublePFTauTightIso25Track5 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTightIsoTrackPt5" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.1 ),
    MinN = cms.int32( 2 )
)
hltL1HLTDoubleIsoPFTau25Trk5JetsMatch = cms.EDProducer( "L1HLTJetsMatching",
    JetSrc = cms.InputTag( "hltConvPFTausTightIsoTrackPt5Isolation" ),
    L1TauTrigger = cms.InputTag( "hltL1sDoubleIsoTau25Trk5" ),
    EtMin = cms.double( 0.0 )
)
hltFilterDoubleIsoPFTau25Trk5LeadTrack5IsolationL1HLTMatched = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltL1HLTDoubleIsoPFTau25Trk5JetsMatch" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.1 ),
    MinN = cms.int32( 2 )
)
hltL1sL1Mu3Jet16Central = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_Mu3_Jet16_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreBTagMuDiJet20Mu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBDiJet20Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltBSoftMuonGetJetsFromDiJet20 = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBDiJet20Central" )
)
hltSelector4JetsDiJet20 = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltBSoftMuonGetJetsFromDiJet20" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBSoftMuonDiJet20L25Jets = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelector4JetsDiJet20" ),
    filter = cms.bool( False ),
    etMin = cms.double( 20.0 )
)
hltBSoftMuonDiJet20L25TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet20L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltL2Muons" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet20L25BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet20L25TagInfos' )
)
hltBSoftMuonDiJet20L25FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet20L25BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBSoftMuonMu5L3 = cms.EDFilter( "RecoTrackRefSelector",
    src = cms.InputTag( "hltL3Muons" ),
    maxChi2 = cms.double( 10000.0 ),
    tip = cms.double( 120.0 ),
    minRapidity = cms.double( -5.0 ),
    lip = cms.double( 300.0 ),
    ptMin = cms.double( 5.0 ),
    maxRapidity = cms.double( 5.0 ),
    quality = cms.vstring(  ),
    algorithm = cms.vstring(  ),
    minHit = cms.int32( 0 ),
    min3DHit = cms.int32( 0 ),
    beamSpot = cms.InputTag( "offlineBeamSpot" )
)
hltBSoftMuonDiJet20Mu5SelL3TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet20L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltBSoftMuonMu5L3" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet20Mu5SelL3TagInfos' )
)
hltBSoftMuonDiJet20Mu5L3FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltL1sL1Mu3Jet20Central = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_Mu3_Jet20_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreBTagMuDiJet40Mu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBDiJet40Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltBSoftMuonGetJetsFromDiJet40 = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBDiJet40Central" )
)
hltSelector4JetsDiJet40 = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltBSoftMuonGetJetsFromDiJet40" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBSoftMuonDiJet40L25Jets = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelector4JetsDiJet40" ),
    filter = cms.bool( False ),
    etMin = cms.double( 40.0 )
)
hltBSoftMuonDiJet40L25TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet40L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltL2Muons" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet40L25BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet40L25TagInfos' )
)
hltBSoftMuonDiJet40L25FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet40L25BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBSoftMuonDiJet40Mu5SelL3TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet40L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltBSoftMuonMu5L3" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet40Mu5SelL3TagInfos' )
)
hltBSoftMuonDiJet40Mu5L3FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltL1sL1Mu3Jet28Central = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_Mu3_Jet28_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreBTagMuDiJet70Mu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBDiJet70Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 70.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltBSoftMuonGetJetsFromDiJet70 = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBDiJet70Central" )
)
hltSelector4JetsDiJet70 = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltBSoftMuonGetJetsFromDiJet70" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBSoftMuonDiJet70L25Jets = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelector4JetsDiJet70" ),
    filter = cms.bool( False ),
    etMin = cms.double( 70.0 )
)
hltBSoftMuonDiJet70L25TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet70L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltL2Muons" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet70L25BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet70L25TagInfos' )
)
hltBSoftMuonDiJet70L25FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet70L25BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBSoftMuonDiJet70Mu5SelL3TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet70L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltBSoftMuonMu5L3" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet70Mu5SelL3TagInfos' )
)
hltBSoftMuonDiJet70Mu5L3FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltPreBTagMuDiJet110Mu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBDiJet110Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 110.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltBSoftMuonGetJetsFromDiJet110 = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBDiJet110Central" )
)
hltSelector4JetsDiJet110 = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltBSoftMuonGetJetsFromDiJet110" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBSoftMuonDiJet110L25Jets = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelector4JetsDiJet110" ),
    filter = cms.bool( False ),
    etMin = cms.double( 110.0 )
)
hltBSoftMuonDiJet110L25TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet110L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltL2Muons" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet110L25BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet110L25TagInfos' )
)
hltBSoftMuonDiJet110L25FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet110L25BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBSoftMuonDiJet110Mu5SelL3TagInfos = cms.EDProducer( "SoftLepton",
    jets = cms.InputTag( "hltBSoftMuonDiJet110L25Jets" ),
    primaryVertex = cms.InputTag( "nominal" ),
    leptons = cms.InputTag( "hltBSoftMuonMu5L3" ),
    leptonCands = cms.InputTag( "" ),
    leptonId = cms.InputTag( "" ),
    refineJetAxis = cms.uint32( 0 ),
    leptonDeltaRCut = cms.double( 0.4 ),
    leptonChi2Cut = cms.double( 0.0 ),
    muonSelection = cms.uint32( 0 )
)
hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPSoftLeptonByDistance" ),
    tagInfos = cms.VInputTag( 'hltBSoftMuonDiJet110Mu5SelL3TagInfos' )
)
hltBSoftMuonDiJet110Mu5L3FilterByDR = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR" ),
    MinTag = cms.double( 0.5 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltL1sL1MuOpenEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_MuOpen_EG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreMu5Ele8CaloIdLTrkIdVLEle8 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenEG5L1Filtered3 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 3.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2MuOpenEG5L2Filtered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L1Filtered3" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltMuOpenEG5L3Filtered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2MuOpenEG5L2Filtered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltEGRegionalL1MuOpenEG5 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG8EtFilterMuOpenEG5 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1MuOpenEG5" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG8CaloIdLClusterShapeFilterMuOpenEG5 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG8CaloIdLHEFilterMuOpenEG5 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8CaloIdLClusterShapeFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG8CaloIdLPixelMatchFilterMuOpenEG5 = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG8CaloIdLHEFilterMuOpenEG5" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilterRegionalMuOpenEG5 = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEG8CaloIdLPixelMatchFilterMuOpenEG5" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle8CaloIdLTrkIdVLDetaFilterMuOpenEG5 = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilterRegionalMuOpenEG5" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle8CaloIdLTrkIdVLDphiFilterMuOpenEG5 = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdLTrkIdVLDetaFilterMuOpenEG5" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLTNonIsoMu5DoubleEle8NoCandHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterUnseeded" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoubleEG8NoCandPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoMu5DoubleEle8NoCandHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltActivityStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( True ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltPreMu8R005MR200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1Mu3Jet20L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1Mu3Jet20Central" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL2Mu8Jet20L2Filtered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu3Jet20L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltL3Mu8Jet20L3Filtered8 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL2Mu8Jet20L2Filtered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 8.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltR005MR200 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.05 ),
    minMR = cms.double( 200.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreMu8R020MR200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR200 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 200.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreMu8R025MR200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR025MR200 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.25 ),
    minMR = cms.double( 200.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreHT250Mu5pfMHT35 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sL1SingleMuOpenCandidate = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( False ),
    L1NrBxInEvent = cms.int32( 1 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMuOpen" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltL1HTT100L1MuFiltered3 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpenCandidate" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 3.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1HTT100singleMuL2PreFiltered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1HTT100L1MuFiltered3" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1HTT100singleMuL3PreFiltered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1HTT100singleMuL2PreFiltered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPFMHT35Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 35.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreHT250Mu15pfMHT20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHTT100L1MuFiltered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpenCandidate" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1HTT100singleMuL2PreFiltered10 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHTT100L1MuFiltered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 10.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1HTT100singleMuL3PreFiltered15 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1HTT100singleMuL2PreFiltered10" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 15.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPFMHT20Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 20.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreHT300Mu5pfMHT40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPFMHT40Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 40.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreHT350Mu5pfMHT45 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPFMHT45Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 45.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreMu3DiJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 2 )
)
hltPreMu3TriJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTripleJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 3 )
)
hltPreMu3QuadJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltQuadJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltJetIDPassedCorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 4 )
)
hltPreMu8Ele17CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenEG5L1Filtered5 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 5.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1MuOpenEG5L2Filtered5 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L1Filtered5" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1MuOpenEG5L3Filtered8 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L2Filtered5" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 8.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltEG17EtFilterL1MuOpenEG5 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1MuOpenEG5" ),
    etcutEB = cms.double( 17.0 ),
    etcutEE = cms.double( 17.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLMu8Ele17ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17EtFilterL1MuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoMu8Ele17HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLMu8Ele17ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoMu8Ele17HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreMu8Photon20CaloIdVTIsoT = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG20EtFilterMuOpenEG5 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1MuOpenEG5" ),
    etcutEB = cms.double( 20.0 ),
    etcutEE = cms.double( 20.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTMu8ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTMu8EcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTMu8ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 5.0 ),
    thrRegularEE = cms.double( 5.0 ),
    thrOverEEB = cms.double( 0.012 ),
    thrOverEEE = cms.double( 0.012 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTMu8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTMu8EcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTMu8HcalIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTMu8HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.0 ),
    thrRegularEE = cms.double( 3.0 ),
    thrOverEEB = cms.double( 0.0050 ),
    thrOverEEE = cms.double( 0.0050 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter = cms.EDFilter( "HLTEgammaGenericQuadraticFilter",
    candTag = cms.InputTag( "hltPhoton20CaloIdVTIsoTMu8HcalIsoFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHollowTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHollowTrackIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( 3.0 ),
    thrRegularEE = cms.double( 3.0 ),
    thrOverEEB = cms.double( 0.0020 ),
    thrOverEEE = cms.double( 0.0020 ),
    thrOverE2EB = cms.double( 0.0 ),
    thrOverE2EE = cms.double( 0.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1SingleMuOpenEG5L1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltSingleMu5EG5L2Filtered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1SingleMuOpenEG5L1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltSingleMu8EG5L3Filtered8 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltSingleMu5EG5L2Filtered3" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 8.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu8Jet40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltJet40 = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltPreMu15Photon20CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenEG5L3Filtered15 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L2Filtered5" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 15.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltMu15Photon20CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG20EtFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltMu15Photon20CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltMu15Photon20CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreMu15DoublePhoton15CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDoubleEG15EtFilterL1MuOpenEG5 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1MuOpenEG5" ),
    etcutEB = cms.double( 15.0 ),
    etcutEE = cms.double( 15.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltMu15DiPhoton15CaloIdLClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG15EtFilterL1MuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltMu15DiPhoton15CaloIdLHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltMu15DiPhoton15CaloIdLClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreMu15IsoPFTau15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTauJet5 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltAntiKT5CaloJetsEt5" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 5.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFJet15 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltAntiKT5ConvPFJetsForTaus" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTaus = cms.EDProducer( "PFRecoTauProducer",
    PFTauTagInfoProducer = cms.InputTag( "hltPFTauTagInfo" ),
    ElectronPreIDProducer = cms.InputTag( "elecpreid" ),
    PVProducer = cms.InputTag( "hltPixelVertices" ),
    Algorithm = cms.string( "ConeBased" ),
    smearedPVsigmaX = cms.double( 0.0015 ),
    smearedPVsigmaY = cms.double( 0.0015 ),
    smearedPVsigmaZ = cms.double( 0.0050 ),
    JetPtMin = cms.double( 0.0 ),
    LeadPFCand_minPt = cms.double( 0.0 ),
    UseChargedHadrCandLeadChargedHadrCand_tksDZconstraint = cms.bool( True ),
    ChargedHadrCandLeadChargedHadrCand_tksmaxDZ = cms.double( 0.4 ),
    LeadTrack_minPt = cms.double( 0.0 ),
    UseTrackLeadTrackDZconstraint = cms.bool( False ),
    TrackLeadTrack_maxDZ = cms.double( 0.4 ),
    MatchingConeMetric = cms.string( "DR" ),
    MatchingConeSizeFormula = cms.string( "0.2" ),
    MatchingConeSize_min = cms.double( 0.0 ),
    MatchingConeSize_max = cms.double( 0.6 ),
    TrackerSignalConeMetric = cms.string( "DR" ),
    TrackerSignalConeSizeFormula = cms.string( "0.2" ),
    TrackerSignalConeSize_min = cms.double( 0.0 ),
    TrackerSignalConeSize_max = cms.double( 0.2 ),
    TrackerIsolConeMetric = cms.string( "DR" ),
    TrackerIsolConeSizeFormula = cms.string( "0.4" ),
    TrackerIsolConeSize_min = cms.double( 0.0 ),
    TrackerIsolConeSize_max = cms.double( 0.6 ),
    ECALSignalConeMetric = cms.string( "DR" ),
    ECALSignalConeSizeFormula = cms.string( "0.2" ),
    ECALSignalConeSize_min = cms.double( 0.0 ),
    ECALSignalConeSize_max = cms.double( 0.6 ),
    ECALIsolConeMetric = cms.string( "DR" ),
    ECALIsolConeSizeFormula = cms.string( "0.5" ),
    ECALIsolConeSize_min = cms.double( 0.0 ),
    ECALIsolConeSize_max = cms.double( 0.6 ),
    HCALSignalConeMetric = cms.string( "DR" ),
    HCALSignalConeSizeFormula = cms.string( "0.1" ),
    HCALSignalConeSize_min = cms.double( 0.0 ),
    HCALSignalConeSize_max = cms.double( 0.6 ),
    HCALIsolConeMetric = cms.string( "DR" ),
    HCALIsolConeSizeFormula = cms.string( "0.5" ),
    HCALIsolConeSize_min = cms.double( 0.0 ),
    HCALIsolConeSize_max = cms.double( 0.6 ),
    Rphi = cms.double( 2.0 ),
    MaxEtInEllipse = cms.double( 2.0 ),
    AddEllipseGammas = cms.bool( False ),
    AreaMetric_recoElements_maxabsEta = cms.double( 2.5 ),
    ChargedHadrCand_IsolAnnulus_minNhits = cms.uint32( 0 ),
    Track_IsolAnnulus_minNhits = cms.uint32( 0 ),
    ElecPreIDLeadTkMatch_maxDR = cms.double( 0.015 ),
    EcalStripSumE_minClusEnergy = cms.double( 0.0 ),
    EcalStripSumE_deltaEta = cms.double( 0.0 ),
    EcalStripSumE_deltaPhiOverQ_minValue = cms.double( 0.0 ),
    EcalStripSumE_deltaPhiOverQ_maxValue = cms.double( 0.0 ),
    maximumForElectrionPreIDOutput = cms.double( 0.0 ),
    DataType = cms.string( "AOD" ),
    emMergingAlgorithm = cms.string( "None" ),
    candOverlapCriterion = cms.string( "None" ),
    doOneProng = cms.bool( True ),
    doOneProngStrip = cms.bool( True ),
    doOneProngTwoStrips = cms.bool( True ),
    doThreeProng = cms.bool( True ),
    tauPtThreshold = cms.double( 0.0 ),
    leadPionThreshold = cms.double( 1.0 ),
    stripPtThreshold = cms.double( 0.5 ),
    chargeHadrIsolationConeSize = cms.double( 0.5 ),
    gammaIsolationConeSize = cms.double( 0.5 ),
    neutrHadrIsolationConeSize = cms.double( 0.5 ),
    useIsolationAnnulus = cms.bool( False ),
    matchingCone = cms.double( 0.2 ),
    coneMetric = cms.string( "DR" ),
    coneSizeFormula = cms.string( "2.8/ET" ),
    minimumSignalCone = cms.double( 0.0 ),
    maximumSignalCone = cms.double( 1.8 ),
    oneProngStripMassWindow = cms.vdouble( 0.0, 0.0 ),
    oneProngTwoStripsMassWindow = cms.vdouble( 0.0, 0.0 ),
    oneProngTwoStripsPi0MassWindow = cms.vdouble( 0.0, 0.0 ),
    threeProngMassWindow = cms.vdouble( 0.0, 0.0 )
)
hltPFTauTrackFindingDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByLeadingObjectPtCut",
    PFTauProducer = cms.InputTag( "hltPFTaus" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) ),
    UseOnlyChargedHadrons = cms.bool( True ),
    MinPtLeadingObject = cms.double( 0.0 )
)
hltPFTauLooseIsolationDiscriminator = cms.EDProducer( "PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByTrackerIsolation = cms.bool( True ),
    ApplyDiscriminationByECALIsolation = cms.bool( False ),
    applyOccupancyCut = cms.bool( True ),
    maximumOccupancy = cms.uint32( 0 ),
    applySumPtCut = cms.bool( False ),
    maximumSumPtCut = cms.double( 6.0 ),
    applyRelativeSumPtCut = cms.bool( False ),
    relativeSumPtCut = cms.double( 0.0 ),
    PVProducer = cms.InputTag( "hltPixelVertices" ),
    customOuterCone = cms.double( -1.0 ),
    qualityCuts = cms.PSet( 
      isolationQualityCuts = cms.PSet( 
        minTrackHits = cms.uint32( 3 ),
        minTrackPt = cms.double( 1.5 ),
        maxTrackChi2 = cms.double( 100.0 ),
        minTrackPixelHits = cms.uint32( 0 ),
        minGammaEt = cms.double( 1.5 ),
        useTracksInsteadOfPFHadrons = cms.bool( False ),
        maxDeltaZ = cms.double( 0.2 ),
        maxTransverseImpactParameter = cms.double( 0.05 )
      ),
      signalQualityCuts = cms.PSet( 
        maxDeltaZ = cms.double( 0.5 ),
        minTrackPt = cms.double( 0.0 ),
        maxTrackChi2 = cms.double( 1000.0 ),
        useTracksInsteadOfPFHadrons = cms.bool( False ),
        minGammaEt = cms.double( 0.5 ),
        minTrackPixelHits = cms.uint32( 0 ),
        minTrackHits = cms.uint32( 3 ),
        maxTransverseImpactParameter = cms.double( 0.2 )
      )
    ),
    PFTauProducer = cms.InputTag( "hltPFTaus" ),
    Prediscriminants = cms.PSet(  BooleanOperator = cms.string( "and" ) )
)
hltSelectedPFTausTrackFinding = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTaus" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTrackFindingDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltSelectedPFTausTrackFindingLooseIsolation = cms.EDFilter( "PFTauSelector",
    src = cms.InputTag( "hltPFTaus" ),
    discriminators = cms.VPSet( 
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauTrackFindingDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      ),
      cms.PSet(  discriminator = cms.InputTag( "hltPFTauLooseIsolationDiscriminator" ),
        selectionCut = cms.double( 0.5 )
      )
    )
)
hltConvPFTausTrackFinding = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTrackFinding" )
)
hltConvPFTausTrackFindingLooseIsolation = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltSelectedPFTausTrackFindingLooseIsolation" )
)
hltConvPFTaus = cms.EDProducer( "PFTauToJetProducer",
    Source = cms.InputTag( "hltPFTaus" )
)
hltPFTau15 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTaus" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau15Track = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTrackFinding" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau15TrackLooseIso = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTrackFindingLooseIsolation" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltOverlapFilterMu15IsoPFTau15 = cms.EDFilter( "HLT2MuonTau",
    inputTag1 = cms.InputTag( "hltSingleMu15L3Filtered15" ),
    inputTag2 = cms.InputTag( "hltPFTau15TrackLooseIso" ),
    saveTags = cms.bool( True ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 1000.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 1000.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 1000.0 ),
    MinN = cms.int32( 1 )
)
hltPreMu17CentralJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1Mu7CenJetL1MuFiltered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu7" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1Mu7CenJetL2MuFiltered7 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu7CenJetL1MuFiltered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 7.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1Mu7CenJetL3MuFiltered17 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu7CenJetL2MuFiltered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 17.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltPreMu17DiCentralJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltDiJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltPreMu17TriCentralJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTriJet30Central = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 3 )
)
hltPreMu17Ele8CaloIdL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenEG5L1Filtered12 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 12.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1MuOpenEG5L2Filtered12 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L1Filtered12" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1MuOpenEG5L3Filtered17 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L2Filtered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 17.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1NonIsoHLTCaloIdLMu17Ele8ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoMu17Ele8HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLMu17Ele8ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoMu17Ele8HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPreMu12BTagIPCenJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltBJet30Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
hltGetJetsfromBJet30Central = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltBJet30Central" )
)
hltSelectorJetsSingleTop = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfromBJet30Central" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBLifetimeL25JetsSingleTop = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelectorJetsSingleTop" ),
    filter = cms.bool( False ),
    etMin = cms.double( 20.0 )
)
hltBLifetimeL25AssociatorSingleTop = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsSingleTop" ),
    tracks = cms.InputTag( "hltPixelTracks" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL25TagInfosSingleTop = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL25AssociatorSingleTop" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL25BJetTagsSingleTop = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL25TagInfosSingleTop' )
)
hltBLifetimeL25FilterSingleTop = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL25BJetTagsSingleTop" ),
    MinTag = cms.double( 0.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltL1Mu7CenJetL3MuFiltered12 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu7CenJetL2MuFiltered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltBLifetimeL3AssociatorSingleTop = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsSingleTop" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracksSingleTop" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfosSingleTop = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3AssociatorSingleTop" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTagsSingleTop = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfosSingleTop' )
)
hltBLifetimeL3FilterSingleTop = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTagsSingleTop" ),
    MinTag = cms.double( 3.3 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltPreMu17BTagIPCenJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sL1Mu0HTT50 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_Mu0_HTT50" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreMu15HT200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1Mu0HTT50L1MuFiltered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1Mu0HTT50" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1Mu0HTT50L2MuFiltered10 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L1MuFiltered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 10.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1Mu0HTT50L3MuFiltered15 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L2MuFiltered10" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 15.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreMu20HT200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1Mu0HTT50L2MuFiltered12 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L1MuFiltered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 12.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1Mu0HTT50L3MuFiltered20 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L2MuFiltered12" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 20.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreIsoMu15IsoPFTau15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltOverlapFilterIsoMu15IsoPFTau15 = cms.EDFilter( "HLT2MuonTau",
    inputTag1 = cms.InputTag( "hltSingleMuIsoL3IsoFiltered15" ),
    inputTag2 = cms.InputTag( "hltPFTau15TrackLooseIso" ),
    saveTags = cms.bool( True ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 1000.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 1000.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 1000.0 ),
    MinN = cms.int32( 1 )
)
hltPreIsoMu17BTagIPCentJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltIsoMu7CenJet30L2IsoFiltered7 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu7CenJetL2MuFiltered7" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    DepTag = cms.VInputTag( 'hltL2MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltIsoMu17CenJet30L3Filtered17 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltIsoMu7CenJet30L2IsoFiltered7" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 17.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( False )
)
hltIsoMu17CenJet30L3IsoFiltered17 = cms.EDFilter( "HLTMuonIsoFilter",
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltIsoMu17CenJet30L3Filtered17" ),
    MinN = cms.int32( 1 ),
    saveTags = cms.bool( True ),
    DepTag = cms.VInputTag( 'hltL3MuonIsolations' ),
    IsolatorPSet = cms.PSet(  )
)
hltPreDoubleMu3HT150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1Mu0HTT50L1DiMuFiltered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1Mu0HTT50" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1Mu0HTT50L2DiMuFiltered0 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L1DiMuFiltered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1Mu0HTT50L3DiMuFiltered3 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1Mu0HTT50L2DiMuFiltered0" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreDoubleMu3HT200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreDoubleMu5Ele8 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuOpenEG5L1DiMuFiltered3 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1MuOpenEG5" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 3.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
hltL1MuOpenEG5L2DiMuFiltered3 = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L1DiMuFiltered3" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltL1MuOpenEG5L3DiMuFiltered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltL1MuOpenEG5L2DiMuFiltered3" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 2.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltEG8HEFilterMuOpenEG5 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilterMuOpenEG5" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG8PixelMatchFilterMuOpenEG5 = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG8HEFilterMuOpenEG5" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltPrePhoton40R014MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR005MR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.05 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPrePhoton40R014MR450 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhoton40R020MR300 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR020MR300 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.2 ),
    minMR = cms.double( 300.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPrePhoton40R025MR200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhoton40R038MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltR038MR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.38 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreDoublePhoton40MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG40HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG40EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltDoubleIsoEG40EtFilterUnseededTight = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEcalActivitySuperClusterWrapper" ),
    etcutEB = cms.double( 40.0 ),
    etcutEE = cms.double( 40.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltDoublePhoton40EgammaLHEDoubleFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleIsoEG40EtFilterUnseededTight" ),
    isoTag = cms.InputTag( "hltActivityPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltRecoEcalSuperClusterActivityCandidate" ),
    L1NonIsoCand = cms.InputTag( "" )
)
hltMR150 = cms.EDFilter( "HLTRFilter",
    inputTag = cms.InputTag( "hltRHemisphere" ),
    inputMetTag = cms.InputTag( "hltMet" ),
    minR = cms.double( 0.0 ),
    minMR = cms.double( 150.0 ),
    doRPrime = cms.bool( False ),
    acceptNJ = cms.bool( True )
)
hltPreDoublePhoton40R014MR150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreEle8CaloIdTTrkIdTDiJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle8CaloIdTTrkIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG8EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdTTrkIdTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdTTrkIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.1 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdTTrkIdTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle8CaloIdTTrkIdTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle8CaloIdTTrkIdTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle8CaloIdTTrkIdTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle8CaloIdTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdTTrkIdTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle8CaloIdTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle8CaloIdTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle8CaloIdLCaloIsoVLJet40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltAntiKT5L2L3CaloJetsEle8CaloIdLCaloIsoVLRemoved = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle8CaloIdLCaloIsoVLPixelMatchFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.5 )
)
hltJet40Ele8CaloIdLCaloIsoVLRemoved = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltAntiKT5L2L3CaloJetsEle8CaloIdLCaloIsoVLRemoved" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 5.0 ),
    MinN = cms.int32( 1 )
)
hltL1sL1EG5HTT75 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_EG5_HTT75" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLHT200 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1EG5HTT75 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1EG5HTT75" ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltEG15EtFilterL1EG5HTT75 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1EG5HTT75" ),
    etcutEB = cms.double( 15.0 ),
    etcutEE = cms.double( 15.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTClusterShapeFilterEG5HTT75 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15EtFilterL1EG5HTT75" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTCaloIsoVLEcalIsoFilterEG5HTT75 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTClusterShapeFilterEG5HTT75" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTCaloIsoVLHEFilterEG5HTT75 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTCaloIsoVLEcalIsoFilterEG5HTT75" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.1 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTCaloIsoVLHcalIsoFilterEG5HTT75 = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTCaloIsoVLHEFilterEG5HTT75" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.2 ),
    thrOverEEE = cms.double( 0.2 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEG15CaloIdTCaloIsoVLPixelMatchFilterEG5HTT75 = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTCaloIsoVLHcalIsoFilterEG5HTT75" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdTCaloIsoVLOneOEMinusOneOPFilterEG5HTT75 = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEG15CaloIdTCaloIsoVLPixelMatchFilterEG5HTT75" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle15CaloIdTCaloIsoVLTrkIdTDetaFilterEG5HTT75 = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdTCaloIsoVLOneOEMinusOneOPFilterEG5HTT75" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle15CaloIdTCaloIsoVLTrkIdTDphiFilterEG5HTT75 = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdTCaloIsoVLTrkIdTDetaFilterEG5HTT75" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLTrackIsoFilterEG5HTT75 = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdTCaloIsoVLTrkIdTDphiFilterEG5HTT75" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.2 ),
    thrOverPtEE = cms.double( 0.2 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLHT250 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreEle15CaloIdVTTrkIdTLooseIsoPFTau15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle15CaloIdVTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG15CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle15CaloIdVTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle15CaloIdVTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle15CaloIdVTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle15CaloIdVTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle15CaloIdVTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltOverlapFilterEle15CaloJet5 = cms.EDFilter( "HLT2ElectronTau",
    inputTag1 = cms.InputTag( "hltEle15CaloIdVTTrkIdTDphiFilter" ),
    inputTag2 = cms.InputTag( "hltTauJet5" ),
    saveTags = cms.bool( False ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 9999.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 9999.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 99999.0 ),
    MinN = cms.int32( 1 )
)
hltOverlapFilterEle15IsoPFTau15 = cms.EDFilter( "HLT2ElectronTau",
    inputTag1 = cms.InputTag( "hltEle15CaloIdVTTrkIdTDphiFilter" ),
    inputTag2 = cms.InputTag( "hltPFTau15TrackLooseIso" ),
    saveTags = cms.bool( True ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 1000.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 1000.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 1000.0 ),
    MinN = cms.int32( 1 )
)
hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTLooseIsoPFTau15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltOverlapFilterIsoEle15CaloJet5 = cms.EDFilter( "HLT2ElectronTau",
    inputTag1 = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter" ),
    inputTag2 = cms.InputTag( "hltTauJet5" ),
    saveTags = cms.bool( False ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 9999.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 9999.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 99999.0 ),
    MinN = cms.int32( 1 )
)
hltOverlapFilterIsoEle15IsoPFTau15 = cms.EDFilter( "HLT2ElectronTau",
    inputTag1 = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter" ),
    inputTag2 = cms.InputTag( "hltPFTau15TrackLooseIso" ),
    saveTags = cms.bool( True ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 1000.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 1000.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 1000.0 ),
    MinN = cms.int32( 1 )
)
hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTLooseIsoPFTau20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPFJet20 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltAntiKT5ConvPFJetsForTaus" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau20 = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTaus" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau20Track = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTrackFinding" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltPFTau20TrackLooseIso = cms.EDFilter( "HLT1Tau",
    inputTag = cms.InputTag( "hltConvPFTausTrackFindingLooseIsolation" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.5 ),
    MinN = cms.int32( 1 )
)
hltOverlapFilterIsoEle15IsoPFTau20 = cms.EDFilter( "HLT2ElectronTau",
    inputTag1 = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter" ),
    inputTag2 = cms.InputTag( "hltPFTau20TrackLooseIso" ),
    saveTags = cms.bool( True ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( 1000.0 ),
    MinDeta = cms.double( 0.0 ),
    MaxDeta = cms.double( 1000.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( 14000.0 ),
    MinDelR = cms.double( 0.3 ),
    MaxDelR = cms.double( 1000.0 ),
    MinN = cms.int32( 1 )
)
hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltCleanEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta3 = cms.EDFilter( "HLTJetVBFFilter",
    inputTag = cms.InputTag( "hltCleanEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    minEtLow = cms.double( 25.0 ),
    minEtHigh = cms.double( 35.0 ),
    etaOpposite = cms.bool( False ),
    minDeltaEta = cms.double( 3.0 ),
    minInvMass = cms.double( 0.0 ),
    maxEta = cms.double( -1.0 )
)
hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta2 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta2 = cms.EDFilter( "HLTJetVBFFilter",
    inputTag = cms.InputTag( "hltCleanEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    minEtLow = cms.double( 25.0 ),
    minEtHigh = cms.double( 35.0 ),
    etaOpposite = cms.bool( False ),
    minDeltaEta = cms.double( 2.0 ),
    minInvMass = cms.double( 0.0 ),
    maxEta = cms.double( -1.0 )
)
hltPreEle15CaloIdVTTrkIdTJet35Jet25Deta2 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltCleanEle15CaloIdVTTrkIdTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle15CaloIdVTTrkIdTDphiFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle15CaloIdVTTrkIdTJet35Jet25Deta2 = cms.EDFilter( "HLTJetVBFFilter",
    inputTag = cms.InputTag( "hltCleanEle15CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    minEtLow = cms.double( 25.0 ),
    minEtHigh = cms.double( 35.0 ),
    etaOpposite = cms.bool( False ),
    minDeltaEta = cms.double( 2.0 ),
    minInvMass = cms.double( 0.0 ),
    maxEta = cms.double( -1.0 )
)
hltPreEle25CaloIdVTTrkIdTCentralJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG25EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1SingleEG12" ),
    etcutEB = cms.double( 25.0 ),
    etcutEE = cms.double( 25.0 ),
    ncandcut = cms.int32( 1 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG25EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle25CaloIdVTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle25CaloIdVTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTDphiFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle25CaloIdVTTrkIdTCentralJet30Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltPreEle25CaloIdVTTrkIdTCentralDiJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle25CaloIdVTTrkIdTCentralDiJet30Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltPreEle25CaloIdVTTrkIdTCentralTriJet30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle25CaloIdVTTrkIdTCentralTriJet30Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 3 )
)
hltPreEle25CaloIdVTTrkIdTCentralJet30BTagIP = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltSingleEleCleanBJet30Central = cms.EDFilter( "HLT1CaloBJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
hltGetJetsfrom1EleCleanBJet30Central = cms.EDProducer( "GetJetsFromHLTobject",
    jets = cms.InputTag( "hltSingleEleCleanBJet30Central" )
)
hltSelectorEleJetsSingleTop = cms.EDFilter( "LargestEtCaloJetSelector",
    src = cms.InputTag( "hltGetJetsfrom1EleCleanBJet30Central" ),
    filter = cms.bool( False ),
    maxNumber = cms.uint32( 4 )
)
hltBLifetimeL25JetsEleJetSingleTop = cms.EDFilter( "EtMinCaloJetSelector",
    src = cms.InputTag( "hltSelectorEleJetsSingleTop" ),
    filter = cms.bool( False ),
    etMin = cms.double( 20.0 )
)
hltBLifetimeL25AssociatorEleJetSingleTop = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsEleJetSingleTop" ),
    tracks = cms.InputTag( "hltPixelTracks" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL25TagInfosEleJetSingleTop = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL25AssociatorEleJetSingleTop" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 3 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 5.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL25BJetTagsEleJetSingleTop = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL25TagInfosEleJetSingleTop' )
)
hltBLifetimeL25FilterEleJetSingleTop = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL25BJetTagsEleJetSingleTop" ),
    MinTag = cms.double( 0.0 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( False )
)
hltBLifetimeL3AssociatorEleJetSingleTop = cms.EDProducer( "JetTracksAssociatorAtVertex",
    jets = cms.InputTag( "hltBLifetimeL25JetsEleJetSingleTop" ),
    tracks = cms.InputTag( "hltBLifetimeRegionalCtfWithMaterialTracksEleJetSingleTop" ),
    coneSize = cms.double( 0.5 )
)
hltBLifetimeL3TagInfosEleJetSingleTop = cms.EDProducer( "TrackIPProducer",
    jetTracks = cms.InputTag( "hltBLifetimeL3AssociatorEleJetSingleTop" ),
    primaryVertex = cms.InputTag( "hltPixelVertices" ),
    computeProbabilities = cms.bool( False ),
    computeGhostTrack = cms.bool( False ),
    ghostTrackPriorDeltaR = cms.double( 0.03 ),
    minimumNumberOfPixelHits = cms.int32( 2 ),
    minimumNumberOfHits = cms.int32( 8 ),
    maximumTransverseImpactParameter = cms.double( 0.2 ),
    minimumTransverseMomentum = cms.double( 1.0 ),
    maximumChiSquared = cms.double( 20.0 ),
    maximumLongitudinalImpactParameter = cms.double( 17.0 ),
    jetDirectionUsingTracks = cms.bool( False ),
    jetDirectionUsingGhostTrack = cms.bool( False ),
    useTrackQuality = cms.bool( False )
)
hltBLifetimeL3BJetTagsEleJetSingleTop = cms.EDProducer( "JetTagProducer",
    jetTagComputer = cms.string( "hltESPTrackCounting3D2nd" ),
    tagInfos = cms.VInputTag( 'hltBLifetimeL3TagInfosEleJetSingleTop' )
)
hltBLifetimeL3FilterEleJetSingleTop = cms.EDFilter( "HLTJetTag",
    JetTag = cms.InputTag( "hltBLifetimeL3BJetTagsEleJetSingleTop" ),
    MinTag = cms.double( 3.3 ),
    MaxTag = cms.double( 99999.0 ),
    MinJets = cms.int32( 1 ),
    saveTags = cms.bool( True )
)
hltPreEle17CaloIdVTTrkIdTCenJet30CenJet25 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG17CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTHEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle17CaloIdVTTrkIdTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTTrkIdTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltCleanEle17CaloIdVTTrkIdTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTDphiFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle17CaloIdVTTrkIdTCentralJet30Cleaned  = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle17CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltEle17CaloIdVTTrkIdTCentralDiJet25Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle17CaloIdVTTrkIdTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltPreEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCenJet30CenJet25PFMHT15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG17CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltCleanEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralJet30Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralDiJet25Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltPFMHT15Filter = cms.EDFilter( "HLTMhtHtFilter",
    inputJetTag = cms.InputTag( "hltAntiKT5ConvPFJets" ),
    saveTags = cms.bool( True ),
    minMht = cms.double( 15.0 ),
    minNJet = cms.int32( 0 ),
    mode = cms.int32( 1 ),
    usePt = cms.bool( True ),
    minPT12 = cms.double( 0.0 ),
    minMeff = cms.double( 0.0 ),
    minHt = cms.double( 0.0 ),
    minAlphaT = cms.double( 0.0 ),
    minPtJet = cms.vdouble( 5.0, 5.0 ),
    etaJet = cms.vdouble( 9999.0, 9999.0 )
)
hltPreEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCenJet30CenJet25PFMHT20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEG25CaloIdTClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG25EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEG25CaloIdTClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonEcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonEcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.05 ),
    thrOverEEE = cms.double( 0.05 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalIsol" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.125 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.0080 ),
    thrRegularEE = cms.double( 0.0080 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.07 ),
    thrRegularEE = cms.double( 0.05 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter" ),
    isoTag = cms.InputTag( "hltL1IsoElectronTrackIsol" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoElectronTrackIsol" ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverPtEB = cms.double( 0.125 ),
    thrOverPtEE = cms.double( 0.075 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 1 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltCleanEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets = cms.EDProducer( "JetCollectionForEleHT",
    HltElectronTag = cms.InputTag( "hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter" ),
    SourceJetTag = cms.InputTag( "hltAntiKT5L2L3CorrCaloJets" ),
    minDeltaR = cms.double( 0.3 )
)
hltEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralJet30Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 1 )
)
hltEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralDiJet25Cleaned = cms.EDFilter( "HLT1CaloJet",
    inputTag = cms.InputTag( "hltCleanEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets" ),
    saveTags = cms.bool( True ),
    MinPt = cms.double( 25.0 ),
    MaxEta = cms.double( 2.6 ),
    MinN = cms.int32( 2 )
)
hltL1sL1DoubleEG5HTT50 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleEG5_HTT50" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleEle8CaloIdLTrkIdVLHT150 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1DoubleEG5HTT50 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1DoubleEG5HTT50" ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltDoubleEG8EtFilterL1DoubleEG5HTT50 = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1DoubleEG5HTT50" ),
    etcutEB = cms.double( 8.0 ),
    etcutEE = cms.double( 8.0 ),
    ncandcut = cms.int32( 2 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8HTT50ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterL1DoubleEG5HTT50" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8HTT50HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8HTT50ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLDoubleEle8HTT50PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8HTT50HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLDoubleEle8HTT50PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50DetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreDoubleEle8CaloIdTTrkIdVLHT160 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1NonIsoHLTCaloIdTDoubleEle8HTT50ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltDoubleEG8EtFilterL1DoubleEG5HTT50" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.011 ),
    thrRegularEE = cms.double( 0.031 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdTDoubleEle8HTT50HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdTDoubleEle8HTT50ClusterShapeFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.1 ),
    thrOverEEE = cms.double( 0.075 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdTDoubleEle8HTT50PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdTDoubleEle8HTT50HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdTDoubleEle8HTT50PixelMatchFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50DetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50DphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50DetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1sL1TripleEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_TripleEG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreDoubleEle10CaloIdLTrkIdVLEle10 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltEGRegionalL1TripleEG5 = cms.EDFilter( "HLTEgammaL1MatchFilterRegional",
    candIsolatedTag = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    l1IsolatedTag = cms.InputTag( 'l1extraParticles','Isolated' ),
    candNonIsolatedTag = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" ),
    l1NonIsolatedTag = cms.InputTag( 'l1extraParticles','NonIsolated' ),
    L1SeedFilterTag = cms.InputTag( "hltL1sL1TripleEG5" ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    region_eta_size = cms.double( 0.522 ),
    region_eta_size_ecap = cms.double( 1.0 ),
    region_phi_size = cms.double( 1.044 ),
    barrel_end = cms.double( 1.4791 ),
    endcap_end = cms.double( 2.65 )
)
hltTripleEG10EtFilter = cms.EDFilter( "HLTEgammaEtFilter",
    inputTag = cms.InputTag( "hltEGRegionalL1TripleEG5" ),
    etcutEB = cms.double( 10.0 ),
    etcutEE = cms.double( 10.0 ),
    ncandcut = cms.int32( 3 ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoTripleElectronEt10HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltTripleEG10EtFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLTNonIsoTripleElectronEt10PixelMatchFilter = cms.EDFilter( "HLTElectronPixelMatchFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoTripleElectronEt10HEFilter" ),
    L1IsoPixelSeedsTag = cms.InputTag( "hltL1IsoStartUpElectronPixelSeeds" ),
    L1NonIsoPixelSeedsTag = cms.InputTag( "hltL1NonIsoStartUpElectronPixelSeeds" ),
    npixelmatchcut = cms.double( 1.0 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLT2CaloIdLTripleElectronEt10HEFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoTripleElectronEt10PixelMatchFilter" ),
    isoTag = cms.InputTag( "hltL1IsolatedPhotonHcalForHE" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsolatedPhotonHcalForHE" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( -1.0 ),
    thrRegularEE = cms.double( -1.0 ),
    thrOverEEB = cms.double( 0.15 ),
    thrOverEEE = cms.double( 0.1 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLT2LegEleIdTripleElectronEt10ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLT2CaloIdLTripleElectronEt10HEFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLT2LegEleIdTripleElectronEt10OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltL1NonIsoHLT2LegEleIdTripleElectronEt10ClusterShapeFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltL1NonIsoHLT2LegEleIdTripleElectronEt10EleIdDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLT2LegEleIdTripleElectronEt10OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLT2LegEleIdTripleElectronEt10EleIdDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLT2LegEleIdTripleElectronEt10EleIdDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 2 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltPreTripleEle10CaloIdLTrkIdVL = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1NonIsoHLT3LegEleIdTripleElectronEt10ClusterShapeFilter = cms.EDFilter( "HLTEgammaGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLTNonIsoTripleElectronEt10PixelMatchFilter" ),
    isoTag = cms.InputTag( "hltL1IsoHLTClusterShape" ),
    nonIsoTag = cms.InputTag( "hltL1NonIsoHLTClusterShape" ),
    lessThan = cms.bool( True ),
    useEt = cms.bool( False ),
    thrRegularEB = cms.double( 0.014 ),
    thrRegularEE = cms.double( 0.035 ),
    thrOverEEB = cms.double( -1.0 ),
    thrOverEEE = cms.double( -1.0 ),
    thrOverE2EB = cms.double( -1.0 ),
    thrOverE2EE = cms.double( -1.0 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltL1IsoRecoEcalCandidate" ),
    L1NonIsoCand = cms.InputTag( "hltL1NonIsoRecoEcalCandidate" )
)
hltL1NonIsoHLT3LegEleIdTripleElectronEt10OneOEMinusOneOPFilter = cms.EDFilter( "HLTElectronOneOEMinusOneOPFilterRegional",
    candTag = cms.InputTag( "hltL1NonIsoHLT3LegEleIdTripleElectronEt10ClusterShapeFilter" ),
    electronIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    electronNonIsolatedProducer = cms.InputTag( "hltPixelMatchElectronsL1NonIso" ),
    barrelcut = cms.double( 999.9 ),
    endcapcut = cms.double( 999.9 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False )
)
hltL1NonIsoHLT3LegEleIdTripleElectronEt10EleIdDetaFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLT3LegEleIdTripleElectronEt10OneOEMinusOneOPFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Deta' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Deta' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.01 ),
    thrRegularEE = cms.double( 0.01 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( False ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1NonIsoHLT3LegEleIdTripleElectronEt10EleIdDphiFilter = cms.EDFilter( "HLTElectronGenericFilter",
    candTag = cms.InputTag( "hltL1NonIsoHLT3LegEleIdTripleElectronEt10EleIdDetaFilter" ),
    isoTag = cms.InputTag( 'hltElectronL1IsoDetaDphi','Dphi' ),
    nonIsoTag = cms.InputTag( 'hltElectronL1NonIsoDetaDphi','Dphi' ),
    lessThan = cms.bool( True ),
    thrRegularEB = cms.double( 0.15 ),
    thrRegularEE = cms.double( 0.1 ),
    thrOverPtEB = cms.double( -1.0 ),
    thrOverPtEE = cms.double( -1.0 ),
    thrTimesPtEB = cms.double( -1.0 ),
    thrTimesPtEE = cms.double( -1.0 ),
    ncandcut = cms.int32( 3 ),
    doIsolated = cms.bool( False ),
    saveTags = cms.bool( True ),
    L1IsoCand = cms.InputTag( "hltPixelMatchElectronsL1Iso" ),
    L1NonIsoCand = cms.InputTag( "hltPixelMatchElectronsL1NonIso" )
)
hltL1sETT220 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ETT220" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPrePixelTracksMultiplicity80 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPixelClusterShapeFilter = cms.EDFilter( "HLTPixelClusterShapeFilter",
    inputTag = cms.InputTag( "hltSiPixelRecHits" ),
    saveTags = cms.bool( False ),
    minZ = cms.double( -10.1 ),
    maxZ = cms.double( 10.1 ),
    zStep = cms.double( 0.2 ),
    nhitsTrunc = cms.int32( 150 ),
    clusterTrunc = cms.double( 2.0 ),
    clusterPars = cms.vdouble( 0.0, 0.0045 )
)
hltPixelVerticesForHighMult = cms.EDProducer( "PixelVertexProducer",
    Verbosity = cms.int32( 0 ),
    Finder = cms.string( "DivisiveVertexFinder" ),
    UseError = cms.bool( True ),
    WtAverage = cms.bool( True ),
    ZOffset = cms.double( 5.0 ),
    ZSeparation = cms.double( 0.05 ),
    NTrkMin = cms.int32( 50 ),
    PtMin = cms.double( 0.4 ),
    TrackCollection = cms.InputTag( "hltPixelTracksForHighMult" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    Method2 = cms.bool( True )
)
hltPixelCandsForHighMult = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltPixelTracksForHighMult" ),
    particleType = cms.string( "pi+" )
)
hlt1HighMult80 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltPixelVerticesForHighMult" ),
    trackCollection = cms.InputTag( "hltPixelCandsForHighMult" ),
    MinPt = cms.double( 0.4 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 10.0 ),
    MinTrks = cms.int32( 80 ),
    MinSep = cms.double( 0.05 )
)
hltPrePixelTracksMultiplicity100 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hlt1HighMult100 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltPixelVerticesForHighMult" ),
    trackCollection = cms.InputTag( "hltPixelCandsForHighMult" ),
    MinPt = cms.double( 0.4 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 10.0 ),
    MinTrks = cms.int32( 100 ),
    MinSep = cms.double( 0.05 )
)
hltL1sL1BeamGasHf = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BeamGas_Hf" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1BeamGasHf = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltHFAsymmetryFilter = cms.EDFilter( "HLTHFAsymmetryFilter",
    HFHitCollection = cms.InputTag( "hltHfreco" ),
    ECut_HF = cms.double( 3.0 ),
    OS_Asym_max = cms.double( 0.2 ),
    SS_Asym_min = cms.double( 0.8 )
)
hltL1sL1BeamGasBsc = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BeamGas_Bsc" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1BeamGasBsc = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPixelActivityFilter = cms.EDFilter( "HLTPixelActivityFilter",
    inputTag = cms.InputTag( "hltSiPixelClusters" ),
    saveTags = cms.bool( False ),
    minClusters = cms.uint32( 3 ),
    maxClusters = cms.uint32( 0 )
)
hltPixelAsymmetryFilter = cms.EDFilter( "HLTPixelAsymmetryFilter",
    inputTag = cms.InputTag( "hltSiPixelClusters" ),
    saveTags = cms.bool( False ),
    MinAsym = cms.double( 0.0 ),
    MaxAsym = cms.double( 1.0 ),
    MinCharge = cms.double( 4000.0 ),
    MinBarrel = cms.double( 10000.0 )
)
hltL1sL1BeamHalo = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BeamHalo" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1BeamHalo = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPixelActivityFilterForHalo = cms.EDFilter( "HLTPixelActivityFilter",
    inputTag = cms.InputTag( "hltSiPixelClusters" ),
    saveTags = cms.bool( False ),
    minClusters = cms.uint32( 0 ),
    maxClusters = cms.uint32( 10 )
)
hltTrackerHaloFilter = cms.EDFilter( "HLTTrackerHaloFilter",
    inputTag = cms.InputTag( "hltSiStripClusters" ),
    saveTags = cms.bool( False ),
    MaxClustersTECp = cms.int32( 50 ),
    MaxClustersTECm = cms.int32( 50 ),
    SignalAccumulation = cms.int32( 5 ),
    MaxClustersTEC = cms.int32( 60 ),
    MaxAccus = cms.int32( 4 ),
    FastProcessing = cms.int32( 1 )
)
hltL1sZeroBias = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ZeroBias" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1TechBSCminBiasthreshold1 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1TechBSCminBiasthreshold1 = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'L1Tech_BSC_minBias_threshold1' ),
    hltResults = cms.InputTag( "" ),
    l1tResults = cms.InputTag( "gtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( True )
)
hltPreL1TechBSChalo = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1TechBSChalo = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'L1Tech_BSC_halo_beam2_inner',
      'L1Tech_BSC_halo_beam2_outer',
      'L1Tech_BSC_halo_beam1_inner',
      'L1Tech_BSC_halo_beam1_outer' ),
    hltResults = cms.InputTag( "" ),
    l1tResults = cms.InputTag( "gtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( True )
)
hltL1sL1PreCollisions = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_PreCollisions" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1PreCollisions = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sL1InterbunchBsc = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_InterBunch_Bsc" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1Interbunch1 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sGlobalRunHPDNoise = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet20_NotBptxOR_NotMuBeamHalo" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreGlobalRunHPDNoise = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sTechTrigHCALNoise = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( True ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "11 OR 12" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreTechTrigHCALNoise = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreZeroBias = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhysics = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPrePhysicsNanoDST = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sTrackerCosmics = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( True ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "25" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreTrackerCosmics = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTrackerCosmicsPattern = cms.EDFilter( "HLTLevel1Pattern",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    triggerBit = cms.string( "L1Tech_RPC_TTU_pointing_Cosmics.v0" ),
    daqPartitions = cms.uint32( 1 ),
    ignoreL1Mask = cms.bool( False ),
    invert = cms.bool( False ),
    throw = cms.bool( True ),
    bunchCrossings = cms.vint32( -2, -1, 0, 1, 2 ),
    triggerPattern = cms.vint32( 1, 1, 1, 0, 0 )
)
hltPreRegionalCosmicTracking = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1MuORL1Filtered0 = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "l1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMuOpenCandidate" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
hltSingleL2MuORL2PreFilteredNoVtx = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidatesNoVtx" ),
    PreviousCandTag = cms.InputTag( "hltL1MuORL1Filtered0" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltRegionalCosmicTrackerSeeds = cms.EDProducer( "SeedGeneratorFromRegionHitsEDProducer",
    ClusterCheckPSet = cms.PSet( 
      MaxNumberOfPixelClusters = cms.uint32( 300000 ),
      PixelClusterCollectionLabel = cms.InputTag( "hltSiPixelClusters" ),
      MaxNumberOfCosmicClusters = cms.uint32( 300000 ),
      ClusterCollectionLabel = cms.InputTag( "hltSiStripClusters" ),
      doClusterCheck = cms.bool( False )
    ),
    RegionFactoryPSet = cms.PSet( 
      CollectionsPSet = cms.PSet( 
        recoMuonsCollection = cms.InputTag( "muons" ),
        recoTrackMuonsCollection = cms.InputTag( "cosmicMuons" ),
        recoL2MuonsCollection = cms.InputTag( "hltL2MuonCandidatesNoVtx" )
      ),
      ComponentName = cms.string( "CosmicRegionalSeedGenerator" ),
      RegionInJetsCheckPSet = cms.PSet( 
        recoCaloJetsCollection = cms.InputTag( "hltIterativeCone5CaloJets" ),
        deltaRExclusionSize = cms.double( 0.3 ),
        jetsPtMin = cms.double( 5.0 ),
        doJetsExclusionCheck = cms.bool( False )
      ),
      ToolsPSet = cms.PSet( 
        regionBase = cms.string( "seedOnL2Muon" ),
        thePropagatorName = cms.string( "hltESPAnalyticalPropagator" )
      ),
      RegionPSet = cms.PSet( 
        precise = cms.bool( False ),
        deltaPhiRegion = cms.double( 0.3 ),
        measurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
        zVertex = cms.double( 5.0 ),
        deltaEtaRegion = cms.double( 0.3 ),
        ptMin = cms.double( 5.0 ),
        rVertex = cms.double( 5.0 )
      )
    ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "GenericPairGenerator" ),
      LayerPSet = cms.PSet( 
        TOB = cms.PSet(  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ) ),
        layerList = cms.vstring( 'TOB6+TOB5',
          'TOB6+TOB4',
          'TOB6+TOB3',
          'TOB5+TOB4',
          'TOB5+TOB3',
          'TOB4+TOB3' )
      )
    ),
    SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
    SeedCreatorPSet = cms.PSet( 
      ComponentName = cms.string( "CosmicSeedCreator" ),
      maxseeds = cms.int32( 100000 ),
      propagator = cms.string( "PropagatorWithMaterial" )
    ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" )
)
hltRegionalCosmicCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltRegionalCosmicTrackerSeeds" ),
    TrajectoryBuilder = cms.string( "hltESPCkfTrajectoryBuilder" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    NavigationSchool = cms.string( "CosmicNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    useHitsSplitting = cms.bool( False ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterial" ),
      numberMeasurementsForFit = cms.int32( 4 ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOpposite" )
    ),
    cleanTrajectoryAfterInOut = cms.bool( False ),
    maxNSeeds = cms.uint32( 200 )
)
hltRegionalCosmicTracks = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "hltRegionalCosmicTracks" ),
    Fitter = cms.string( "hltESPKFFittingSmoother" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    src = cms.InputTag( "hltRegionalCosmicCkfTrackCandidates" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltCosmicTrackSelector = cms.EDProducer( "CosmicTrackSelector",
    src = cms.InputTag( "hltRegionalCosmicTracks" ),
    beamspot = cms.InputTag( "offlineBeamSpot" ),
    keepAllTracks = cms.bool( False ),
    chi2n_par = cms.double( 10.0 ),
    max_d0 = cms.double( 999.0 ),
    max_z0 = cms.double( 999.0 ),
    min_pt = cms.double( 5.0 ),
    max_eta = cms.double( 2.0 ),
    min_nHit = cms.uint32( 6 ),
    min_nPixelHit = cms.uint32( 0 ),
    minNumberLayers = cms.uint32( 0 ),
    minNumber3DLayers = cms.uint32( 0 ),
    maxNumberLostLayers = cms.uint32( 999 ),
    qualityBit = cms.string( "" )
)
hltPreL3MuonsCosmicTracking = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL3TrajectorySeedNoVtx = cms.EDProducer( "TSGFromL2Muon",
    PtCut = cms.double( 1.0 ),
    PCut = cms.double( 2.5 ),
    MuonCollectionLabel = cms.InputTag( "hltL2Muons" ),
    ServiceParameters = cms.PSet( 
      Propagators = cms.untracked.vstring( 'hltESPSteppingHelixPropagatorOpposite',
        'hltESPSteppingHelixPropagatorAlong' ),
      RPCLayers = cms.bool( True ),
      UseMuonNavigation = cms.untracked.bool( True )
    ),
    MuonTrackingRegionBuilder = cms.PSet(  ),
    TkSeedGenerator = cms.PSet( 
      propagatorCompatibleName = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
      option = cms.uint32( 3 ),
      maxChi2 = cms.double( 40.0 ),
      errorMatrixPset = cms.PSet( 
        atIP = cms.bool( True ),
        action = cms.string( "use" ),
        errorMatrixValuesPSet = cms.PSet( 
          pf3_V12 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V13 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V11 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V14 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          yAxis = cms.vdouble( 0.0, 1.0, 1.4, 10.0 ),
          pf3_V34 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V15 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V33 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V45 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V44 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          xAxis = cms.vdouble( 0.0, 13.0, 30.0, 70.0, 1000.0 ),
          pf3_V23 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V22 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          pf3_V55 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 3.0, 3.0, 3.0, 5.0, 4.0, 5.0, 10.0, 7.0, 10.0, 10.0, 10.0, 10.0 )
          ),
          zAxis = cms.vdouble( -3.14159, 3.14159 ),
          pf3_V35 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V25 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          ),
          pf3_V24 = cms.PSet( 
            action = cms.string( "scale" ),
            values = cms.vdouble( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 )
          )
        )
      ),
      propagatorName = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
      manySeeds = cms.bool( False ),
      copyMuonRecHit = cms.bool( False ),
      ComponentName = cms.string( "TSGForRoadSearch" )
    ),
    TrackerSeedCleaner = cms.PSet(  ),
    TSGFromMixedPairs = cms.PSet(  ),
    TSGFromPixelTriplets = cms.PSet(  ),
    TSGFromPixelPairs = cms.PSet(  ),
    TSGForRoadSearchOI = cms.PSet(  ),
    TSGForRoadSearchIOpxl = cms.PSet(  ),
    TSGFromPropagation = cms.PSet(  ),
    TSGFromCombinedHits = cms.PSet(  )
)
hltL3TkTracksFromL2NoVtx = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( False ),
    useHitsSplitting = cms.bool( False ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "" ),
    Fitter = cms.string( "hltESPKFFittingSmoother" ),
    Propagator = cms.string( "PropagatorWithMaterial" ),
    src = cms.InputTag( "hltL3TrackCandidateFromL2NoVtx" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    AlgorithmName = cms.string( "undefAlgorithm" ),
    NavigationSchool = cms.string( "" )
)
hltL3MuonsNoVtx = cms.EDProducer( "L3TkMuonProducer",
    InputObjects = cms.InputTag( "hltL3TkTracksFromL2NoVtx" )
)
hltL3MuonCandidatesNoVtx = cms.EDProducer( "L3MuonCandidateProducer",
    InputObjects = cms.InputTag( "hltL3MuonsNoVtx" )
)
hltMu5NoVertexL3PreFiltered5 = cms.EDFilter( "HLTMuonL3PreFilter",
    BeamSpotTag = cms.InputTag( "offlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL3MuonCandidatesNoVtx" ),
    PreviousCandTag = cms.InputTag( "hltSingleL2MuORL2PreFilteredNoVtx" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 2.5 ),
    MinNhits = cms.int32( 6 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True )
)
hltPreLogMonitor = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltLogMonitorFilter = cms.EDFilter( "HLTLogMonitorFilter",
    default_threshold = cms.uint32( 10 ),
    categories = cms.VPSet( 
    )
)
hltPreL1ETM30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltPreL1DoubleJet36Central = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltL1sL1MultiJet = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HTT50 OR L1_TripleJet28_Central OR L1_QuadJet20_Central" ),
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "gtDigis" ),
    L1CollectionsTag = cms.InputTag( "l1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "l1extraParticles" ),
    saveTags = cms.bool( True )
)
hltPreL1MultiJet = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "gtDigis" ),
    offset = cms.uint32( 0 )
)
hltTriggerSummaryAOD = cms.EDProducer( "TriggerSummaryProducerAOD",
    processName = cms.string( "@" )
)
hltTriggerSummaryRAW = cms.EDProducer( "TriggerSummaryProducerRAW",
    processName = cms.string( "@" )
)
hltBoolTrue = cms.EDFilter( "HLTBool",
    result = cms.bool( True )
)
hltPixelVertices = cms.EDProducer( "PixelVertexProducer",
    Verbosity = cms.int32( 0 ),
    Finder = cms.string( "DivisiveVertexFinder" ),
    UseError = cms.bool( True ),
    WtAverage = cms.bool( True ),
    ZOffset = cms.double( 5.0 ),
    ZSeparation = cms.double( 0.05 ),
    NTrkMin = cms.int32( 2 ),
    PtMin = cms.double( 1.0 ),
    TrackCollection = cms.InputTag( "hltPixelTracks" ),
    beamSpot = cms.InputTag( "offlineBeamSpot" ),
    Method2 = cms.bool( True )
)

HLTEcalActivitySequence = cms.Sequence( hltEcalRawToRecHitFacility + hltESRawToRecHitFacility + hltEcalRegionalRestFEDs + hltEcalRegionalESRestFEDs + hltEcalRecHitAll + hltESRecHitAll + hltHybridSuperClustersActivity + hltCorrectedHybridSuperClustersActivity + hltMulti5x5BasicClustersActivity + hltMulti5x5SuperClustersActivity + hltMulti5x5SuperClustersWithPreshowerActivity + hltCorrectedMulti5x5SuperClustersWithPreshowerActivity + hltRecoEcalSuperClusterActivityCandidate + hltEcalActivitySuperClusterWrapper )
HLTDoLocalHcalSequence = cms.Sequence( hltHcalDigis + hltHbhereco + hltHfreco + hltHoreco )
HLTDoCaloSequence = cms.Sequence( hltEcalRawToRecHitFacility + hltEcalRegionalRestFEDs + hltEcalRecHitAll + HLTDoLocalHcalSequence + hltTowerMakerForAll )
HLTRecoJetSequenceAK5Uncorrected = cms.Sequence( HLTDoCaloSequence + hltAntiKT5CaloJets )
HLTRecoJetSequenceAK5Corrected = cms.Sequence( HLTRecoJetSequenceAK5Uncorrected + hltAntiKT5L2L3CorrCaloJets + hltJetIDPassedCorrJets )
HLTDoRegionalJetEcalSequence = cms.Sequence( hltEcalRawToRecHitFacility + hltEcalRegionalJetsFEDs + hltEcalRegionalJetsRecHit )
HLTRegionalTowerMakerForJetsSequence = cms.Sequence( HLTDoRegionalJetEcalSequence + HLTDoLocalHcalSequence + hltTowerMakerForJets )
HLTRegionalRecoJetSequenceAK5Corrected = cms.Sequence( HLTRegionalTowerMakerForJetsSequence + hltAntiKT5CaloJetsRegional + hltAntiKT5L2L3CorrCaloJetsRegional + hltL1MatchedJetsRegional + hltJetIDPassedJetsRegional )
HLTRecoMETSequence = cms.Sequence( HLTDoCaloSequence + hltMet )
HLTBtagIPSequenceL25Hbb = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltGetJetsfromBJetHbb + hltSelectorJetsHbb + hltBLifetimeL25JetsHbb + hltBLifetimeL25AssociatorHbb + hltBLifetimeL25TagInfosHbb + hltBLifetimeL25BJetTagsHbb )
HLTBtagIPSequenceL3Hbb = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGeneratorHbb + hltBLifetimeRegionalCkfTrackCandidatesHbb + hltBLifetimeRegionalCtfWithMaterialTracksHbb + hltBLifetimeL3AssociatorHbb + hltBLifetimeL3TagInfosHbb + hltBLifetimeL3BJetTagsHbb )
HLTRecoJetSequencePrePF = cms.Sequence( HLTRecoJetSequenceAK5Uncorrected + hltAntiKT5CaloJetsEt5 )
HLTTrackReconstructionForJets = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + HLTDoLocalStripSequence + hltPFJetPixelSeeds + hltPFJetCkfTrackCandidates + hltPFJetCtfWithMaterialTracks + hltPFlowTrackSelectionHighPurity )
HLTPreshowerSequence = cms.Sequence( hltESRawToRecHitFacility + hltEcalRegionalESRestFEDs + hltESRecHitAll )
HLTParticleFlowSequenceForTaus = cms.Sequence( HLTPreshowerSequence + hltParticleFlowRecHitECAL + hltParticleFlowRecHitHCAL + hltParticleFlowRecHitPS + hltParticleFlowClusterECALForTaus + hltParticleFlowClusterHCAL + hltParticleFlowClusterHFEM + hltParticleFlowClusterHFHAD + hltParticleFlowClusterPS + hltLightPFTracksForTaus + hltParticleFlowBlockForTaus + hltParticleFlowForTaus )
HLTPFJetsSequenceForTaus = cms.Sequence( hltAntiKT5PFJetsForTaus + hltAntiKT5ConvPFJetsForTaus )
HLTPFTauTightIsoSequence = cms.Sequence( hltPFTauJetTracksAssociator + hltPFTauTagInfo + hltPFTausTightIso + hltPFTauTightIsoTrackFindingDiscriminator + hltPFTauTightIsoTrackPt5Discriminator + hltPFTauTightIsoIsolationDiscriminator + hltSelectedPFTausTightIsoTrackFinding + hltSelectedPFTausTightIsoTrackPt5 + hltSelectedPFTausTightIsoTrackFindingIsolation + hltSelectedPFTausTightIsoTrackPt5Isolation + hltConvPFTausTightIsoTrackFinding + hltConvPFTausTightIsoTrackFindingIsolation + hltConvPFTausTightIso + hltConvPFTausTightIsoTrackPt5 + hltConvPFTausTightIsoTrackPt5Isolation )
HLTBTagIPSequenceL25 = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltSelector4Jets + hltBLifetimeL25Jets + hltBLifetimeL25Associator + hltBLifetimeL25TagInfos + hltBLifetimeL25BJetTags )
HLTBTagIPSequenceL3 = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGenerator + hltBLifetimeRegionalCkfTrackCandidates + hltBLifetimeRegionalCtfWithMaterialTracks + hltBLifetimeL3Associator + hltBLifetimeL3TagInfos + hltBLifetimeL3BJetTags )
HLTParticleFlowSequence = cms.Sequence( HLTPreshowerSequence + hltParticleFlowRecHitECAL + hltParticleFlowRecHitHCAL + hltParticleFlowRecHitPS + hltParticleFlowClusterECAL + hltParticleFlowClusterHCAL + hltParticleFlowClusterHFEM + hltParticleFlowClusterHFHAD + hltParticleFlowClusterPS + hltLightPFTracks + hltParticleFlowBlock + hltParticleFlow )
HLTPFJetsSequence = cms.Sequence( hltAntiKT5PFJets + hltAntiKT5ConvPFJets )
HLTBTagIPSequenceL25SlimRA2b = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltGetJetsfromBJetRA2b + hltSelectorJetsRA2b + hltBLifetimeL25JetsRA2b )
HLTBTagIPSequenceL3RA2b = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGeneratorRA2b + hltBLifetimeRegionalCkfTrackCandidatesRA2b + hltBLifetimeRegionalCtfWithMaterialTracksRA2b + hltBLifetimeL3AssociatorRA2b + hltBLifetimeL3TagInfosRA2b + hltBLifetimeL3BJetTagsRA2b )
HLTRSequenceDiJet56 = cms.Sequence( HLTRecoJetSequenceAK5Corrected + hltDiJet56NoID + HLTRecoMETSequence + hltRHemisphere )
HLTBTagIPSequenceL25SlimRAzr = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltGetJetsfromBJetRAzr + hltSelectorJetsRAzr + hltBLifetimeL25JetsRAzr )
HLTBTagIPSequenceL3RAzr = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGeneratorRAzr + hltBLifetimeRegionalCkfTrackCandidatesRAzr + hltBLifetimeRegionalCtfWithMaterialTracksRAzr + hltBLifetimeL3AssociatorRAzr + hltBLifetimeL3TagInfosRAzr + hltBLifetimeL3BJetTagsRAzr )
HLTMuonLocalRecoSequence = cms.Sequence( cms.SequencePlaceholder( "simMuonDTDigis" ) + hltDt1DRecHits + hltDt4DSegments + cms.SequencePlaceholder( "simMuonCSCDigis" ) + hltCsc2DRecHits + hltCscSegments + cms.SequencePlaceholder( "simMuonRPCDigis" ) + hltRpcRecHits )
HLTL2muonrecoNocandSequence = cms.Sequence( HLTMuonLocalRecoSequence + hltL2MuonSeeds + hltL2Muons )
HLTL2muonrecoSequence = cms.Sequence( HLTL2muonrecoNocandSequence + hltL2MuonCandidates )
HLTL3muonTkCandidateSequence = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL3TrajSeedOIState + hltL3TrackCandidateFromL2OIState + hltL3TkTracksFromL2OIState + hltL3MuonsOIState + hltL3TrajSeedOIHit + hltL3TrackCandidateFromL2OIHit + hltL3TkTracksFromL2OIHit + hltL3MuonsOIHit + hltL3TkFromL2OICombination + hltL3TrajSeedIOHit + hltL3TrackCandidateFromL2IOHit + hltL3TkTracksFromL2IOHit + hltL3MuonsIOHit + hltL3TrajectorySeed + hltL3TrackCandidateFromL2 )
HLTL3muonrecoNocandSequence = cms.Sequence( HLTL3muonTkCandidateSequence + hltL3TkTracksFromL2 + hltL3MuonsLinksCombination + hltL3Muons )
HLTL3muonrecoSequence = cms.Sequence( HLTL3muonrecoNocandSequence + hltL3MuonCandidates )
HLTL2muonisorecoSequence = cms.Sequence( hltEcalRawToRecHitFacility + hltEcalRegionalMuonsFEDs + hltEcalRegionalMuonsRecHit + HLTDoLocalHcalSequence + hltTowerMakerForMuons + hltL2MuonIsolations )
HLTL3muonisorecoSequence = cms.Sequence( hltPixelTracks + hltL3MuonIsolations )
HLTL2muonrecoSequenceNoVtx = cms.Sequence( HLTL2muonrecoNocandSequence + hltL2MuonCandidatesNoVtx )
HLTDoRegionalEgammaEcalSequence = cms.Sequence( hltESRawToRecHitFacility + hltEcalRawToRecHitFacility + hltEcalRegionalEgammaFEDs + hltEcalRegionalEgammaRecHit + hltESRegionalEgammaRecHit )
HLTMulti5x5SuperClusterL1Isolated = cms.Sequence( hltMulti5x5BasicClustersL1Isolated + hltMulti5x5SuperClustersL1Isolated + hltMulti5x5EndcapSuperClustersWithPreshowerL1Isolated + hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1Isolated )
HLTL1IsolatedEcalClustersSequence = cms.Sequence( hltHybridSuperClustersL1Isolated + hltCorrectedHybridSuperClustersL1Isolated + HLTMulti5x5SuperClusterL1Isolated )
HLTMulti5x5SuperClusterL1NonIsolated = cms.Sequence( hltMulti5x5BasicClustersL1NonIsolated + hltMulti5x5SuperClustersL1NonIsolated + hltMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated + hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolatedTemp + hltCorrectedMulti5x5EndcapSuperClustersWithPreshowerL1NonIsolated )
HLTL1NonIsolatedEcalClustersSequence = cms.Sequence( hltHybridSuperClustersL1NonIsolated + hltCorrectedHybridSuperClustersL1NonIsolatedTemp + hltCorrectedHybridSuperClustersL1NonIsolated + HLTMulti5x5SuperClusterL1NonIsolated )
HLTDoLocalHcalWithoutHOSequence = cms.Sequence( hltHcalDigis + hltHbhereco + hltHfreco )
HLTPhoton20CaloIdVLIsoLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG20EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG20CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton20CaloIdVLIsoLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton20CaloIdVLIsoLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton20CaloIdVLIsoLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton20CaloIdVLIsoLTrackIsoFilter )
HLTEgammaR9IDSequence = cms.Sequence( hltL1IsoR9ID + hltL1NonIsoR9ID )
HLTPhoton20R9IdPhoton18R9IdSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG20EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton20R9IdPhoton18R9IdEgammaLHEFilter + HLTEgammaR9IDSequence + hltPhoton20R9IdPhoton18R9IdEgammaR9IDFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityR9ID + hltPhoton20R9IdPhoton18R9IdEgammaR9IDDoubleFilter )
HLTPhoton20CaloIdVTIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG20EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltPhoton20CaloIdVTIsoTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton20CaloIdVTIsoTEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton20CaloIdVTIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton20CaloIdVTIsoTHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton20CaloIdVTIsoTTrackIsoFilter )
HLTEle8CaloIdLCaloIsoVLNoL1SeedSequence = cms.Sequence( HLTEcalActivitySequence + hltEG8EtFilterUnseeded + hltActivityPhotonClusterShape + hltEle8CaloIdLCaloIsoVLNoL1SeedClusterShapeFilter + hltActivityPhotonEcalIsol + hltEle8CaloIdLCaloIsoVLNoL1SeedEcalIsolFilter + hltActivityPhotonHcalForHE + hltEle8CaloIdLCaloIsoVLNoL1SeedHEFilter + hltActivityPhotonHcalIsol + hltEle8CaloIdLCaloIsoVLNoL1SeedHcalIsolFilter + hltActivityStartUpElectronPixelSeeds + hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter )
HLTPhoton26Photon18Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton26Photon18EgammaLHEFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HELastFilter )
HLTPhoton26IsoVLPhoton18Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton26IsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton26IsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton26IsoVLTrackIsoFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HELastFilter )
HLTPhoton26IsoVLPhoton18IsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton26IsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton26IsoVLHcalIsoLastFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityPhotonEcalIsol + hltPhoton26IsoVLPhoton18IsoVLEgammaEcalIsolDoubleFilter + hltActivityPhotonHcalIsol + hltPhoton26IsoVLPhoton18IsoVLEgammaHcalIsolDoubleFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTEcalActivityEgammaRegionalRecoTrackerSequence + hltActivityPhotonHollowTrackIsol + hltPhoton26IsoVLPhoton18IsoVLEgammaTrackIsolDoubleFilter )
HLTPhoton26CaloIdLIsoVLPhoton18Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG26CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG26CaloIdLIsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG26CaloIdLIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltEG26CaloIdLIsoVLTrackIsoFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HELastFilter )
HLTPhoton26CaloIdLIsoVLPhoton18R9IdSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG26CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG26CaloIdLIsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG26CaloIdLIsoVLHcalIsoFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityR9ID + hltR9IdFilterUnseeded + hltActivityPhotonClusterShape + hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaClusterShapeFilterUnseeded + hltActivityPhotonEcalIsol + hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaEcalIsolFilterUnseeded + hltActivityPhotonHcalIsol + hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaHcalIsolFilterUnseeded + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTEcalActivityEgammaRegionalRecoTrackerSequence + hltActivityPhotonHollowTrackIsol + hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaTrackIsolFilterUnseeded + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltEG26CaloIdLIsoVLTrackIsoFilter + hltPhoton26CaloIdLIsoVLPhoton18R9IdEgammaDoubleLegCombFilter )
HLTPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG26CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG26CaloIdLIsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG26CaloIdLIsoVLHcalIsoLastFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityPhotonClusterShape + hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaClusterShapeDoubleFilter + hltActivityPhotonEcalIsol + hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaEcalIsolDoubleFilter + hltActivityPhotonHcalIsol + hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaHcalIsolDoubleFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTEcalActivityEgammaRegionalRecoTrackerSequence + hltActivityPhotonHollowTrackIsol + hltPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLEgammaTrackIsolDoubleFilter )
HLTPhoton26R9IdPhoton18CaloIdLIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG26HEFilter + hltL1IsoR9ID + hltL1NonIsoR9ID + hltEG26R9IdFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityR9ID + hltR9IdFilterUnseeded + hltActivityPhotonClusterShape + hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaClusterShapeFilterUnseeded + hltActivityPhotonEcalIsol + hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaEcalIsolFilterUnseeded + hltActivityPhotonHcalIsol + hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaHcalIsolFilterUnseeded + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTEcalActivityEgammaRegionalRecoTrackerSequence + hltActivityPhotonHollowTrackIsol + hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaTrackIsolFilterUnseeded + hltPhoton26R9IdPhoton18CaloIdLIsoVLEgammaDoubleLegCombFilter )
HLTPhoton26R9IdPhoton18R9IdSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG26EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton26R9IdPhoton18R9IdEgammaLHEFilter + HLTEgammaR9IDSequence + hltPhoton26R9IdPhoton18R9IdEgammaR9IDFilter + HLTEcalActivitySequence + hltDoubleIsoEG18EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG18HEFilter + hltActivityR9ID + hltPhoton26R9IdPhoton18R9IdEgammaR9IDDoubleFilter )
HLTPhoton30CaloIdVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG15 + hltEG30EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG30CaloIdVLClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG30CaloIdVLHEFilter )
HLTPhoton30CaloIdVLIsoLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG15 + hltEG30EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG30CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton30CaloIdVLIsoLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton30CaloIdVLIsoLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton30CaloIdVLIsoLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton30CaloIdVLIsoLTrackIsoFilter )
HLTPhoton32CaloIdLPhoton26CaloIdLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG32EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG32CaloIdLHEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG32CaloIdLClusterShapeFilter + HLTEcalActivitySequence + hltDoubleIsoEG26EtFilterUnseeded + hltActivityPhotonHcalForHE + hltPhoton32CaloIdLPhoton26CaloIdLEgammaLHEDoubleFilter + hltActivityPhotonClusterShape + hltPhoton32CaloIdLPhoton26CaloIdLEgammaClusterShapeDoubleFilter )
HLTPhoton36IsoVLPhoton22Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG36EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG36HEFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton36IsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton36IsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton36IsoVLTrackIsoFilter + HLTEcalActivitySequence + hltDoubleIsoEG22EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG22HELastFilter )
HLTPhoton36CaloIdLPhoton22CaloIdLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG36EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG36CaloIdLHEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG36CaloIdLClusterShapeLastFilter + HLTEcalActivitySequence + hltDoubleIsoEG22EtFilterUnseeded + hltActivityPhotonHcalForHE + hltPhoton36CaloIdLPhoton22CaloIdLEgammaLHEDoubleFilter + hltActivityPhotonClusterShape + hltPhoton36CaloIdLPhoton22CaloIdLEgammaClusterShapeDoubleFilter )
HLTPhoton36CaloIdLIsoVLPhoton22Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG36EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG36CaloIdLHEFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG36CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG36CaloIdLIsoVLEcalIsoFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG36CaloIdLIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltEG36CaloIdLIsoVLTrackIsoFilter + HLTEcalActivitySequence + hltDoubleIsoEG22EtFilterUnseeded + hltActivityPhotonHcalForHE + hltDoubleIsoEG22HELastFilter )
HLTPhoton40CaloIdLPhoton28CaloIdLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG40EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG40CaloIdLClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG40CaloIdLHEFilter + HLTEcalActivitySequence + hltDoubleIsoEG28EtFilterUnseeded + hltActivityPhotonHcalForHE + hltPhoton40CaloIdLPhoton28CaloIdLEgammaLHEDoubleFilter + hltActivityPhotonClusterShape + hltPhoton40CaloIdLPhoton28CaloIdLEgammaClusterShapeDoubleFilter )
HLTPhoton50CaloIdVLIsoLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG50EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG50CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton50CaloIdVLIsoLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton50CaloIdVLIsoLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton50CaloIdVLIsoLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton50CaloIdVLIsoLTrackIsoFilter )
HLTSinglePhoton70CaloIdLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG70EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG70CaloIdLClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG70CaloIdLHEFilter )
HLTPhoton75CaloIdVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG75EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG75CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonHcalForHE + HLTDoLocalHcalWithoutHOSequence + hltL1NonIsolatedPhotonHcalForHE + hltPhoton75CaloIdVLHEFilter )
HLTPhoton75CaloIdVLIsoLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG75EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG75CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton75CaloIdVLIsoLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton75CaloIdVLIsoLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton75CaloIdVLIsoLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton75CaloIdVLIsoLTrackIsoFilter )
HLTPhoton90CaloIdVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG90EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG90CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonHcalForHE + HLTDoLocalHcalWithoutHOSequence + hltL1NonIsolatedPhotonHcalForHE + hltPhoton90CaloIdVLHEFilter )
HLTPhoton90CaloIdVLIsoLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG90EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG90CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton90CaloIdVLIsoLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton90CaloIdVLIsoLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton90CaloIdVLIsoLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton90CaloIdVLIsoLTrackIsoFilter )
HLTSinglePhoton125Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG125EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton125HEFilter )
HLTSinglePhoton200NoHESequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG200EtFilter )
HLTDoublePhoton33Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG33EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG33HEFilter + HLTEcalActivitySequence + hltDoubleIsoEG33EtFilterUnseededTight + hltActivityPhotonHcalForHE + hltDoublePhoton33EgammaLHEDoubleFilter )
HLTDoublePhoton33HEVTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG33EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG33HEVTFilter + HLTEcalActivitySequence + hltDoubleEG33EtDoubleFilter + hltActivityPhotonHcalForHE + hltDoubleEG33HEVTDoubleFilter )
HLTPhoton50Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG50EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG50HEFilter )
HLTDoublePhoton50UnseededLegSequence = cms.Sequence( HLTEcalActivitySequence + hltDoubleEG50EtDoubleFilter + hltActivityPhotonHcalForHE + hltDoubleEG50HEDoubleFilter )
HLTPhoton60Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG60EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG60HEFilter )
HLTDoublePhoton60UnseededLegSequence = cms.Sequence( HLTEcalActivitySequence + hltDoubleEG60EtDoubleFilter + hltActivityPhotonHcalForHE + hltDoubleEG60HEDoubleFilter )
HLTDoublePhoton5IsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1DoubleEG2FwdVeto + hltDoublePhoton5IsoVLEtPhiFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltDoublePhoton5IsoVLEgammaHEFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltDoublePhoton5IsoVLEgammaEcalIsolFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltDoublePhoton5IsoVLEgammaHcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltDoublePhoton5IsoVLEgammaTrackIsolFilter )
HLTEle8Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG5 + hltEG8EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle8HEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle8PixelMatchFilter )
HLTEle8CaloIdLCaloIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG5 + hltEG8EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle8CaloIdLCaloIsoVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle8CaloIdLCaloIsoVLEcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle8CaloIdLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle8CaloIdLCaloIsoVLHcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle8CaloIdLCaloIsoVLPixelMatchFilter )
HLTPixelMatchElectronL1IsoTrackingSequence = cms.Sequence( hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso )
HLTPixelMatchElectronL1NonIsoTrackingSequence = cms.Sequence( hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso )
HLTDoElectronDetaDphiSequence = cms.Sequence( hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi )
HLTEle8CaloIdLTrkIdVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG5 + hltEG8EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle8CaloIdLCaloIsoVLClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle8CaloIdLTrkIdVLHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle8CaloIdLTrkIdVLPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltEle8CaloIdLTrkIdVLDetaFilter + hltEle8CaloIdLTrkIdVLDphiFilter )
HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG15EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG15CaloIdTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle15CaloIdTCaloIsoTEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle15CaloIdVTCaloIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle15CaloIdVTCaloIsoTHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle15CaloIdVTCaloIsoTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle15CaloIdVTCaloIsoTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle15CaloIdVTCaloIsoTTrkIdTDetaFilter + hltEle15CaloIdVTCaloIsoTTrkIdTDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter )
HLTEle17CaloIdLCaloIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG17CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG17CaloIdLCaloIsoVLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG17CaloIdLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG17CaloIdLCaloIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdLCaloIsoVLPixelMatchFilter )
HLTEle17CaloIdIsoEle8CaloIdIsoSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG17CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG17CaloIdLCaloIsoVLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG17CaloIdLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG17CaloIdLCaloIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdLCaloIsoVLPixelMatchFilter + HLTEcalActivitySequence + hltDoubleEG8EtFilterUnseeded + hltActivityPhotonClusterShape + hltEle17CaloIdIsoEle8CaloIdIsoClusterShapeDoubleFilter + hltActivityPhotonEcalIsol + hltEle17CaloIdIsoEle8CaloIdIsoEcalIsolDoubleFilter + hltActivityPhotonHcalForHE + hltEle17CaloIdIsoEle8CaloIdIsoHEDoubleFilter + hltActivityPhotonHcalIsol + hltEle17CaloIdIsoEle8CaloIdIsoHcalIsolDoubleFilter + hltActivityStartUpElectronPixelSeeds + hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter )
HLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8ClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8EcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8OneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DetaFilter + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter + HLTEcalActivitySequence + hltDoubleEG8EtFilterUnseeded + hltActivityPhotonHcalForHE + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8HEDoubleFilter + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter )
HLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass30Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8L1MatchFilterRegional + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8ClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8EcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8OneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8DetaFilter + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter + HLTEcalActivitySequence + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8DoubleEtFilter + hltActivityPhotonHcalForHE + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8HEDoubleFilter + hltActivityStartUpElectronPixelSeeds + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PixelMatchDoubleFilter + hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter )
HLTSingleElectronEt17CaloIdIsoSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG17CaloIdLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG17CaloIdLCaloIsoVLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG17CaloIdLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG17CaloIdLCaloIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdLCaloIsoVLPixelMatchFilter )
HLTHFEM15Sequence = cms.Sequence( hltHFEMClusters + hltHFRecoEcalCandidate + hltHFEMFilter )
HLTEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG15 + hltEG27EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG27CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle27CaloIdVLCaloIsoVLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle27CaloIdVLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle27CaloIdVLCaloIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle27CaloIdVLCaloIsoVLPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle27CaloIdVLCaloIsoVLOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle27CaloIdVLCaloIsoVLTrkIdVLDetaFilter + hltEle27CaloIdVLCaloIsoVLTrkIdVLDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter )
HLTEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG15 + hltEG27EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG27CaloIdTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle27CaloIdTCaloIsoTEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle27CaloIdVTCaloIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle27CaloIdVTCaloIsoTHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle27CaloIdVTCaloIsoTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle27CaloIdVTCaloIsoTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle27CaloIdVTCaloIsoTTrkIdTDetaFilter + hltEle27CaloIdVTCaloIsoTTrkIdTDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter )
HLTEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG32EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG32CaloIdVLClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle32CaloIdVLCaloIsoVLEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle32CaloIdVLCaloIsoVLHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle32CaloIdVLCaloIsoVLHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle32CaloIdVLCaloIsoVLPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle32CaloIdVLCaloIsoVLOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle32CaloIdVLCaloIsoVLTrkIdVLDetaFilter + hltEle32CaloIdVLCaloIsoVLTrkIdVLDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter )
HLTEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG32EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG32CaloIdTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle32CaloIdTCaloIsoTEcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle32CaloIdVTCaloIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle32CaloIdVTCaloIsoTHcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle32CaloIdVTCaloIsoTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle32CaloIdVTCaloIsoTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle32CaloIdVTCaloIsoTTrkIdTDetaFilter + hltEle32CaloIdVTCaloIsoTTrkIdTDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter )
HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17L1MatchFilterRegional + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17ClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17EcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17OneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DetaFilter + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter + HLTEcalActivitySequence + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17DoubleEtFilter + hltActivityPhotonHcalForHE + hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter )
HLTEle45CaloIdVTTrkIdTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG45EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG45CaloIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG45CaloIdVTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle45CaloIdVTPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle45CaloIdVTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle45CaloIdVTTrkIdTDetaFilter + hltEle45CaloIdVTTrkIdTDphiFilter )
HLTEle52CaloIdVTTrkIdTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG52EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG52CaloIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG52CaloIdVTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle52CaloIdVTPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle52CaloIdVTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle52CaloIdVTTrkIdTDetaFilter + hltEle52CaloIdVTTrkIdTDphiFilter )
HLTDoEGammaStartupSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate )
HLTDoEgammaClusterShapeSequence = cms.Sequence( hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape )
HLTDoEGammaHESequence = cms.Sequence( HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE )
HLTDoEGammaPixelSequence = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds )
HLTDoubleEle8L1NonIsoHLTCaloIdLSequence = cms.Sequence( HLTDoEGammaStartupSequence + hltEGRegionalL1DoubleEG5 + hltDoubleEG8EtFilterL1DoubleEG5 + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLTCaloIdLDoubleEle8ClusterShapeFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTCaloIdLDoubleEle8HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTCaloIdLDoubleEle8PixelMatchFilter )
HLTPhoton33Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG33EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG33HEFilter )
HLTDoublePhoton33UnseededLegSequence = cms.Sequence( HLTEcalActivitySequence + hltDoubleEG33EtDoubleFilter + hltActivityPhotonHcalForHE + hltDoubleEG33HEDoubleFilter )
HLTActivityPixelMatchSequence = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltActivityStartUpElectronPixelSeeds )
HLTCaloTausCreatorRegionalSequence = cms.Sequence( HLTDoRegionalJetEcalSequence + HLTDoLocalHcalSequence + hltTowerMakerForJets + hltCaloTowersTau1Regional + hltIconeTau1Regional + hltCaloTowersTau2Regional + hltIconeTau2Regional + hltCaloTowersTau3Regional + hltIconeTau3Regional + hltCaloTowersTau4Regional + hltIconeTau4Regional + hltCaloTowersCentral1Regional + hltIconeCentral1Regional + hltCaloTowersCentral2Regional + hltIconeCentral2Regional + hltCaloTowersCentral3Regional + hltIconeCentral3Regional + hltCaloTowersCentral4Regional + hltIconeCentral4Regional )
HLTL2TauJetsSequence = cms.Sequence( HLTCaloTausCreatorRegionalSequence + hltL2TauJets )
HLTBTagMuDiJet20SequenceL25 = cms.Sequence( HLTL2muonrecoNocandSequence + hltBSoftMuonGetJetsFromDiJet20 + hltSelector4JetsDiJet20 + hltBSoftMuonDiJet20L25Jets + hltBSoftMuonDiJet20L25TagInfos + hltBSoftMuonDiJet20L25BJetTagsByDR )
HLTBTagMuDiJet20Mu5SelSequenceL3 = cms.Sequence( HLTL3muonrecoNocandSequence + hltBSoftMuonMu5L3 + hltBSoftMuonDiJet20Mu5SelL3TagInfos + hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR )
HLTBTagMuDiJet40SequenceL25 = cms.Sequence( HLTL2muonrecoNocandSequence + hltBSoftMuonGetJetsFromDiJet40 + hltSelector4JetsDiJet40 + hltBSoftMuonDiJet40L25Jets + hltBSoftMuonDiJet40L25TagInfos + hltBSoftMuonDiJet40L25BJetTagsByDR )
HLTBTagMuDiJet40Mu5SelSequenceL3 = cms.Sequence( HLTL3muonrecoNocandSequence + hltBSoftMuonMu5L3 + hltBSoftMuonDiJet40Mu5SelL3TagInfos + hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR )
HLTBTagMuDiJet70SequenceL25 = cms.Sequence( HLTL2muonrecoNocandSequence + hltBSoftMuonGetJetsFromDiJet70 + hltSelector4JetsDiJet70 + hltBSoftMuonDiJet70L25Jets + hltBSoftMuonDiJet70L25TagInfos + hltBSoftMuonDiJet70L25BJetTagsByDR )
HLTBTagMuDiJet70Mu5SelSequenceL3 = cms.Sequence( HLTL3muonrecoNocandSequence + hltBSoftMuonMu5L3 + hltBSoftMuonDiJet70Mu5SelL3TagInfos + hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR )
HLTBTagMuDiJet110SequenceL25 = cms.Sequence( HLTL2muonrecoNocandSequence + hltBSoftMuonGetJetsFromDiJet110 + hltSelector4JetsDiJet110 + hltBSoftMuonDiJet110L25Jets + hltBSoftMuonDiJet110L25TagInfos + hltBSoftMuonDiJet110L25BJetTagsByDR )
HLTBTagMuDiJet110Mu5SelSequenceL3 = cms.Sequence( HLTL3muonrecoNocandSequence + hltBSoftMuonMu5L3 + hltBSoftMuonDiJet110Mu5SelL3TagInfos + hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR )
HLTMu5Ele8CaloIdLTrkIdVLEle8L1NonIsoHLTnonIsoSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1MuOpenEG5 + hltEG8EtFilterMuOpenEG5 + HLTDoEgammaClusterShapeSequence + hltEG8CaloIdLClusterShapeFilterMuOpenEG5 + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG8CaloIdLHEFilterMuOpenEG5 + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEG8CaloIdLPixelMatchFilterMuOpenEG5 + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle8CaloIdLTrkIdVLOneOEMinusOneOPFilterRegionalMuOpenEG5 + HLTDoElectronDetaDphiSequence + hltEle8CaloIdLTrkIdVLDetaFilterMuOpenEG5 + hltEle8CaloIdLTrkIdVLDphiFilterMuOpenEG5 + HLTEcalActivitySequence + hltDoubleEG8EtFilterUnseeded + hltActivityPhotonHcalForHE + hltL1NonIsoHLTNonIsoMu5DoubleEle8NoCandHEFilter + hltActivityStartUpElectronPixelSeeds + hltDoubleEG8NoCandPixelMatchFilter )
HLTRSequenceNoJetFilter = cms.Sequence( HLTRecoJetSequenceAK5Corrected + HLTRecoMETSequence + hltRHemisphere )
HLTPhoton20CaloIdVTIsoTMu8Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1MuOpenEG5 + hltEG20EtFilterMuOpenEG5 + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltPhoton20CaloIdVTIsoTMu8ClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltPhoton20CaloIdVTIsoTMu8EcalIsoFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltPhoton20CaloIdVTIsoTMu8HEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltPhoton20CaloIdVTIsoTMu8HcalIsoFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsolatedPhotonHollowTrackIsol + hltL1NonIsolatedPhotonHollowTrackIsol + hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter )
HLTPFJetTriggerSequenceForTaus = cms.Sequence( HLTTrackReconstructionForJets + HLTParticleFlowSequenceForTaus + HLTPFJetsSequenceForTaus )
HLTPFTauSequence = cms.Sequence( hltPFTauJetTracksAssociator + hltPFTauTagInfo + hltPFTaus + hltPFTauTrackFindingDiscriminator + hltPFTauLooseIsolationDiscriminator + hltSelectedPFTausTrackFinding + hltSelectedPFTausTrackFindingLooseIsolation + hltConvPFTausTrackFinding + hltConvPFTausTrackFindingLooseIsolation + hltConvPFTaus )
HLTBTagIPSequenceL25SingleTop = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltGetJetsfromBJet30Central + hltSelectorJetsSingleTop + hltBLifetimeL25JetsSingleTop + hltBLifetimeL25AssociatorSingleTop + hltBLifetimeL25TagInfosSingleTop + hltBLifetimeL25BJetTagsSingleTop )
HLTBTagIPSequenceL3SingleTop = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGeneratorSingleTop + hltBLifetimeRegionalCkfTrackCandidatesSingleTop + hltBLifetimeRegionalCtfWithMaterialTracksSingleTop + hltBLifetimeL3AssociatorSingleTop + hltBLifetimeL3TagInfosSingleTop + hltBLifetimeL3BJetTagsSingleTop )
HLTDoubleMu5Ele8L1NonIsoHLTnonIsoSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1MuOpenEG5 + hltEG8EtFilterMuOpenEG5 + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG8HEFilterMuOpenEG5 + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEG8PixelMatchFilterMuOpenEG5 )
HLTSinglePhoton40CaloIdLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG40EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG40CaloIdLClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG40CaloIdLHEFilter )
HLTDoublePhoton40Sequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG20 + hltEG40EtFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG40HEFilter + HLTEcalActivitySequence + hltDoubleIsoEG40EtFilterUnseededTight + hltActivityPhotonHcalForHE + hltDoublePhoton40EgammaLHEDoubleFilter )
HLTEle8CaloIdTTrkIdTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG5 + hltEG8EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle8CaloIdTTrkIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle8CaloIdTTrkIdTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle8CaloIdTTrkIdTPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle8CaloIdTTrkIdTOneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltEle8CaloIdTTrkIdTDetaFilter + hltEle8CaloIdTTrkIdTDphiFilter )
HLTEle15L1EG5HTT75CaloIdTCaloIsoVLTrkIdTTrkIsoVLSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1EG5HTT75 + hltEG15EtFilterL1EG5HTT75 + HLTDoEgammaClusterShapeSequence + hltEG15CaloIdTClusterShapeFilterEG5HTT75 + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEG15CaloIdTCaloIsoVLEcalIsoFilterEG5HTT75 + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEG15CaloIdTCaloIsoVLHEFilterEG5HTT75 + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEG15CaloIdTCaloIsoVLHcalIsoFilterEG5HTT75 + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEG15CaloIdTCaloIsoVLPixelMatchFilterEG5HTT75 + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle15CaloIdTCaloIsoVLOneOEMinusOneOPFilterEG5HTT75 + HLTDoElectronDetaDphiSequence + hltEle15CaloIdTCaloIsoVLTrkIdTDetaFilterEG5HTT75 + hltEle15CaloIdTCaloIsoVLTrkIdTDphiFilterEG5HTT75 + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLTrackIsoFilterEG5HTT75 )
HLTEle15CaloIdVTTrkIdTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG15EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG15CaloIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle15CaloIdVTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle15CaloIdVTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle15CaloIdVTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle15CaloIdVTTrkIdTDetaFilter + hltEle15CaloIdVTTrkIdTDphiFilter )
HLTEle25CaloIdVTCaloTrkIdSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG25EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEle25CaloIdVTTrkIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle25CaloIdVTTrkIdTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle25CaloIdVTTrkIdTPixelMatchFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltEle25CaloIdVTTrkIdTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle25CaloIdVTTrkIdTDetaFilter + hltEle25CaloIdVTTrkIdTDphiFilter )
HLTBTagIPSequenceL25EleJetSingleTop = cms.Sequence( HLTDoLocalPixelSequence + HLTRecopixelvertexingSequence + hltGetJetsfrom1EleCleanBJet30Central + hltSelectorEleJetsSingleTop + hltBLifetimeL25JetsEleJetSingleTop + hltBLifetimeL25AssociatorEleJetSingleTop + hltBLifetimeL25TagInfosEleJetSingleTop + hltBLifetimeL25BJetTagsEleJetSingleTop )
HLTBTagIPSequenceL3EleJetSingleTop = cms.Sequence( HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltBLifetimeRegionalPixelSeedGeneratorEleJetSingleTop + hltBLifetimeRegionalCkfTrackCandidatesEleJetSingleTop + hltBLifetimeRegionalCtfWithMaterialTracksEleJetSingleTop + hltBLifetimeL3AssociatorEleJetSingleTop + hltBLifetimeL3TagInfosEleJetSingleTop + hltBLifetimeL3BJetTagsEleJetSingleTop )
HLTEle17CaloIdVTTrkIdTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG17CaloIdTClusterShapeFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle17CaloIdVTTrkIdTHEFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdVTTrkIdTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle17CaloIdVTTrkIdTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle17CaloIdVTTrkIdTDetaFilter + hltEle17CaloIdVTTrkIdTDphiFilter )
HLTEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG17EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG17CaloIdTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle17CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter )
HLTPFJetTriggerSequence = cms.Sequence( HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence )
HLTEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence = cms.Sequence( HLTDoRegionalEgammaEcalSequence + HLTL1IsolatedEcalClustersSequence + HLTL1NonIsolatedEcalClustersSequence + hltL1IsoRecoEcalCandidate + hltL1NonIsoRecoEcalCandidate + hltEGRegionalL1SingleEG12 + hltEG25EtFilter + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG25CaloIdTClusterShapeFilter + hltL1IsolatedPhotonEcalIsol + hltL1NonIsolatedPhotonEcalIsol + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTEcalIsolFilter + HLTDoLocalHcalWithoutHOSequence + hltL1IsolatedPhotonHcalForHE + hltL1NonIsolatedPhotonHcalForHE + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHEFilter + hltL1IsolatedPhotonHcalIsol + hltL1NonIsolatedPhotonHcalIsol + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTHcalIsolFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTPixelMatchFilter + hltCkfL1IsoTrackCandidates + hltCtfL1IsoWithMaterialTracks + hltPixelMatchElectronsL1Iso + hltCkfL1NonIsoTrackCandidates + hltCtfL1NonIsoWithMaterialTracks + hltPixelMatchElectronsL1NonIso + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTOneOEMinusOneOPFilter + hltElectronL1IsoDetaDphi + hltElectronL1NonIsoDetaDphi + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDetaFilter + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter + HLTL1IsoEgammaRegionalRecoTrackerSequence + HLTL1NonIsoEgammaRegionalRecoTrackerSequence + hltL1IsoElectronTrackIsol + hltL1NonIsoElectronTrackIsol + hltEle25CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter )
HLTDoubleEle8HTT50L1NonIsoHLTCaloIdLSequence = cms.Sequence( HLTDoEGammaStartupSequence + hltEGRegionalL1DoubleEG5HTT50 + hltDoubleEG8EtFilterL1DoubleEG5HTT50 + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLTCaloIdLDoubleEle8HTT50ClusterShapeFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTCaloIdLDoubleEle8HTT50HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTCaloIdLDoubleEle8HTT50PixelMatchFilter )
HLTDoubleEle8HTT50L1NonIsoHLTCaloIdTSequence = cms.Sequence( HLTDoEGammaStartupSequence + hltEGRegionalL1DoubleEG5HTT50 + hltDoubleEG8EtFilterL1DoubleEG5HTT50 + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLTCaloIdTDoubleEle8HTT50ClusterShapeFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTCaloIdTDoubleEle8HTT50HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTCaloIdTDoubleEle8HTT50PixelMatchFilter )
HLTTripleElectronEt10L1NonIsoHLTNonIsoSequence = cms.Sequence( HLTDoEGammaStartupSequence + hltEGRegionalL1TripleEG5 + hltTripleEG10EtFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTNonIsoTripleElectronEt10HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTNonIsoTripleElectronEt10PixelMatchFilter )
HLTRecopixelvertexingForHighMultSequence = cms.Sequence( hltPixelTracksForHighMult + hltPixelVerticesForHighMult )
HLTDoLocalPixelLight = cms.Sequence( hltSiPixelDigis + hltSiPixelClusters )

HLTriggerFirstPath = cms.Path( hltGetRaw + hltBoolFalse )
HLT_Activity_Ecal_SC7_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1BscMinBiasORBptxPlusANDMinus + hltPreActivityEcalSC7 + HLTEcalActivitySequence + hltEgammaSelectEcalSuperClustersActivityFilterSC7 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleJet16_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet16 + hltPreL1SingleJet16 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleJet36_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet36 + hltPreL1SingleJet36 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet30_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet16 + hltPreJet30 + HLTRecoJetSequenceAK5Corrected + hltSingleJet30 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet60_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet36 + hltPreJet60 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet60Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet80_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet52 + hltPreJet80 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet80Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet110_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet68 + hltPreJet110 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet110Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet150_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet150 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet150Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet190_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet190 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet190Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet240_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet240 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet240Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet300_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet300 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet300Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet370_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet370 + HLTRegionalRecoJetSequenceAK5Corrected + hltSingleJet370Regional + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Jet370_NoJetID_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreJet370NoJetID + HLTRegionalTowerMakerForJetsSequence + hltAntiKT5CaloJetsRegional + hltAntiKT5L2L3CorrCaloJetsRegional + hltL1MatchedJetsRegional + hltSingleJet370RegionalNoJetID + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve30_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet16 + hltPreDiJetAve30 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve30 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve60_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet36 + hltPreDiJetAve60 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve60 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve80_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet52 + hltPreDiJetAve80 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve80 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve110_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet68 + hltPreDiJetAve110 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve110 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve150_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJetAve150 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve190_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJetAve190 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve190 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve240_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJetAve240 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve240 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve300_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJetAve300 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve300 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJetAve370_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJetAve370 + HLTRecoJetSequenceAK5Corrected + hltDiJetAve370 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleJet30_ForwardBackward_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleForJet32EtaOpp + hltPreDoubleJet30ForwardBackward + HLTRecoJetSequenceAK5Corrected + hltDoubleJet30ForwardBackward + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleJet60_ForwardBackward_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleForJet32EtaOpp + hltPreDoubleJet60ForwardBackward + HLTRecoJetSequenceAK5Corrected + hltDoubleJet60ForwardBackward + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleJet70_ForwardBackward_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleForJet32EtaOpp + hltPreDoubleJet70ForwardBackward + HLTRecoJetSequenceAK5Corrected + hltDoubleJet70ForwardBackward + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleJet80_ForwardBackward_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleForJet44EtaOpp + hltPreDoubleJet80ForwardBackward + HLTRecoJetSequenceAK5Corrected + hltDoubleJet80ForwardBackward + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJet130_PT130_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet68 + hltPreDiJet130PT130 + HLTRecoJetSequenceAK5Corrected + hltDijet130PT130 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJet160_PT160_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet92 + hltPreDiJet160PT160 + HLTRecoJetSequenceAK5Corrected + hltDijet160PT160 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_CentralJet80_MET65_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreCenJet80MET65 + HLTRegionalRecoJetSequenceAK5Corrected + hltCenJet80CentralRegional + HLTRecoMETSequence + hltMET65 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_CentralJet80_MET80HF_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreCenJet80MET80HF + HLTRegionalRecoJetSequenceAK5Corrected + hltJetIDMPassedJetsRegional + hltCenJet80MCentralRegional + HLTRecoMETSequence + hltMET80 + hltMetWithHF + hltMETWithHF80 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_CentralJet80_MET100_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreCenJet80MET100 + HLTRegionalRecoJetSequenceAK5Corrected + hltCenJet80CentralRegional + HLTRecoMETSequence + hltMET100 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_CentralJet80_MET160_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreCenJet80MET160 + HLTRegionalRecoJetSequenceAK5Corrected + hltCenJet80CentralRegional + HLTRecoMETSequence + hltMET160 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiJet60_MET45_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM20 + hltPreDiJet60MET45 + HLTRecoJetSequenceAK5Corrected + hltDiJet60 + HLTRecoMETSequence + hltMet45 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiCentralJet20_MET80_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPre2CenJet20MET80 + HLTRegionalRecoJetSequenceAK5Corrected + hlt2CenJet20CentralRegional + HLTRecoMETSequence + hltMET80 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DiCentralJet20_BTagIP_MET65_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPre2CenJet20BtagIPMET65 + HLTRecoMETSequence + hltMET65 + HLTRecoJetSequenceAK5Corrected + hltBJetHbb + HLTBtagIPSequenceL25Hbb + hltBLifetimeL25FilterHbb + HLTBtagIPSequenceL3Hbb + hltBLifetimeL3FilterHbb + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet40_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet40 + HLTRecoJetSequenceAK5Corrected + hltQuadJet40Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet40_IsoPFTau40_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet40IsoPFTau40 + HLTRecoJetSequenceAK5Corrected + hltQuadJet40IsoPFTau40 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequenceForTaus + HLTPFJetsSequenceForTaus + HLTPFTauTightIsoSequence + hltPFTau5Track + hltPFTau5Track5 + hltFilterPFTauTrack5TightIsoL1QuadJet20Central + hltFilterPFTauTrack5TightIsoL1QuadJet20CentralPFTau40 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet50_BTagIP_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet50 + HLTRecoJetSequenceAK5Corrected + hltQuadJet50Central + HLTBTagIPSequenceL25 + hltBLifetimeL25Filter + HLTBTagIPSequenceL3 + hltBLifetimeL3Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet50_Jet40_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet50Jet40 + HLTRecoJetSequenceAK5Corrected + hltPentaJet40Central + hltQuadJet50Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet60_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet60 + HLTRecoJetSequenceAK5Corrected + hltQuadJet60 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_QuadJet70_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1QuadJet20Central + hltPreQuadJet70 + HLTRecoJetSequenceAK5Corrected + hltQuadJet70 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_ExclDiJet60_HFOR_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet36 + hltPreExclDiJet60HFOR + HLTRecoJetSequenceAK5Corrected + hltExclDiJet60HFOR + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_ExclDiJet60_HFAND_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleJet36FwdVeto + hltPreExclDiJet60HFAND + HLTRecoJetSequenceAK5Corrected + hltExclDiJet60HFAND + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT150_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT50 + hltPreHT150 + HLTRecoJetSequenceAK5Corrected + hltHT150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT150_AlphaT0p60_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT75 + hltPreHLTHT150AlphaT0p6 + HLTRecoJetSequenceAK5Corrected + hltHT150AlphaT0p6 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT200_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT75 + hltPreHT200 + HLTRecoJetSequenceAK5Corrected + hltHT200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT200_AlphaT0p53_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT200AlphaT0p53 + HLTRecoJetSequenceAK5Corrected + hltHT200AlphaT0p53 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT200_AlphaT0p60_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT75 + hltPreHLTHT200AlphaT0p6 + HLTRecoJetSequenceAK5Corrected + hltHT200AlphaT0p6 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT250 + HLTRecoJetSequenceAK5Corrected + hltHT250 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_AlphaT0p53_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT250AlphaT0p53 + HLTRecoJetSequenceAK5Corrected + hltHT250AlphaT0p53 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_AlphaT0p54_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT250AlphaT0p54 + HLTRecoJetSequenceAK5Corrected + hltHT250AlphaT0p54 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_MHT60_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT250MHT60 + HLTRecoJetSequenceAK5Corrected + hltHT250 + hltMHT60 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_MHT70_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT250MHT70 + HLTRecoJetSequenceAK5Corrected + hltHT250 + hltMHT70 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300 + HLTRecoJetSequenceAK5Corrected + hltHT300 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_MHT75_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300MHT75 + HLTRecoJetSequenceAK5Corrected + hltHT300 + hltMHT75 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_PFMHT55_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300PFMHT55 + HLTRecoJetSequenceAK5Corrected + hltHT300 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT55Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_CentralJet30_BTagIP_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300BTagIP + HLTRecoJetSequenceAK5Corrected + hltHT300 + hltBJetRA2b + HLTBTagIPSequenceL25SlimRA2b + HLTBTagIPSequenceL3RA2b + hltBLifetimeL3FilterRA2b + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_CentralJet30_BTagIP_PFMHT55_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300BTagIPPFMHT55 + HLTRecoJetSequenceAK5Corrected + hltHT300 + hltBJetRA2b + HLTBTagIPSequenceL25SlimRA2b + HLTBTagIPSequenceL3RA2b + hltBLifetimeL3FilterRA2b + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT55Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_CentralJet30_BTagIP_PFMHT75_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300BTagIPPFMHT75 + HLTRecoJetSequenceAK5Corrected + hltHT300 + hltBJetRA2b + HLTBTagIPSequenceL25SlimRA2b + HLTBTagIPSequenceL3RA2b + hltBLifetimeL3FilterRA2b + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT75Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_AlphaT0p52_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT300AlphaT0p52 + HLTRecoJetSequenceAK5Corrected + hltHT300AlphaT0p52 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_AlphaT0p53_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT300AlphaT0p53 + HLTRecoJetSequenceAK5Corrected + hltHT300AlphaT0p53 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT350_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT350 + HLTRecoJetSequenceAK5Corrected + hltHT350 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT350_AlphaT0p51_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT350AlphaT0p51 + HLTRecoJetSequenceAK5Corrected + hltHT350AlphaT0p51 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT350_AlphaT0p53_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT350AlphaT0p53 + HLTRecoJetSequenceAK5Corrected + hltHT350AlphaT0p53 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT400_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT400 + HLTRecoJetSequenceAK5Corrected + hltHT400 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT400_AlphaT0p51_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHLTHT400AlphaT0p51 + HLTRecoJetSequenceAK5Corrected + hltHT400AlphaT0p51 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT450_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT450 + HLTRecoJetSequenceAK5Corrected + hltHT450 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT500_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT500 + HLTRecoJetSequenceAK5Corrected + hltHT500 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT550_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT550 + HLTRecoJetSequenceAK5Corrected + hltHT550 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_PFMHT150_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPrePFMHT150 + HLTRecoMETSequence + hltMET80 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT150Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_MET100_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreMET100 + HLTRecoMETSequence + hltMET100 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_MET120_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreMET120 + HLTRecoMETSequence + hltMET120 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_MET200_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreMET200 + HLTRecoMETSequence + hltMET200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R014_MR150_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR014MR150 + HLTRSequenceDiJet56 + hltR014MR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R014_MR150_CentralJet40_BTagIP_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR014MR150BTag + HLTRSequenceDiJet56 + hltR014MR150 + hltBJetRAzr + HLTBTagIPSequenceL25SlimRAzr + HLTBTagIPSequenceL3RAzr + hltBLifetimeL3FilterRAzr + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R014_MR450_CentralJet40_BTagIP_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR014MR450BTag + HLTRSequenceDiJet56 + hltR014MR450 + hltBJetRAzr + HLTBTagIPSequenceL25SlimRAzr + HLTBTagIPSequenceL3RAzr + hltBLifetimeL3FilterRAzr + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R014_MR600_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR014MR600 + HLTRSequenceDiJet56 + hltR014MR600 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R020_MR150_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR020MR150 + HLTRSequenceDiJet56 + hltR020MR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R020_MR350_CentralJet40_BTagIP_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR020MR350BTag + HLTRSequenceDiJet56 + hltR020MR350 + hltBJetRAzr + HLTBTagIPSequenceL25SlimRAzr + HLTBTagIPSequenceL3RAzr + hltBLifetimeL3FilterRAzr + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R020_MR500_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR020MR500 + HLTRSequenceDiJet56 + hltR020MR500 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R020_MR550_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR020MR550 + HLTRSequenceDiJet56 + hltR020MR550 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R025_MR150_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR025MR150 + HLTRSequenceDiJet56 + hltR025MR150 )
HLT_R025_MR250_CentralJet40_BTagIP_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR025MR250BTag + HLTRSequenceDiJet56 + hltR025MR250 + hltBJetRAzr + HLTBTagIPSequenceL25SlimRAzr + HLTBTagIPSequenceL3RAzr + hltBLifetimeL3FilterRAzr + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R025_MR400_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR025MR400 + HLTRSequenceDiJet56 + hltR025MR400 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R025_MR450_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR025MR450 + HLTRSequenceDiJet56 + hltR025MR450 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R033_MR300_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR033MR300 + HLTRSequenceDiJet56 + hltR033MR300 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R033_MR350_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR033MR350 + HLTRSequenceDiJet56 + hltR033MR350 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R038_MR200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR038MR200 + HLTRSequenceDiJet56 + hltR038MR200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_R038_MR250_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreR038MR250 + HLTRSequenceDiJet56 + hltR038MR250 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleMuOpen_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreL1SingleMuOpen + hltL1MuOpenL1Filtered0 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleMuOpen_DT_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreL1SingleMuOpenDT + hltL1MuOpenL1FilteredDT + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleMu10_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreL1Mu10 + hltL1SingleMu10L1Filtered0 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleMu20_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu20 + hltPreL1Mu20 + hltL1SingleMu20L1Filtered0 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1DoubleMu0_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreL1DoubleMu0 + hltDiMuonL1Filtered0 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L2Mu10_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreL2Mu10 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L2Mu20_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreL2Mu20 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu20L2Filtered20 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L2DoubleMu0_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreL2DoubleMu0 + hltDiMuonL1Filtered0 + HLTL2muonrecoSequence + hltDiMuonL2PreFiltered0 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu3_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreMu3 + hltSingleMuOpenL1Filtered + HLTL2muonrecoSequence + hltSingleMu3L2Filtered0 + HLTL3muonrecoSequence + hltSingleMu3L3Filtered3 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu5_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu3 + hltPreMu5 + hltL1SingleMu3L1Filtered0 + HLTL2muonrecoSequence + hltSingleMu5L2Filtered3 + HLTL3muonrecoSequence + hltSingleMu5L3Filtered5 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu3 + hltPreMu8 + hltL1SingleMu3L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu3L2Filtered3 + HLTL3muonrecoSequence + hltSingleMu8L3Filtered8 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu12_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu12 + hltL1SingleMu7L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu7L2Filtered7 + HLTL3muonrecoSequence + hltSingleMu12L3Filtered12 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu15_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreMu15 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + HLTL3muonrecoSequence + hltSingleMu15L3Filtered15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu20_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreMu20 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltSingleMu12L2Filtered12 + HLTL3muonrecoSequence + hltSingleMu20L3Filtered20 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu24_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreMu24 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu12L2Filtered12 + HLTL3muonrecoSequence + hltSingleMu24L3Filtered24 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu30_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreMu30 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu12L2Filtered12 + HLTL3muonrecoSequence + hltSingleMu30L3Filtered30 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu40_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu16 + hltPreMu40 + hltL1SingleMu16L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu16L2Filtered16 + HLTL3muonrecoSequence + hltSingleMu40L3Filtered40 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu12_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreIsoMu12 + hltL1SingleMu7L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu7L2Filtered7 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered7 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered12 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered12 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu15_v7 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreIsoMu15 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered10 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered15 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu17_v7 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreIsoMu17 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered10 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered17 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered17 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu24_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreIsoMu24 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu12L2Filtered12 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered12 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered24 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered24 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu30_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu12 + hltPreIsoMu30 + hltL1SingleMu12L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu12L2Filtered12 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered12 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered30 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered30 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L2DoubleMu23_NoVertex_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreL2DoubleMu23NoVertex + hltL1DoubleMuon3L1Filtered0 + HLTL2muonrecoSequenceNoVtx + hltL2DoubleMu23NoVertexL2PreFiltered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu3_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDoubleMu3 + hltDiMuonL1Filtered0 + HLTL2muonrecoSequence + hltDiMuonL2PreFiltered0 + HLTL3muonrecoSequence + hltDiMuonL3PreFiltered3 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu6_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreDoubleMu6 + hltL1DoubleMuon3L1Filtered0 + HLTL2muonrecoSequence + hltDiMuon3L2PreFiltered0 + HLTL3muonrecoSequence + hltDiMuonL3PreFiltered6 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu7_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreDoubleMu7 + hltL1DoubleMuon3L1Filtered0 + HLTL2muonrecoSequence + hltDiMuon3L2PreFiltered0 + HLTL3muonrecoSequence + hltDiMuonL3PreFiltered7 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu2_Bs_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDoubleMu0Bs + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltDimuonL3PreFiltered2Bs + hltDoubleMu2BsL3Filtered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu4_Acoplanarity03_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreDoubleMu4Excl + hltL1DoubleMuon3L1Filtered3 + HLTL2muonrecoSequence + hltL2DoubleMu3L2Filtered + HLTL3muonrecoSequence + hltDiMuonL3PreFiltered4 + hltDoubleMu4ExclL3PreFiltered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon0_Jpsi_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon0Jpsi + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltJpsiL3Filtered + hltDisplacedmumuVtxProducerJpsi0 + hltVertexmumuFilterJpsi + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon0_Upsilon_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon0Upsilon + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltUpsilonL3Filtered + hltDisplacedmumuVtxProducerUpsilon + hltVertexmumuFilterUpsilon + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon4_Bs_Barrel_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDoubleMu2BarrelBs + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltDoubleMu2BarrelBsL3Filtered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon5_Upsilon_Barrel_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon5BarrelUpsilon + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltBarrelUpsilonL3Filtered + hltDisplacedmumuVtxProducerUpsilonBarrel + hltVertexmumuFilterUpsilonBarrel + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon6_Bs_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDoubleMu2Dimuon6Bs + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltDoubleMu2Dimuon6BsL3Filtered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon7_LowMass_Displaced_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon7LowMassDisplaced + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltLowMassDisplacedL3Filtered + hltDisplacedmumuVtxProducerLowMass + hltDisplacedmumuFilterLowMass + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon7_Jpsi_Displaced_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon7JpsiDisplaced + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltJpsiDisplacedL3Filtered + hltDisplacedmumuVtxProducerJpsi + hltDisplacedmumuFilterJpsi + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon7_Jpsi_X_Barrel_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon7JpsiXBarrel + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltJpsiXBarrelL3Filtered + hltDisplacedmumuVtxProducerJpsiXBarrel + hltVertexmumuFilterJpsiXBarrel + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon7_PsiPrime_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon7PsiPrime + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltPsiPrimeL3Filtered + hltDisplacedmumuVtxProducerPsiPrime + hltVertexmumuFilterPsiPrime + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon10_Jpsi_Barrel_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon10BarrelJpsi + hltDimuonL1Filtered0 + HLTL2muonrecoSequence + hltDimuonL2PreFiltered0 + HLTL3muonrecoSequence + hltBarrelJpsiL3Filtered + hltDisplacedmumuVtxProducerJpsiBarrel + hltVertexmumuFilterJpsiBarrel + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon0_Jpsi_Muon_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon0JpsiMuon + hltTripleMuonL1Filtered0 + HLTL2muonrecoSequence + hltTripleMuonL2PreFiltered0 + HLTL3muonrecoSequence + hltTripleMuL3PreFiltered0 + hltJpsiMuonL3Filtered + hltDisplacedmumuVtxProducerJpsiMuon + hltVertexmumuFilterJpsiMuon + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Dimuon0_Upsilon_Muon_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreDimuon0UpsilonMuon + hltTripleMuonL1Filtered0 + HLTL2muonrecoSequence + hltTripleMuonL2PreFiltered0 + HLTL3muonrecoSequence + hltTripleMuL3PreFiltered0 + hltUpsilonMuonL3Filtered + hltDisplacedmumuVtxProducerUpsilonMuon + hltVertexmumuFilterUpsilonMuon + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu13_Mu8_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreMu13Mu8 + hltL1DoubleMuon3L1Filtered0 + HLTL2muonrecoSequence + hltDiMuon3L2PreFiltered0 + hltL1DoubleMuon3L2Filtered7 + HLTL3muonrecoSequence + hltDiMuonL3PreFiltered8 + hltSingleMu13L3Filtered13 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_TripleMu5_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu3 + hltPreTripleMu3 + hltL1DoubleMu3L1TriMuFiltered3 + HLTL2muonrecoSequence + hltL1DoubleMu3L2TriMuFiltered3 + HLTL3muonrecoSequence + hltL1DoubleMu3L3TriMuFiltered5 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu5_L2Mu2_Jpsi_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleMu0 + hltPreMu5L2Mu2Jpsi + hltMu5L2Mu2L1Filtered0 + HLTL2muonrecoSequence + hltMu5L2Mu2L2PreFiltered0 + HLTL3muonrecoSequence + hltMu5L2Mu2L3Filtered5 + hltMu5L2Mu2JpsiTrackMassFiltered + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon20_CaloIdVL_IsoL_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton20CaloIdVLIsoL + HLTPhoton20CaloIdVLIsoLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon20_R9Id_Photon18_R9Id_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton20R9IdPhoton18R9Id + HLTPhoton20R9IdPhoton18R9IdSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVL + HLTPhoton20CaloIdVTIsoTSequence + HLTEle8CaloIdLCaloIsoVLNoL1SeedSequence + hltPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLDoubleLegCombFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_Photon18_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26Photon18 + HLTPhoton26Photon18Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_IsoVL_Photon18_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26IsoVLPhoton18 + HLTPhoton26IsoVLPhoton18Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_IsoVL_Photon18_IsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26IsoVLPhoton18IsoVL + HLTPhoton26IsoVLPhoton18IsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_CaloIdL_IsoVL_Photon18_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26CaloIdLIsoVLPhoton18 + HLTPhoton26CaloIdLIsoVLPhoton18Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26CaloIdLIsoVLPhoton18R9Id + HLTPhoton26CaloIdLIsoVLPhoton18R9IdSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVL + HLTPhoton26CaloIdLIsoVLPhoton18CaloIdLIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26R9IdPhoton18CaloIdLIsoVL + HLTPhoton26R9IdPhoton18CaloIdLIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon26_R9Id_Photon18_R9Id_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPrePhoton26R9IdPhoton18R9Id + HLTPhoton26R9IdPhoton18R9IdSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon30_CaloIdVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG15 + hltPrePhoton30CaloIdVL + HLTPhoton30CaloIdVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon30_CaloIdVL_IsoL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG15 + hltPrePhoton30CaloIdVLIsoL + HLTPhoton30CaloIdVLIsoLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton32CaloIdLPhoton26CaloIdL + HLTPhoton32CaloIdLPhoton26CaloIdLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon36_IsoVL_Photon22_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton36IsoVLPhoton22 + HLTPhoton36IsoVLPhoton22Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton36CaloIdLPhoton22CaloIdL + HLTPhoton36CaloIdLPhoton22CaloIdLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon36_CaloIdL_IsoVL_Photon22_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton36CaloIdLIsoVLPhoton22 + HLTPhoton36CaloIdLIsoVLPhoton22Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40CaloIdLPhoton28CaloIdL + HLTPhoton40CaloIdLPhoton28CaloIdLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon50_CaloIdVL_IsoL_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton50CaloIdVLIsoL + HLTPhoton50CaloIdVLIsoLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon70_CaloIdL_HT300_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton70CaloIdLHT300 + HLTSinglePhoton70CaloIdLSequence + HLTRecoJetSequenceAK5Corrected + hltHT300 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon70_CaloIdL_HT350_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton70CaloIdLHT350 + HLTSinglePhoton70CaloIdLSequence + HLTRecoJetSequenceAK5Corrected + hltHT350 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon70_CaloIdL_MHT50_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton70CaloIdLMHT50 + HLTSinglePhoton70CaloIdLSequence + HLTRecoJetSequenceAK5Corrected + hltMHT50 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon70_CaloIdL_MHT70_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton70CaloIdLMHT70 + HLTSinglePhoton70CaloIdLSequence + HLTRecoJetSequenceAK5Corrected + hltMHT70 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon75_CaloIdVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton75CaloIdVL + HLTPhoton75CaloIdVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon75_CaloIdVL_IsoL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton75CaloIdVLIsoL + HLTPhoton75CaloIdVLIsoLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon90_CaloIdVL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton90CaloIdVL + HLTPhoton90CaloIdVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon90_CaloIdVL_IsoL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton90CaloIdVLIsoL + HLTPhoton90CaloIdVLIsoLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon125_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton125 + HLTSinglePhoton125Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon200_NoHE_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton200NoHE + HLTSinglePhoton200NoHESequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton33_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton33 + HLTDoublePhoton33Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton33_HEVT_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton33HEVT + HLTDoublePhoton33HEVTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton50_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton50 + HLTPhoton50Sequence + HLTDoublePhoton50UnseededLegSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton60_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton60 + HLTPhoton60Sequence + HLTDoublePhoton60UnseededLegSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton5_IsoVL_CEP_v3 = cms.Path( HLTBeginSequence + hltL1sL1DoubleEG2FwdVeto + hltPreDoublePhoton5IsoVLCEP + HLTDoublePhoton5IsoVLSequence + hltTowerMakerForHcal + hltHcalTowerFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleEG5_v2 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreL1SingleEG5 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1SingleEG12_v2 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreL1SingleEG12 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8 + HLTEle8Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdL_CaloIsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdLCaloIsoVL + HLTEle8CaloIdLCaloIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdL_TrkIdVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdLTrkIdVL + HLTEle8CaloIdLTrkIdVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoT + HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdL_CaloIsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdLCaloIsoVL + HLTEle17CaloIdLCaloIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVL + HLTEle17CaloIdIsoEle8CaloIdIsoSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30 + HLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass30 + HLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass30Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdLCaloIsoVLEle15HFL + HLTSingleElectronEt17CaloIdIsoSequence + HLTHFEM15Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele27_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG15 + hltPreEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVL + HLTEle27CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG15 + hltPreEle27CaloIdVTCaloIsoTTrkIdTTrkIsoT + HLTEle27CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVL + HLTEle32CaloIdVLCaloIsoVLTrkIdVLTrkIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreEle32CaloIdVTCaloIsoTTrkIdTTrkIsoT + HLTEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleEG20 + hltPreEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17 + HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Sequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele45_CaloIdVT_TrkIdT_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreEle45CaloIdVTTrkIdT + HLTEle45CaloIdVTTrkIdTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele52_CaloIdVT_TrkIdT_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreEle52CaloIdVTTrkIdT + HLTEle52CaloIdVTTrkIdTSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle8_CaloIdL_TrkIdVL_v1 = cms.Path( HLTBeginSequence + hltL1sL1DoubleEG5 + hltPreDoubleEle8CaloIdLTrkIdVL + HLTDoubleEle8L1NonIsoHLTCaloIdLSequence + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8OneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8DetaFilter + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8DphiFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle33_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoubleEle33 + HLTPhoton33Sequence + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle33PixelMatchFilter + HLTDoublePhoton33UnseededLegSequence + HLTActivityPixelMatchSequence + hltDiEle33PixelMatchDoubleFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle33_CaloIdL_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoubleEle33CaloIdL + HLTPhoton33Sequence + hltL1IsoHLTClusterShape + hltL1NonIsoHLTClusterShape + hltEG33CaloIdLClusterShapeFilter + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL1IsoStartUpElectronPixelSeeds + hltL1NonIsoStartUpElectronPixelSeeds + hltEle33CaloIdLPixelMatchFilter + HLTDoublePhoton33UnseededLegSequence + hltActivityPhotonClusterShape + hltDoubleEG33CaloIdLClusterShapeDoubleFilter + HLTActivityPixelMatchSequence + hltDiEle33CaloIdLPixelMatchDoubleFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoPFTau35_Trk20_MET45_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sSingleIsoTau35Trk20MET45 + hltPreSingleIsoTau35Trk20MET45 + HLTL2TauJetsSequence + hltFilterL2EtCutSingleIsoPFTau35Trk20MET45 + HLTRecoMETSequence + hltMet45 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequenceForTaus + HLTPFJetsSequenceForTaus + HLTPFTauTightIsoSequence + hltPFTauTightIso35 + hltPFTauTightIso35Track + hltPFTauTightIsoTrackPt20Discriminator + hltSelectedPFTauTightIsoTrackPt20 + hltConvPFTauTightIsoTrackPt20 + hltFilterSingleIsoPFTau35Trk20LeadTrackPt20 + hltSelectedPFTauTightIsoTrackPt20Isolation + hltConvPFTauTightIsoTrackPt20Isolation + hltL1HLTSingleIsoPFTau35Trk20Met45JetsMatch + hltFilterSingleIsoPFTau35Trk20MET45LeadTrack20MET45IsolationL1HLTMatched + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleIsoPFTau25_Trk5_eta2p1_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sDoubleIsoTau25Trk5 + hltPreDoubleIsoTau25Trk5 + HLTL2TauJetsSequence + hltFilterL2EtCutDoublePFIsoTau25Trk5 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequenceForTaus + HLTPFJetsSequenceForTaus + HLTPFTauTightIsoSequence + hltDoublePFTauTightIso25Track + hltDoublePFTauTightIso25Track5 + hltL1HLTDoubleIsoPFTau25Trk5JetsMatch + hltFilterDoubleIsoPFTau25Trk5LeadTrack5IsolationL1HLTMatched + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BTagMu_DiJet20_Mu5_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet16Central + hltPreBTagMuDiJet20Mu5 + HLTRecoJetSequenceAK5Corrected + hltBDiJet20Central + HLTBTagMuDiJet20SequenceL25 + hltBSoftMuonDiJet20L25FilterByDR + HLTBTagMuDiJet20Mu5SelSequenceL3 + hltBSoftMuonDiJet20Mu5L3FilterByDR + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BTagMu_DiJet40_Mu5_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet20Central + hltPreBTagMuDiJet40Mu5 + HLTRecoJetSequenceAK5Corrected + hltBDiJet40Central + HLTBTagMuDiJet40SequenceL25 + hltBSoftMuonDiJet40L25FilterByDR + HLTBTagMuDiJet40Mu5SelSequenceL3 + hltBSoftMuonDiJet40Mu5L3FilterByDR + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BTagMu_DiJet70_Mu5_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet28Central + hltPreBTagMuDiJet70Mu5 + HLTRecoJetSequenceAK5Corrected + hltBDiJet70Central + HLTBTagMuDiJet70SequenceL25 + hltBSoftMuonDiJet70L25FilterByDR + HLTBTagMuDiJet70Mu5SelSequenceL3 + hltBSoftMuonDiJet70Mu5L3FilterByDR + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BTagMu_DiJet110_Mu5_v4 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet28Central + hltPreBTagMuDiJet110Mu5 + HLTRecoJetSequenceAK5Corrected + hltBDiJet110Central + HLTBTagMuDiJet110SequenceL25 + hltBSoftMuonDiJet110L25FilterByDR + HLTBTagMuDiJet110Mu5SelSequenceL3 + hltBSoftMuonDiJet110Mu5L3FilterByDR + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v5 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu5Ele8CaloIdLTrkIdVLEle8 + hltL1MuOpenEG5L1Filtered3 + HLTL2muonrecoSequence + hltL2MuOpenEG5L2Filtered3 + HLTL3muonrecoSequence + hltMuOpenEG5L3Filtered5 + HLTMu5Ele8CaloIdLTrkIdVLEle8L1NonIsoHLTnonIsoSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_R005_MR200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet20Central + hltL1sL1DoubleJet36Central + hltPreMu8R005MR200 + hltL1Mu3Jet20L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu8Jet20L2Filtered3 + HLTL3muonrecoSequence + hltL3Mu8Jet20L3Filtered8 + HLTRSequenceNoJetFilter + hltR005MR200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_R020_MR200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet20Central + hltL1sL1DoubleJet36Central + hltPreMu8R020MR200 + hltL1Mu3Jet20L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu8Jet20L2Filtered3 + HLTL3muonrecoSequence + hltL3Mu8Jet20L3Filtered8 + HLTRSequenceNoJetFilter + hltR020MR200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_R025_MR200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu3Jet20Central + hltL1sL1DoubleJet36Central + hltPreMu8R025MR200 + hltL1Mu3Jet20L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu8Jet20L2Filtered3 + HLTL3muonrecoSequence + hltL3Mu8Jet20L3Filtered8 + HLTRSequenceNoJetFilter + hltR025MR200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_Mu5_PFMHT35_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT250Mu5pfMHT35 + cms.ignore(hltL1sL1SingleMuOpenCandidate) + HLTRecoJetSequenceAK5Corrected + hltHT250 + hltL1HTT100L1MuFiltered3 + HLTL2muonrecoSequence + hltL1HTT100singleMuL2PreFiltered3 + HLTL3muonrecoSequence + hltL1HTT100singleMuL3PreFiltered5 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT35Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT250_Mu15_PFMHT20_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT250Mu15pfMHT20 + cms.ignore(hltL1sL1SingleMuOpenCandidate) + HLTRecoJetSequenceAK5Corrected + hltHT250 + hltHTT100L1MuFiltered0 + HLTL2muonrecoSequence + hltL1HTT100singleMuL2PreFiltered10 + HLTL3muonrecoSequence + hltL1HTT100singleMuL3PreFiltered15 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT20Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT300_Mu5_PFMHT40_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT300Mu5pfMHT40 + cms.ignore(hltL1sL1SingleMuOpenCandidate) + HLTRecoJetSequenceAK5Corrected + hltHT300 + hltL1HTT100L1MuFiltered3 + HLTL2muonrecoSequence + hltL1HTT100singleMuL2PreFiltered3 + HLTL3muonrecoSequence + hltL1HTT100singleMuL3PreFiltered5 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT40Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_HT350_Mu5_PFMHT45_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1HTT100 + hltPreHT350Mu5pfMHT45 + cms.ignore(hltL1sL1SingleMuOpenCandidate) + HLTRecoJetSequenceAK5Corrected + hltHT350 + hltL1HTT100L1MuFiltered3 + HLTL2muonrecoSequence + hltL1HTT100singleMuL2PreFiltered3 + HLTL3muonrecoSequence + hltL1HTT100singleMuL3PreFiltered5 + HLTRecoJetSequencePrePF + HLTTrackReconstructionForJets + HLTParticleFlowSequence + HLTPFJetsSequence + hltPFMHT45Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu3_DiJet30_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreMu3DiJet30 + hltSingleMuOpenL1Filtered + HLTL2muonrecoSequence + hltSingleMu3L2Filtered0 + HLTL3muonrecoSequence + hltSingleMu3L3Filtered3 + HLTRecoJetSequenceAK5Corrected + hltDoubleJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu3_TriJet30_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreMu3TriJet30 + hltSingleMuOpenL1Filtered + HLTL2muonrecoSequence + hltSingleMu3L2Filtered0 + HLTL3muonrecoSequence + hltSingleMu3L3Filtered3 + HLTRecoJetSequenceAK5Corrected + hltTripleJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu3_QuadJet30_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMuOpen + hltPreMu3QuadJet30 + hltSingleMuOpenL1Filtered + HLTL2muonrecoSequence + hltSingleMu3L2Filtered0 + HLTL3muonrecoSequence + hltSingleMu3L3Filtered3 + HLTRecoJetSequenceAK5Corrected + hltQuadJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_Ele17_CaloIdL_v4 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu8Ele17CaloIdL + hltL1MuOpenEG5L1Filtered5 + HLTL2muonrecoSequence + hltL1MuOpenEG5L2Filtered5 + HLTL3muonrecoSequence + hltL1MuOpenEG5L3Filtered8 + HLTDoEGammaStartupSequence + hltEGRegionalL1MuOpenEG5 + hltEG17EtFilterL1MuOpenEG5 + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLTCaloIdLMu8Ele17ClusterShapeFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTNonIsoMu8Ele17HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_Photon20_CaloIdVT_IsoT_v4 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu8Photon20CaloIdVTIsoT + HLTPhoton20CaloIdVTIsoTMu8Sequence + hltL1SingleMuOpenEG5L1Filtered0 + HLTL2muonrecoSequence + hltSingleMu5EG5L2Filtered3 + HLTL3muonrecoSequence + hltSingleMu8EG5L3Filtered8 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu8_Jet40_v5 = cms.Path( HLTBeginSequence + hltL1sL1Mu3Jet20Central + hltPreMu8Jet40 + hltL1Mu3Jet20L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu8Jet20L2Filtered3 + HLTL3muonrecoSequence + hltL3Mu8Jet20L3Filtered8 + HLTRecoJetSequenceAK5Corrected + hltJet40 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu15_Photon20_CaloIdL_v5 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu15Photon20CaloIdL + hltL1MuOpenEG5L1Filtered5 + HLTL2muonrecoSequence + hltL1MuOpenEG5L2Filtered5 + HLTL3muonrecoSequence + hltL1MuOpenEG5L3Filtered15 + HLTDoEGammaStartupSequence + hltEGRegionalL1MuOpenEG5 + hltEG20EtFilterMuOpenEG5 + HLTDoEgammaClusterShapeSequence + hltMu15Photon20CaloIdLClusterShapeFilter + HLTDoEGammaHESequence + hltMu15Photon20CaloIdLHEFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu15_DoublePhoton15_CaloIdL_v5 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu15DoublePhoton15CaloIdL + hltL1MuOpenEG5L1Filtered5 + HLTL2muonrecoSequence + hltL1MuOpenEG5L2Filtered5 + HLTL3muonrecoSequence + hltL1MuOpenEG5L3Filtered15 + HLTDoEGammaStartupSequence + hltEGRegionalL1MuOpenEG5 + hltDoubleEG15EtFilterL1MuOpenEG5 + HLTDoEgammaClusterShapeSequence + hltMu15DiPhoton15CaloIdLClusterShapeFilter + HLTDoEGammaHESequence + hltMu15DiPhoton15CaloIdLHEFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu15_LooseIsoPFTau15_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreMu15IsoPFTau15 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + HLTL3muonrecoSequence + hltSingleMu15L3Filtered15 + HLTRecoJetSequencePrePF + hltTauJet5 + HLTPFJetTriggerSequenceForTaus + hltPFJet15 + HLTPFTauSequence + hltPFTau15 + hltPFTau15Track + hltPFTau15TrackLooseIso + hltOverlapFilterMu15IsoPFTau15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu17_CentralJet30_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu17CentralJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTL3muonrecoSequence + hltL1Mu7CenJetL3MuFiltered17 + HLTRecoJetSequenceAK5Corrected + hltJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu17_DiCentralJet30_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu17DiCentralJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTL3muonrecoSequence + hltL1Mu7CenJetL3MuFiltered17 + HLTRecoJetSequenceAK5Corrected + hltDiJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu17_TriCentralJet30_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu17TriCentralJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTL3muonrecoSequence + hltL1Mu7CenJetL3MuFiltered17 + HLTRecoJetSequenceAK5Corrected + hltTriJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu17_Ele8_CaloIdL_v4 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreMu17Ele8CaloIdL + hltL1MuOpenEG5L1Filtered12 + HLTL2muonrecoSequence + hltL1MuOpenEG5L2Filtered12 + HLTL3muonrecoSequence + hltL1MuOpenEG5L3Filtered17 + HLTDoEGammaStartupSequence + hltEGRegionalL1MuOpenEG5 + hltEG8EtFilterMuOpenEG5 + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLTCaloIdLMu17Ele8ClusterShapeFilter + HLTDoEGammaHESequence + hltL1NonIsoHLTNonIsoMu17Ele8HEFilter + HLTDoEGammaPixelSequence + hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu12_CentralJet30_BTagIP_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu12BTagIPCenJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTRecoJetSequenceAK5Corrected + hltBJet30Central + HLTBTagIPSequenceL25SingleTop + hltBLifetimeL25FilterSingleTop + HLTL3muonrecoSequence + hltL1Mu7CenJetL3MuFiltered12 + HLTBTagIPSequenceL3SingleTop + hltBLifetimeL3FilterSingleTop + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu17_CentralJet30_BTagIP_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreMu17BTagIPCenJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTRecoJetSequenceAK5Corrected + hltBJet30Central + HLTBTagIPSequenceL25SingleTop + hltBLifetimeL25FilterSingleTop + HLTL3muonrecoSequence + hltL1Mu7CenJetL3MuFiltered17 + HLTBTagIPSequenceL3SingleTop + hltBLifetimeL3FilterSingleTop + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Mu15_HT200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu0HTT50 + hltPreMu15HT200 + hltL1Mu0HTT50L1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu0HTT50L2MuFiltered10 + HLTL3muonrecoSequence + hltL1Mu0HTT50L3MuFiltered15 )
HLT_Mu20_HT200_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu0HTT50 + hltPreMu20HT200 + hltL1Mu0HTT50L1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu0HTT50L2MuFiltered12 + HLTL3muonrecoSequence + hltL1Mu0HTT50L3MuFiltered20 )
HLT_IsoMu15_LooseIsoPFTau15_v1 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu10 + hltPreIsoMu15IsoPFTau15 + hltL1SingleMu10L1Filtered0 + HLTL2muonrecoSequence + hltL2Mu10L2Filtered10 + HLTL2muonisorecoSequence + hltSingleMuIsoL2IsoFiltered10 + HLTL3muonrecoSequence + hltSingleMuIsoL3PreFiltered15 + HLTL3muonisorecoSequence + hltSingleMuIsoL3IsoFiltered15 + HLTRecoJetSequencePrePF + hltTauJet5 + HLTPFJetTriggerSequenceForTaus + HLTPFTauSequence + hltPFTau15 + hltPFTau15Track + hltPFTau15TrackLooseIso + hltOverlapFilterIsoMu15IsoPFTau15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_IsoMu17_CentralJet30_BTagIP_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1SingleMu7 + hltPreIsoMu17BTagIPCentJet30 + hltL1Mu7CenJetL1MuFiltered0 + HLTL2muonrecoSequence + hltL1Mu7CenJetL2MuFiltered7 + HLTL2muonisorecoSequence + hltIsoMu7CenJet30L2IsoFiltered7 + HLTRecoJetSequenceAK5Corrected + hltBJet30Central + HLTBTagIPSequenceL25SingleTop + hltBLifetimeL25FilterSingleTop + HLTL3muonrecoSequence + hltIsoMu17CenJet30L3Filtered17 + HLTL3muonisorecoSequence + hltIsoMu17CenJet30L3IsoFiltered17 + HLTBTagIPSequenceL3SingleTop + hltBLifetimeL3FilterSingleTop + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu3_HT150_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu0HTT50 + hltPreDoubleMu3HT150 + hltL1Mu0HTT50L1DiMuFiltered0 + HLTL2muonrecoSequence + hltL1Mu0HTT50L2DiMuFiltered0 + HLTL3muonrecoSequence + hltL1Mu0HTT50L3DiMuFiltered3 + HLTRecoJetSequenceAK5Corrected + hltHT150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu3_HT200_v5 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1Mu0HTT50 + hltPreDoubleMu3HT200 + hltL1Mu0HTT50L1DiMuFiltered0 + HLTL2muonrecoSequence + hltL1Mu0HTT50L2DiMuFiltered0 + HLTL3muonrecoSequence + hltL1Mu0HTT50L3DiMuFiltered3 + HLTRecoJetSequenceAK5Corrected + hltHT200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleMu5_Ele8_v5 = cms.Path( HLTBeginSequence + hltL1sL1MuOpenEG5 + hltPreDoubleMu5Ele8 + hltL1MuOpenEG5L1DiMuFiltered3 + HLTL2muonrecoSequence + hltL1MuOpenEG5L2DiMuFiltered3 + HLTL3muonrecoSequence + hltL1MuOpenEG5L3DiMuFiltered5 + HLTDoubleMu5Ele8L1NonIsoHLTnonIsoSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_R005_MR150_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40R014MR150 + HLTSinglePhoton40CaloIdLSequence + HLTRSequenceNoJetFilter + hltR005MR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_R014_MR450_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40R014MR450 + HLTSinglePhoton40CaloIdLSequence + HLTRSequenceNoJetFilter + hltR014MR450 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_R020_MR300_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40R020MR300 + HLTSinglePhoton40CaloIdLSequence + HLTRSequenceNoJetFilter + hltR020MR300 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_R025_MR200_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40R025MR200 + HLTSinglePhoton40CaloIdLSequence + HLTRSequenceNoJetFilter + hltR025MR200 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Photon40_R038_MR150_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPrePhoton40R038MR150 + HLTSinglePhoton40CaloIdLSequence + HLTRSequenceNoJetFilter + hltR038MR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton40_MR150_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton40MR150 + HLTDoublePhoton40Sequence + HLTRSequenceNoJetFilter + hltMR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoublePhoton40_R014_MR150_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG20 + hltPreDoublePhoton40R014MR150 + HLTDoublePhoton40Sequence + HLTRSequenceNoJetFilter + hltR014MR150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdT_TrkIdT_DiJet30_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdTTrkIdTDiJet30 + HLTEle8CaloIdTTrkIdTSequence + HLTRecoJetSequenceAK5Corrected + hltDoubleJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdT_TrkIdT_TriJet30_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdTTrkIdTDiJet30 + HLTEle8CaloIdTTrkIdTSequence + HLTRecoJetSequenceAK5Corrected + hltTripleJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdT_TrkIdT_QuadJet30_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdTTrkIdTDiJet30 + HLTEle8CaloIdTTrkIdTSequence + HLTRecoJetSequenceAK5Corrected + hltQuadJet30Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG5 + hltPreEle8CaloIdLCaloIsoVLJet40 + HLTEle8CaloIdLCaloIsoVLSequence + HLTRecoJetSequenceAK5Corrected + hltAntiKT5L2L3CaloJetsEle8CaloIdLCaloIsoVLRemoved + hltJet40Ele8CaloIdLCaloIsoVLRemoved + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v4 = cms.Path( HLTBeginSequence + hltL1sL1EG5HTT75 + hltPreEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLHT200 + HLTRecoJetSequenceAK5Corrected + hltHT200 + HLTEle15L1EG5HTT75CaloIdTCaloIsoVLTrkIdTTrkIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT250_v4 = cms.Path( HLTBeginSequence + hltL1sL1EG5HTT75 + hltPreEle15CaloIdTCaloIsoVLTrkIdTTrkIsoVLHT250 + HLTRecoJetSequenceAK5Corrected + hltHT250 + HLTEle15L1EG5HTT75CaloIdTCaloIsoVLTrkIdTTrkIsoVLSequence + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v5 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTTrkIdTLooseIsoPFTau15 + HLTEle15CaloIdVTTrkIdTSequence + HLTRecoJetSequencePrePF + hltTauJet5 + hltOverlapFilterEle15CaloJet5 + HLTPFJetTriggerSequenceForTaus + hltPFJet15 + HLTPFTauSequence + hltPFTau15 + hltPFTau15Track + hltPFTau15TrackLooseIso + hltOverlapFilterEle15IsoPFTau15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v5 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTLooseIsoPFTau15 + HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequencePrePF + hltTauJet5 + hltOverlapFilterIsoEle15CaloJet5 + HLTPFJetTriggerSequenceForTaus + hltPFJet15 + HLTPFTauSequence + hltPFTau15 + hltPFTau15Track + hltPFTau15TrackLooseIso + hltOverlapFilterIsoEle15IsoPFTau15 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v5 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTLooseIsoPFTau20 + HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequencePrePF + hltTauJet5 + hltOverlapFilterIsoEle15CaloJet5 + HLTPFJetTriggerSequenceForTaus + hltPFJet20 + HLTPFTauSequence + hltPFTau20 + hltPFTau20Track + hltPFTau20TrackLooseIso + hltOverlapFilterIsoEle15IsoPFTau20 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta3_v2 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta3 + HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets + hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta3 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta2_v2 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta2 + HLTEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets + hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTJet35Jet25Deta2 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele15_CaloIdVT_TrkIdT_Jet35_Jet25_Deta2_v2 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle15CaloIdVTTrkIdTJet35Jet25Deta2 + HLTEle15CaloIdVTTrkIdTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle15CaloIdVTTrkIdTFromAK5CorrJets + hltEle15CaloIdVTTrkIdTJet35Jet25Deta2 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v4 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle25CaloIdVTTrkIdTCentralJet30 + HLTEle25CaloIdVTCaloTrkIdSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets + hltEle25CaloIdVTTrkIdTCentralJet30Cleaned + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele25_CaloIdVT_TrkIdT_DiCentralJet30_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle25CaloIdVTTrkIdTCentralDiJet30 + HLTEle25CaloIdVTCaloTrkIdSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets + hltEle25CaloIdVTTrkIdTCentralDiJet30Cleaned + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle25CaloIdVTTrkIdTCentralTriJet30 + HLTEle25CaloIdVTCaloTrkIdSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets + hltEle25CaloIdVTTrkIdTCentralTriJet30Cleaned + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v3 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle25CaloIdVTTrkIdTCentralJet30BTagIP + HLTEle25CaloIdVTCaloTrkIdSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle25CaloIdVTTrkIdTFromAK5CorrJets + hltSingleEleCleanBJet30Central + HLTBTagIPSequenceL25EleJetSingleTop + hltBLifetimeL25FilterEleJetSingleTop + HLTBTagIPSequenceL3EleJetSingleTop + hltBLifetimeL3FilterEleJetSingleTop + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdVT_TrkIdT_CentralJet30_CentralJet25_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdVTTrkIdTCenJet30CenJet25 + HLTEle17CaloIdVTTrkIdTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle17CaloIdVTTrkIdTFromAK5CorrJets + hltEle17CaloIdVTTrkIdTCentralJet30Cleaned  + hltEle17CaloIdVTTrkIdTCentralDiJet25Cleaned + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele17_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT15_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCenJet30CenJet25PFMHT15 + HLTEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets + hltEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralJet30Cleaned + hltEle17CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralDiJet25Cleaned + HLTRecoJetSequencePrePF + HLTPFJetTriggerSequence + hltPFMHT15Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v1 = cms.Path( HLTBeginSequence + hltL1sL1SingleEG12 + hltPreEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCenJet30CenJet25PFMHT20 + HLTEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTSequence + HLTRecoJetSequenceAK5Corrected + hltCleanEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTFromAK5CorrJets + hltEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralJet30Cleaned + hltEle25CaloIdVTCaloIsoTTrkIdTTrkIsoTCentralDiJet25Cleaned + HLTRecoJetSequencePrePF + HLTPFJetTriggerSequence + hltPFMHT20Filter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v2 = cms.Path( HLTBeginSequence + hltL1sL1DoubleEG5HTT50 + hltPreDoubleEle8CaloIdLTrkIdVLHT150 + HLTDoubleEle8HTT50L1NonIsoHLTCaloIdLSequence + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50DetaFilter + hltL1NonIsoHLTCaloIdLTrkIdVLDoubleEle8HTT50DphiFilter + HLTRecoJetSequenceAK5Corrected + hltHT150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v2 = cms.Path( HLTBeginSequence + hltL1sL1DoubleEG5HTT50 + hltPreDoubleEle8CaloIdTTrkIdVLHT160 + HLTDoubleEle8HTT50L1NonIsoHLTCaloIdTSequence + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50OneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50DetaFilter + hltL1NonIsoHLTCaloIdTTrkIdVLDoubleEle8HTT50DphiFilter + HLTRecoJetSequenceAK5Corrected + hltHT150 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v4 = cms.Path( HLTBeginSequence + hltL1sL1TripleEG5 + hltPreDoubleEle10CaloIdLTrkIdVLEle10 + HLTTripleElectronEt10L1NonIsoHLTNonIsoSequence + hltL1NonIsoHLT2CaloIdLTripleElectronEt10HEFilter + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLT2LegEleIdTripleElectronEt10ClusterShapeFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoHLT2LegEleIdTripleElectronEt10OneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltL1NonIsoHLT2LegEleIdTripleElectronEt10EleIdDetaFilter + hltL1NonIsoHLT2LegEleIdTripleElectronEt10EleIdDphiFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_TripleEle10_CaloIdL_TrkIdVL_v4 = cms.Path( HLTBeginSequence + hltL1sL1TripleEG5 + hltPreTripleEle10CaloIdLTrkIdVL + HLTTripleElectronEt10L1NonIsoHLTNonIsoSequence + HLTDoEgammaClusterShapeSequence + hltL1NonIsoHLT3LegEleIdTripleElectronEt10ClusterShapeFilter + HLTPixelMatchElectronL1IsoTrackingSequence + HLTPixelMatchElectronL1NonIsoTrackingSequence + hltL1NonIsoHLT3LegEleIdTripleElectronEt10OneOEMinusOneOPFilter + HLTDoElectronDetaDphiSequence + hltL1NonIsoHLT3LegEleIdTripleElectronEt10EleIdDetaFilter + hltL1NonIsoHLT3LegEleIdTripleElectronEt10EleIdDphiFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_PixelTracks_Multiplicity80_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sETT220 + hltPrePixelTracksMultiplicity80 + HLTDoLocalPixelSequence + hltPixelClusterShapeFilter + HLTRecopixelvertexingForHighMultSequence + hltPixelCandsForHighMult + hlt1HighMult80 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_PixelTracks_Multiplicity100_v3 = cms.Path( HLTBeginSequenceBPTX + hltL1sETT220 + hltPrePixelTracksMultiplicity100 + HLTDoLocalPixelSequence + hltPixelClusterShapeFilter + HLTRecopixelvertexingForHighMultSequence + hltPixelCandsForHighMult + hlt1HighMult100 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BeamGas_HF_v4 = cms.Path( HLTBeginSequence + hltL1sL1BeamGasHf + hltPreL1BeamGasHf + hltHcalDigis + hltHfreco + hltHFAsymmetryFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BeamGas_BSC_v3 = cms.Path( HLTBeginSequence + hltL1sL1BeamGasBsc + hltPreL1BeamGasBsc + HLTDoLocalPixelLight + hltPixelActivityFilter + hltPixelAsymmetryFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_BeamHalo_v3 = cms.Path( HLTBeginSequence + hltL1sL1BeamHalo + hltPreL1BeamHalo + HLTDoLocalPixelLight + hltPixelActivityFilterForHalo + HLTDoLocalStripSequence + hltTrackerHaloFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1Tech_BSC_minBias_threshold1_v4 = cms.Path( HLTBeginSequence + hltL1sZeroBias + hltPreL1TechBSCminBiasthreshold1 + hltL1TechBSCminBiasthreshold1 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1Tech_BSC_halo_v4 = cms.Path( HLTBeginSequence + hltL1sZeroBias + hltPreL1TechBSChalo + hltL1TechBSChalo + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1_PreCollisions_v2 = cms.Path( HLTBeginSequence + hltL1sL1PreCollisions + hltPreL1PreCollisions + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1_Interbunch_BSC_v2 = cms.Path( HLTBeginSequence + hltL1sL1InterbunchBsc + hltPreL1Interbunch1 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_GlobalRunHPDNoise_v3 = cms.Path( HLTBeginSequence + hltL1sGlobalRunHPDNoise + hltPreGlobalRunHPDNoise + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1Tech_HBHEHO_totalOR_v2 = cms.Path( HLTBeginSequence + hltL1sTechTrigHCALNoise + hltPreTechTrigHCALNoise + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_ZeroBias_v3 = cms.Path( HLTBeginSequence + hltL1sZeroBias + hltPreZeroBias + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Physics_v1 = cms.Path( HLTBeginSequence + hltPrePhysics + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_Physics_NanoDST_v1 = cms.Path( HLTBeginSequence + hltPrePhysicsNanoDST + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1TrackerCosmics_v3 = cms.Path( HLTBeginSequence + hltL1sTrackerCosmics + hltPreTrackerCosmics + hltTrackerCosmicsPattern + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_RegionalCosmicTracking_v3 = cms.Path( HLTBeginSequence + hltL1sTrackerCosmics + hltPreRegionalCosmicTracking + hltTrackerCosmicsPattern + hltL1sL1SingleMuOpenCandidate + hltL1MuORL1Filtered0 + HLTL2muonrecoSequenceNoVtx + hltSingleL2MuORL2PreFilteredNoVtx + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltRegionalCosmicTrackerSeeds + hltRegionalCosmicCkfTrackCandidates + hltRegionalCosmicTracks + hltCosmicTrackSelector + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L3MuonsCosmicTracking_v3 = cms.Path( HLTBeginSequence + hltL1sTrackerCosmics + hltPreL3MuonsCosmicTracking + hltTrackerCosmicsPattern + hltL1sL1SingleMuOpenCandidate + hltL1MuORL1Filtered0 + HLTL2muonrecoSequenceNoVtx + hltSingleL2MuORL2PreFilteredNoVtx + HLTDoLocalPixelSequence + HLTDoLocalStripSequence + hltL3TrajectorySeedNoVtx + hltL3TrackCandidateFromL2NoVtx + hltL3TkTracksFromL2NoVtx + hltL3MuonsNoVtx + hltL3MuonCandidatesNoVtx + hltMu5NoVertexL3PreFiltered5 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_LogMonitor_v1 = cms.Path( HLTBeginSequence + hltPreLogMonitor + hltLogMonitorFilter + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1ETM30_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1ETM30 + hltPreL1ETM30 + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1DoubleJet36Central_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1DoubleJet36Central + hltPreL1DoubleJet36Central + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLT_L1MultiJet_v2 = cms.Path( HLTBeginSequenceBPTX + hltL1sL1MultiJet + hltPreL1MultiJet + cms.SequencePlaceholder( "HLTEndSequence" ) )
HLTriggerFinalPath = cms.Path( HLTBeginSequence + hltFEDSelector + hltTriggerSummaryAOD + hltTriggerSummaryRAW + hltBoolTrue )


HLTSchedule = cms.Schedule( *(HLTriggerFirstPath, HLT_Activity_Ecal_SC7_v3, HLT_L1SingleJet16_v2, HLT_L1SingleJet36_v2, HLT_Jet30_v3, HLT_Jet60_v3, HLT_Jet80_v3, HLT_Jet110_v3, HLT_Jet150_v3, HLT_Jet190_v3, HLT_Jet240_v3, HLT_Jet300_v2, HLT_Jet370_v3, HLT_Jet370_NoJetID_v3, HLT_DiJetAve30_v3, HLT_DiJetAve60_v3, HLT_DiJetAve80_v3, HLT_DiJetAve110_v3, HLT_DiJetAve150_v3, HLT_DiJetAve190_v3, HLT_DiJetAve240_v3, HLT_DiJetAve300_v3, HLT_DiJetAve370_v3, HLT_DoubleJet30_ForwardBackward_v4, HLT_DoubleJet60_ForwardBackward_v4, HLT_DoubleJet70_ForwardBackward_v4, HLT_DoubleJet80_ForwardBackward_v4, HLT_DiJet130_PT130_v2, HLT_DiJet160_PT160_v2, HLT_CentralJet80_MET65_v3, HLT_CentralJet80_MET80HF_v2, HLT_CentralJet80_MET100_v3, HLT_CentralJet80_MET160_v3, HLT_DiJet60_MET45_v3, HLT_DiCentralJet20_MET80_v1, HLT_DiCentralJet20_BTagIP_MET65_v1, HLT_QuadJet40_v4, HLT_QuadJet40_IsoPFTau40_v4, HLT_QuadJet50_BTagIP_v3, HLT_QuadJet50_Jet40_v3, HLT_QuadJet60_v3, HLT_QuadJet70_v3, HLT_ExclDiJet60_HFOR_v3, HLT_ExclDiJet60_HFAND_v3, HLT_HT150_v4, HLT_HT150_AlphaT0p60_v3, HLT_HT200_v4, HLT_HT200_AlphaT0p53_v2, HLT_HT200_AlphaT0p60_v3, HLT_HT250_v4, HLT_HT250_AlphaT0p53_v2, HLT_HT250_AlphaT0p54_v2, HLT_HT250_MHT60_v4, HLT_HT250_MHT70_v1, HLT_HT300_v5, HLT_HT300_MHT75_v5, HLT_HT300_PFMHT55_v1, HLT_HT300_CentralJet30_BTagIP_v1, HLT_HT300_CentralJet30_BTagIP_PFMHT55_v1, HLT_HT300_CentralJet30_BTagIP_PFMHT75_v1, HLT_HT300_AlphaT0p52_v3, HLT_HT300_AlphaT0p53_v2, HLT_HT350_v4, HLT_HT350_AlphaT0p51_v3, HLT_HT350_AlphaT0p53_v3, HLT_HT400_v4, HLT_HT400_AlphaT0p51_v3, HLT_HT450_v4, HLT_HT500_v4, HLT_HT550_v4, HLT_PFMHT150_v5, HLT_MET100_v3, HLT_MET120_v3, HLT_MET200_v3, HLT_R014_MR150_v1, HLT_R014_MR150_CentralJet40_BTagIP_v1, HLT_R014_MR450_CentralJet40_BTagIP_v1, HLT_R014_MR600_v1, HLT_R020_MR150_v1, HLT_R020_MR350_CentralJet40_BTagIP_v1, HLT_R020_MR500_v1, HLT_R020_MR550_v1, HLT_R025_MR150_v1, HLT_R025_MR250_CentralJet40_BTagIP_v1, HLT_R025_MR400_v1, HLT_R025_MR450_v1, HLT_R033_MR300_v1, HLT_R033_MR350_v1, HLT_R038_MR200_v1, HLT_R038_MR250_v1, HLT_L1SingleMuOpen_v2, HLT_L1SingleMuOpen_DT_v2, HLT_L1SingleMu10_v2, HLT_L1SingleMu20_v2, HLT_L1DoubleMu0_v2, HLT_L2Mu10_v3, HLT_L2Mu20_v3, HLT_L2DoubleMu0_v4, HLT_Mu3_v5, HLT_Mu5_v5, HLT_Mu8_v3, HLT_Mu12_v3, HLT_Mu15_v4, HLT_Mu20_v3, HLT_Mu24_v3, HLT_Mu30_v3, HLT_Mu40_v1, HLT_IsoMu12_v3, HLT_IsoMu15_v7, HLT_IsoMu17_v7, HLT_IsoMu24_v3, HLT_IsoMu30_v3, HLT_L2DoubleMu23_NoVertex_v3, HLT_DoubleMu3_v5, HLT_DoubleMu6_v3, HLT_DoubleMu7_v3, HLT_DoubleMu2_Bs_v3, HLT_DoubleMu4_Acoplanarity03_v3, HLT_Dimuon0_Jpsi_v1, HLT_Dimuon0_Upsilon_v1, HLT_Dimuon4_Bs_Barrel_v1, HLT_Dimuon5_Upsilon_Barrel_v1, HLT_Dimuon6_Bs_v1, HLT_Dimuon7_LowMass_Displaced_v1, HLT_Dimuon7_Jpsi_Displaced_v1, HLT_Dimuon7_Jpsi_X_Barrel_v1, HLT_Dimuon7_PsiPrime_v1, HLT_Dimuon10_Jpsi_Barrel_v1, HLT_Dimuon0_Jpsi_Muon_v1, HLT_Dimuon0_Upsilon_Muon_v1, HLT_Mu13_Mu8_v2, HLT_TripleMu5_v4, HLT_Mu5_L2Mu2_Jpsi_v4, HLT_Photon20_CaloIdVL_IsoL_v3, HLT_Photon20_R9Id_Photon18_R9Id_v4, HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4, HLT_Photon26_Photon18_v4, HLT_Photon26_IsoVL_Photon18_v4, HLT_Photon26_IsoVL_Photon18_IsoVL_v4, HLT_Photon26_CaloIdL_IsoVL_Photon18_v4, HLT_Photon26_CaloIdL_IsoVL_Photon18_R9Id_v3, HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v4, HLT_Photon26_R9Id_Photon18_CaloIdL_IsoVL_v3, HLT_Photon26_R9Id_Photon18_R9Id_v1, HLT_Photon30_CaloIdVL_v4, HLT_Photon30_CaloIdVL_IsoL_v4, HLT_Photon32_CaloIdL_Photon26_CaloIdL_v4, HLT_Photon36_IsoVL_Photon22_v1, HLT_Photon36_CaloIdL_Photon22_CaloIdL_v3, HLT_Photon36_CaloIdL_IsoVL_Photon22_v1, HLT_Photon40_CaloIdL_Photon28_CaloIdL_v1, HLT_Photon50_CaloIdVL_IsoL_v3, HLT_Photon70_CaloIdL_HT300_v4, HLT_Photon70_CaloIdL_HT350_v3, HLT_Photon70_CaloIdL_MHT50_v4, HLT_Photon70_CaloIdL_MHT70_v3, HLT_Photon75_CaloIdVL_v4, HLT_Photon75_CaloIdVL_IsoL_v4, HLT_Photon90_CaloIdVL_v1, HLT_Photon90_CaloIdVL_IsoL_v1, HLT_Photon125_v1, HLT_Photon200_NoHE_v1, HLT_DoublePhoton33_v4, HLT_DoublePhoton33_HEVT_v1, HLT_DoublePhoton50_v1, HLT_DoublePhoton60_v1, HLT_DoublePhoton5_IsoVL_CEP_v3, HLT_L1SingleEG5_v2, HLT_L1SingleEG12_v2, HLT_Ele8_v4, HLT_Ele8_CaloIdL_CaloIsoVL_v4, HLT_Ele8_CaloIdL_TrkIdVL_v4, HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4, HLT_Ele17_CaloIdL_CaloIsoVL_v4, HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4, HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4, HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v1, HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v4, HLT_Ele27_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1, HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4, HLT_Ele32_CaloIdVL_CaloIsoVL_TrkIdVL_TrkIsoVL_v1, HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3, HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1, HLT_Ele45_CaloIdVT_TrkIdT_v4, HLT_Ele52_CaloIdVT_TrkIdT_v1, HLT_DoubleEle8_CaloIdL_TrkIdVL_v1, HLT_DoubleEle33_v1, HLT_DoubleEle33_CaloIdL_v1, HLT_IsoPFTau35_Trk20_MET45_v5, HLT_DoubleIsoPFTau25_Trk5_eta2p1_v1, HLT_BTagMu_DiJet20_Mu5_v4, HLT_BTagMu_DiJet40_Mu5_v4, HLT_BTagMu_DiJet70_Mu5_v4, HLT_BTagMu_DiJet110_Mu5_v4, HLT_Mu5_Ele8_CaloIdL_TrkIdVL_Ele8_v5, HLT_Mu8_R005_MR200_v1, HLT_Mu8_R020_MR200_v1, HLT_Mu8_R025_MR200_v1, HLT_HT250_Mu5_PFMHT35_v3, HLT_HT250_Mu15_PFMHT20_v1, HLT_HT300_Mu5_PFMHT40_v1, HLT_HT350_Mu5_PFMHT45_v1, HLT_Mu3_DiJet30_v1, HLT_Mu3_TriJet30_v1, HLT_Mu3_QuadJet30_v1, HLT_Mu8_Ele17_CaloIdL_v4, HLT_Mu8_Photon20_CaloIdVT_IsoT_v4, HLT_Mu8_Jet40_v5, HLT_Mu15_Photon20_CaloIdL_v5, HLT_Mu15_DoublePhoton15_CaloIdL_v5, HLT_Mu15_LooseIsoPFTau15_v1, HLT_Mu17_CentralJet30_v5, HLT_Mu17_DiCentralJet30_v5, HLT_Mu17_TriCentralJet30_v5, HLT_Mu17_Ele8_CaloIdL_v4, HLT_Mu12_CentralJet30_BTagIP_v3, HLT_Mu17_CentralJet30_BTagIP_v3, HLT_Mu15_HT200_v1, HLT_Mu20_HT200_v1, HLT_IsoMu15_LooseIsoPFTau15_v1, HLT_IsoMu17_CentralJet30_BTagIP_v3, HLT_DoubleMu3_HT150_v2, HLT_DoubleMu3_HT200_v5, HLT_DoubleMu5_Ele8_v5, HLT_Photon40_R005_MR150_v1, HLT_Photon40_R014_MR450_v1, HLT_Photon40_R020_MR300_v1, HLT_Photon40_R025_MR200_v1, HLT_Photon40_R038_MR150_v1, HLT_DoublePhoton40_MR150_v1, HLT_DoublePhoton40_R014_MR150_v1, HLT_Ele8_CaloIdT_TrkIdT_DiJet30_v1, HLT_Ele8_CaloIdT_TrkIdT_TriJet30_v1, HLT_Ele8_CaloIdT_TrkIdT_QuadJet30_v1, HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4, HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT200_v4, HLT_Ele15_CaloIdT_CaloIsoVL_TrkIdT_TrkIsoVL_HT250_v4, HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v5, HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v5, HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v5, HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta3_v2, HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_Jet35_Jet25_Deta2_v2, HLT_Ele15_CaloIdVT_TrkIdT_Jet35_Jet25_Deta2_v2, HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_v4, HLT_Ele25_CaloIdVT_TrkIdT_DiCentralJet30_v3, HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30_v3, HLT_Ele25_CaloIdVT_TrkIdT_CentralJet30_BTagIP_v3, HLT_Ele17_CaloIdVT_TrkIdT_CentralJet30_CentralJet25_v1, HLT_Ele17_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT15_v1, HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v1, HLT_DoubleEle8_CaloIdL_TrkIdVL_HT150_v2, HLT_DoubleEle8_CaloIdT_TrkIdVL_HT150_v2, HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v4, HLT_TripleEle10_CaloIdL_TrkIdVL_v4, HLT_PixelTracks_Multiplicity80_v3, HLT_PixelTracks_Multiplicity100_v3, HLT_BeamGas_HF_v4, HLT_BeamGas_BSC_v3, HLT_BeamHalo_v3, HLT_L1Tech_BSC_minBias_threshold1_v4, HLT_L1Tech_BSC_halo_v4, HLT_L1_PreCollisions_v2, HLT_L1_Interbunch_BSC_v2, HLT_GlobalRunHPDNoise_v3, HLT_L1Tech_HBHEHO_totalOR_v2, HLT_ZeroBias_v3, HLT_Physics_v1, HLT_Physics_NanoDST_v1, HLT_L1TrackerCosmics_v3, HLT_RegionalCosmicTracking_v3, HLT_L3MuonsCosmicTracking_v3, HLT_LogMonitor_v1, HLT_L1ETM30_v2, HLT_L1DoubleJet36Central_v2, HLT_L1MultiJet_v2, HLTriggerFinalPath ))

