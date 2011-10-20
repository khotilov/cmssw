# /dev/CMSSW_4_2_0/HIon/V324 (CMSSW_4_2_0_HLT34)

import FWCore.ParameterSet.Config as cms

process = cms.Process( "HLT" )

process.HLTConfigVersion = cms.PSet(
  tableName = cms.string('/dev/CMSSW_4_2_0/HIon/V324')
)

process.streams = cms.PSet( 
  A = cms.vstring( 'HIDiMuon',
    'HIHighPt',
    'HIMinBiasUPC',
    'LogMonitor' ),
  Calibration = cms.vstring( 'TestEnablesEcalHcalDT' ),
  DQM = cms.vstring( 'OnlineMonitor',
    'OnlineMonitorHI' ),
  EcalCalibration = cms.vstring( 'EcalLaser' ),
  Express = cms.vstring( 'HIExpressPhysics' ),
  HLTDQM = cms.vstring( 'OnlineHltMonitor',
    'OnlineHltMonitorHI' ),
  HLTDQMResults = cms.vstring( 'OnlineHltResults' ),
  HLTMON = cms.vstring( 'OfflineMonitor',
    'OfflineMonitorHI' ),
  TrackerCalibration = cms.vstring( 'TestEnablesTracker' )
)
process.datasets = cms.PSet( 
  EcalLaser = cms.vstring( 'HLT_EcalCalibration_v2' ),
  HIDiMuon = cms.vstring( 'HLT_HIL1DoubleMuOpen_v1',
    'HLT_HIL1SingleMu3_v1',
    'HLT_HIL1SingleMu5_v1',
    'HLT_HIL1SingleMu7_v1',
    'HLT_HIL2DoubleMu0_v1',
    'HLT_HIL2DoubleMu3_v1',
    'HLT_HIL2Mu20_v1',
    'HLT_HIL2Mu3_v1',
    'HLT_HIL2Mu5Tight_v1' ),
  HIExpressPhysics = cms.vstring( 'HLT_HIL2DoubleMu3_v1',
    'HLT_HIMinBiasHfOrBSC_v1' ),
  HIHighPt = cms.vstring( 'HLT_HIDoublePhoton10_v1',
    'HLT_HIDoublePhoton15_v1',
    'HLT_HIDoublePhoton20_v1',
    'HLT_HIFullTrack12_L1Central_v1',
    'HLT_HIFullTrack14_L1Central_v1',
    'HLT_HIFullTrack20_L1Central_v1',
    'HLT_HIFullTrack25_L1Central_v1',
    'HLT_HIPhoton10_Photon15_v1',
    'HLT_HIPhoton15_Photon20_v1',
    'HLT_HISinglePhoton15_v1',
    'HLT_HISinglePhoton20_v1',
    'HLT_HISinglePhoton30_v1',
    'HLT_HISinglePhoton40_v1',
    'HLT_HIStoppedHSCP35_v1' ),
  HIMinBiasUPC = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
    'HLT_HIActivityHF_Single3_v1',
    'HLT_HIBptxXOR_v1',
    'HLT_HICentralityVeto_v1',
    'HLT_HIClusterVertexCompatibility_v1',
    'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
    'HLT_HIMinBiasBSC_OR_v1',
    'HLT_HIMinBiasBSC_v1',
    'HLT_HIMinBiasHF_v1',
    'HLT_HIMinBiasHfOrBSC_v1',
    'HLT_HIMinBiasHf_OR_v1',
    'HLT_HIMinBiasPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
    'HLT_HIMinBiasZDC_Calo_v1',
    'HLT_HIMinBiasZDC_Scint_v1',
    'HLT_HIRandom_v1',
    'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
    'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
    'HLT_HIUpcEcal_v1',
    'HLT_HIUpcMu_v1',
    'HLT_HIZeroBiasPixel_SingleTrack_v1',
    'HLT_HIZeroBiasXOR_v1',
    'HLT_HIZeroBias_v1' ),
  LogMonitor = cms.vstring( 'HLT_LogMonitor_v1' ),
  OfflineMonitor = cms.vstring( 'HLT_HcalCalibration_v2',
    'HLT_LogMonitor_v1' ),
  OfflineMonitorHI = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
    'HLT_HIActivityHF_Single3_v1',
    'HLT_HIBptxXOR_v1',
    'HLT_HICentralityVeto_v1',
    'HLT_HIClusterVertexCompatibility_v1',
    'HLT_HIDoublePhoton10_v1',
    'HLT_HIDoublePhoton15_v1',
    'HLT_HIDoublePhoton20_v1',
    'HLT_HIFullTrack12_L1Central_v1',
    'HLT_HIFullTrack14_L1Central_v1',
    'HLT_HIFullTrack20_L1Central_v1',
    'HLT_HIFullTrack25_L1Central_v1',
    'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
    'HLT_HIL1DoubleMuOpen_v1',
    'HLT_HIL1SingleMu3_v1',
    'HLT_HIL1SingleMu5_v1',
    'HLT_HIL1SingleMu7_v1',
    'HLT_HIL2DoubleMu0_v1',
    'HLT_HIL2DoubleMu3_v1',
    'HLT_HIL2Mu20_v1',
    'HLT_HIL2Mu3_v1',
    'HLT_HIL2Mu5Tight_v1',
    'HLT_HIMinBiasBSC_OR_v1',
    'HLT_HIMinBiasBSC_v1',
    'HLT_HIMinBiasHF_v1',
    'HLT_HIMinBiasHfOrBSC_v1',
    'HLT_HIMinBiasHf_OR_v1',
    'HLT_HIMinBiasPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
    'HLT_HIMinBiasZDC_Calo_v1',
    'HLT_HIMinBiasZDC_Scint_v1',
    'HLT_HIPhoton10_Photon15_v1',
    'HLT_HIPhoton15_Photon20_v1',
    'HLT_HIRandom_v1',
    'HLT_HISinglePhoton15_v1',
    'HLT_HISinglePhoton20_v1',
    'HLT_HISinglePhoton30_v1',
    'HLT_HISinglePhoton40_v1',
    'HLT_HIStoppedHSCP35_v1',
    'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
    'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
    'HLT_HIUpcEcal_v1',
    'HLT_HIUpcMu_v1',
    'HLT_HIZeroBiasPixel_SingleTrack_v1',
    'HLT_HIZeroBiasXOR_v1',
    'HLT_HIZeroBias_v1' ),
  OnlineHltMonitor = cms.vstring( 'HLT_LogMonitor_v1' ),
  OnlineHltMonitorHI = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
    'HLT_HIActivityHF_Single3_v1',
    'HLT_HIBptxXOR_v1',
    'HLT_HICentralityVeto_v1',
    'HLT_HIClusterVertexCompatibility_v1',
    'HLT_HIDoublePhoton10_v1',
    'HLT_HIDoublePhoton15_v1',
    'HLT_HIDoublePhoton20_v1',
    'HLT_HIFullTrack12_L1Central_v1',
    'HLT_HIFullTrack14_L1Central_v1',
    'HLT_HIFullTrack20_L1Central_v1',
    'HLT_HIFullTrack25_L1Central_v1',
    'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
    'HLT_HIL1DoubleMuOpen_v1',
    'HLT_HIL1SingleMu3_v1',
    'HLT_HIL1SingleMu5_v1',
    'HLT_HIL1SingleMu7_v1',
    'HLT_HIL2DoubleMu0_v1',
    'HLT_HIL2DoubleMu3_v1',
    'HLT_HIL2Mu20_v1',
    'HLT_HIL2Mu3_v1',
    'HLT_HIL2Mu5Tight_v1',
    'HLT_HIMinBiasBSC_OR_v1',
    'HLT_HIMinBiasBSC_v1',
    'HLT_HIMinBiasHF_v1',
    'HLT_HIMinBiasHfOrBSC_v1',
    'HLT_HIMinBiasHf_OR_v1',
    'HLT_HIMinBiasPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
    'HLT_HIMinBiasZDC_Calo_v1',
    'HLT_HIMinBiasZDC_Scint_v1',
    'HLT_HIPhoton10_Photon15_v1',
    'HLT_HIPhoton15_Photon20_v1',
    'HLT_HIRandom_v1',
    'HLT_HISinglePhoton15_v1',
    'HLT_HISinglePhoton20_v1',
    'HLT_HISinglePhoton30_v1',
    'HLT_HISinglePhoton40_v1',
    'HLT_HIStoppedHSCP35_v1',
    'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
    'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
    'HLT_HIUpcEcal_v1',
    'HLT_HIUpcMu_v1',
    'HLT_HIZeroBiasPixel_SingleTrack_v1',
    'HLT_HIZeroBiasXOR_v1',
    'HLT_HIZeroBias_v1' ),
  OnlineHltResults = cms.vstring( 'HLTriggerFinalPath' ),
  OnlineMonitor = cms.vstring( 'HLT_DTCalibration_v1',
    'HLT_EcalCalibration_v2',
    'HLT_HcalCalibration_v2',
    'HLT_LogMonitor_v1',
    'HLT_TrackerCalibration_v2' ),
  OnlineMonitorHI = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
    'HLT_HIActivityHF_Single3_v1',
    'HLT_HIBptxXOR_v1',
    'HLT_HICentralityVeto_v1',
    'HLT_HIClusterVertexCompatibility_v1',
    'HLT_HIDoublePhoton10_v1',
    'HLT_HIDoublePhoton15_v1',
    'HLT_HIDoublePhoton20_v1',
    'HLT_HIFullTrack12_L1Central_v1',
    'HLT_HIFullTrack14_L1Central_v1',
    'HLT_HIFullTrack20_L1Central_v1',
    'HLT_HIFullTrack25_L1Central_v1',
    'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
    'HLT_HIL1DoubleMuOpen_v1',
    'HLT_HIL1SingleMu3_v1',
    'HLT_HIL1SingleMu5_v1',
    'HLT_HIL1SingleMu7_v1',
    'HLT_HIL2DoubleMu0_v1',
    'HLT_HIL2DoubleMu3_v1',
    'HLT_HIL2Mu20_v1',
    'HLT_HIL2Mu3_v1',
    'HLT_HIL2Mu5Tight_v1',
    'HLT_HIMinBiasBSC_OR_v1',
    'HLT_HIMinBiasBSC_v1',
    'HLT_HIMinBiasHF_v1',
    'HLT_HIMinBiasHfOrBSC_v1',
    'HLT_HIMinBiasHf_OR_v1',
    'HLT_HIMinBiasPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
    'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
    'HLT_HIMinBiasZDC_Calo_v1',
    'HLT_HIMinBiasZDC_Scint_v1',
    'HLT_HIPhoton10_Photon15_v1',
    'HLT_HIPhoton15_Photon20_v1',
    'HLT_HIRandom_v1',
    'HLT_HISinglePhoton15_v1',
    'HLT_HISinglePhoton20_v1',
    'HLT_HISinglePhoton30_v1',
    'HLT_HISinglePhoton40_v1',
    'HLT_HIStoppedHSCP35_v1',
    'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
    'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
    'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
    'HLT_HIUpcEcal_v1',
    'HLT_HIUpcMu_v1',
    'HLT_HIZeroBiasPixel_SingleTrack_v1',
    'HLT_HIZeroBiasXOR_v1',
    'HLT_HIZeroBias_v1' ),
  TestEnablesEcalHcalDT = cms.vstring( 'HLT_DTCalibration_v1',
    'HLT_EcalCalibration_v2',
    'HLT_HcalCalibration_v2' ),
  TestEnablesTracker = cms.vstring( 'HLT_TrackerCalibration_v2' )
)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring( '/store/data/Run2011A/MinimumBias/RAW/v1/000/165/205/6C8BA6D0-F680-E011-B467-003048F118AC.root' )
)

process.GlobalTag = cms.ESSource( "PoolDBESSource",
    appendToDataLabel = cms.string( "" ),
    timetype = cms.string( "runnumber" ),
    connect = cms.string( "frontier://(proxyurl=http://localhost:3128)(serverurl=http://localhost:8000/FrontierOnProd)(serverurl=http://localhost:8000/FrontierOnProd)(retrieve-ziplevel=0)/CMS_COND_31X_GLOBALTAG" ),
    DumpStat = cms.untracked.bool( False ),
    BlobStreamerName = cms.untracked.string( "TBufferBlobStreamingService" ),
    globaltag = cms.string( "GR_H_V22::All" ),
    DBParameters = cms.PSet( 
      authenticationPath = cms.untracked.string( "." ),
      connectionRetrialTimeOut = cms.untracked.int32( 60 ),
      idleConnectionCleanupPeriod = cms.untracked.int32( 10 ),
      messageLevel = cms.untracked.int32( 0 ),
      enablePoolAutomaticCleanUp = cms.untracked.bool( False ),
      enableConnectionSharing = cms.untracked.bool( True ),
      enableReadOnlySessionOnUpdateConnection = cms.untracked.bool( False ),
      connectionTimeOut = cms.untracked.int32( 0 ),
      connectionRetrialPeriod = cms.untracked.int32( 10 )
    ),
    toGet = cms.VPSet( 
    ),
    RefreshEachRun = cms.untracked.bool( True )
)
process.HepPDTESSource = cms.ESSource( "HepPDTESSource",
    pdtFileName = cms.FileInPath( "SimGeneral/HepPDTESSource/data/pythiaparticle.tbl" ),
    appendToDataLabel = cms.string( "" )
)
process.eegeom = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "EcalMappingRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    appendToDataLabel = cms.string( "" ),
    firstValid = cms.vuint32( 1 )
)
process.es_hardcode = cms.ESSource( "HcalHardcodeCalibrations",
    toGet = cms.untracked.vstring( 'GainWidths' ),
    appendToDataLabel = cms.string( "" )
)
process.hltESSAK5CaloL1L2L3 = cms.ESSource( "JetCorrectionServiceChain",
    appendToDataLabel = cms.string( "" ),
    correctors = cms.vstring( 'hltESSL1FastJetCorrectionService',
      'hltESSL2RelativeCorrectionService',
      'hltESSL3AbsoluteCorrectionService' ),
    label = cms.string( "hltESSAK5CaloL1L2L3" )
)
process.hltESSAK5CaloL2L3 = cms.ESSource( "JetCorrectionServiceChain",
    appendToDataLabel = cms.string( "" ),
    correctors = cms.vstring( 'hltESSL2RelativeCorrectionService',
      'hltESSL3AbsoluteCorrectionService' ),
    label = cms.string( "hltESSAK5CaloL2L3" )
)
process.hltESSBTagRecord = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "JetTagComputerRecord" ),
    iovIsRunNotTime = cms.bool( True ),
    appendToDataLabel = cms.string( "" ),
    firstValid = cms.vuint32( 1 )
)
process.hltESSEcalSeverityLevel = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "EcalSeverityLevelAlgoRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    appendToDataLabel = cms.string( "" ),
    firstValid = cms.vuint32( 1 )
)
process.hltESSHcalSeverityLevel = cms.ESSource( "EmptyESSource",
    recordName = cms.string( "HcalSeverityLevelComputerRcd" ),
    iovIsRunNotTime = cms.bool( True ),
    appendToDataLabel = cms.string( "" ),
    firstValid = cms.vuint32( 1 )
)
process.hltESSL1FastJetCorrectionService = cms.ESSource( "L1FastjetCorrectionService",
    era = cms.string( "Jec10V1" ),
    level = cms.string( "L1FastJet" ),
    algorithm = cms.string( "AK5Calo" ),
    section = cms.string( "" ),
    srcRho = cms.InputTag( 'hltKT6CaloJets','rho' ),
    useCondDB = cms.untracked.bool( True )
)
process.hltESSL2RelativeCorrectionService = cms.ESSource( "LXXXCorrectionService",
    appendToDataLabel = cms.string( "" ),
    level = cms.string( "L2Relative" ),
    algorithm = cms.string( "AK5Calo" ),
    section = cms.string( "" ),
    era = cms.string( "" ),
    useCondDB = cms.untracked.bool( True )
)
process.hltESSL3AbsoluteCorrectionService = cms.ESSource( "LXXXCorrectionService",
    appendToDataLabel = cms.string( "" ),
    level = cms.string( "L3Absolute" ),
    algorithm = cms.string( "AK5Calo" ),
    section = cms.string( "" ),
    era = cms.string( "" ),
    useCondDB = cms.untracked.bool( True )
)
process.magfield = cms.ESSource( "XMLIdealGeometryESSource",
    rootNodeName = cms.string( "cmsMagneticField:MAGF" ),
    appendToDataLabel = cms.string( "" ),
    geomXMLFiles = cms.vstring( 'Geometry/CMSCommonData/data/normal/cmsextent.xml',
      'Geometry/CMSCommonData/data/cms.xml',
      'Geometry/CMSCommonData/data/cmsMagneticField.xml',
      'MagneticField/GeomBuilder/data/MagneticFieldVolumes_1103l.xml',
      'MagneticField/GeomBuilder/data/MagneticFieldParameters_07_2pi.xml' )
)

process.AnyDirectionAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  ComponentName = cms.string( "AnyDirectionAnalyticalPropagator" ),
  PropagationDirection = cms.string( "anyDirection" ),
  MaxDPhi = cms.double( 1.6 ),
  appendToDataLabel = cms.string( "" )
)
process.AutoMagneticFieldESProducer = cms.ESProducer( "AutoMagneticFieldESProducer",
  label = cms.untracked.string( "" ),
  valueOverride = cms.int32( -1 ),
  appendToDataLabel = cms.string( "" ),
  nominalCurrents = cms.untracked.vint32( -1, 0, 9558, 14416, 16819, 18268, 19262 ),
  mapLabels = cms.untracked.vstring( '090322_3_8t',
    '0t',
    '071212_2t',
    '071212_3t',
    '071212_3_5t',
    '090322_3_8t',
    '071212_4t' )
)
process.CSCGeometryESModule = cms.ESProducer( "CSCGeometryESModule",
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  useRealWireGeometry = cms.bool( True ),
  useOnlyWiresInME1a = cms.bool( False ),
  useGangedStripsInME1a = cms.bool( True ),
  useCentreTIOffsets = cms.bool( False ),
  useDDD = cms.bool( False ),
  applyAlignment = cms.bool( True )
)
process.CaloGeometryBuilder = cms.ESProducer( "CaloGeometryBuilder",
  appendToDataLabel = cms.string( "" ),
  SelectedCalos = cms.vstring( 'HCAL',
    'ZDC',
    'EcalBarrel',
    'EcalEndcap',
    'EcalPreshower',
    'TOWER' )
)
process.CaloTopologyBuilder = cms.ESProducer( "CaloTopologyBuilder",
  appendToDataLabel = cms.string( "" )
)
process.CaloTowerConstituentsMapBuilder = cms.ESProducer( "CaloTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string( "Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz" ),
  appendToDataLabel = cms.string( "" )
)
process.CaloTowerGeometryFromDBEP = cms.ESProducer( "CaloTowerGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
process.CastorDbProducer = cms.ESProducer( "CastorDbProducer",
  appendToDataLabel = cms.string( "" )
)
process.ClusterShapeHitFilterESProducer = cms.ESProducer( "ClusterShapeHitFilterESProducer",
  ComponentName = cms.string( "ClusterShapeHitFilter" ),
  appendToDataLabel = cms.string( "" )
)
process.DTGeometryESModule = cms.ESProducer( "DTGeometryESModule",
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  fromDDD = cms.bool( False ),
  applyAlignment = cms.bool( True )
)
process.EcalBarrelGeometryFromDBEP = cms.ESProducer( "EcalBarrelGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
process.EcalElectronicsMappingBuilder = cms.ESProducer( "EcalElectronicsMappingBuilder",
  appendToDataLabel = cms.string( "" )
)
process.EcalEndcapGeometryFromDBEP = cms.ESProducer( "EcalEndcapGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
process.EcalLaserCorrectionService = cms.ESProducer( "EcalLaserCorrectionService",
  appendToDataLabel = cms.string( "" )
)
process.EcalPreshowerGeometryFromDBEP = cms.ESProducer( "EcalPreshowerGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True )
)
process.EcalUnpackerWorkerESProducer = cms.ESProducer( "EcalUnpackerWorkerESProducer",
  ComponentName = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  DCCDataUnpacker = cms.PSet( 
    orderedDCCIdList = cms.vint32( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54 ),
    tccUnpacking = cms.bool( False ),
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
    flagsMapDBReco = cms.vint32( 0, 0, 0, 0, 4, -1, -1, -1, 4, 4, 7, 7, 7, 8, 9 ),
    Type = cms.string( "EcalRecHitWorkerSimple" ),
    killDeadChannels = cms.bool( True ),
    ChannelStatusToBeExcluded = cms.vint32( 10, 11, 12, 13, 14 ),
    laserCorrection = cms.bool( False )
  )
)
process.HcalGeometryFromDBEP = cms.ESProducer( "HcalGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
process.HcalTopologyIdealEP = cms.ESProducer( "HcalTopologyIdealEP",
  appendToDataLabel = cms.string( "" )
)
process.MaterialPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "PropagatorWithMaterial" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Mass = cms.double( 0.105 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
process.MaterialPropagatorForHI = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "PropagatorWithMaterialForHI" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Mass = cms.double( 0.139 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
process.MuonNumberingInitialization = cms.ESProducer( "MuonNumberingInitialization",
  appendToDataLabel = cms.string( "" )
)
process.OppositeMaterialPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "PropagatorWithMaterialOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Mass = cms.double( 0.105 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
process.OppositeMaterialPropagatorForHI = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "PropagatorWithMaterialOppositeForHI" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Mass = cms.double( 0.139 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( False ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
process.RPCGeometryESModule = cms.ESProducer( "RPCGeometryESModule",
  useDDD = cms.untracked.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.SiStripGainESProducer = cms.ESProducer( "SiStripGainESProducer",
  AutomaticNormalization = cms.bool( False ),
  NormalizationFactor = cms.double( 1.0 ),
  printDebug = cms.untracked.bool( False ),
  APVGain = cms.VPSet( 
    cms.PSet(  Record = cms.string( "SiStripApvGainRcd" ),
      NormalizationFactor = cms.untracked.double( 1.0 ),
      Label = cms.untracked.string( "" )
    ),
    cms.PSet(  Record = cms.string( "SiStripApvGain2Rcd" ),
      NormalizationFactor = cms.untracked.double( 1.0 ),
      Label = cms.untracked.string( "" )
    )
  )
)
process.SiStripQualityESProducer = cms.ESProducer( "SiStripQualityESProducer",
  appendToDataLabel = cms.string( "" ),
  PrintDebugOutput = cms.bool( False ),
  ThresholdForReducedGranularity = cms.double( 0.3 ),
  UseEmptyRunInfo = cms.bool( False ),
  ReduceGranularity = cms.bool( False ),
  ListOfRecordToMerge = cms.VPSet( 
    cms.PSet(  record = cms.string( "SiStripDetVOffRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripDetCablingRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadChannelRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadFiberRcd" ),
      tag = cms.string( "" )
    ),
    cms.PSet(  record = cms.string( "SiStripBadModuleRcd" ),
      tag = cms.string( "" )
    )
  )
)
process.SiStripRecHitMatcherESProducer = cms.ESProducer( "SiStripRecHitMatcherESProducer",
  ComponentName = cms.string( "StandardMatcher" ),
  NSigmaInside = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
process.SlaveField0 = cms.ESProducer( "UniformMagneticFieldESProducer",
  ZFieldInTesla = cms.double( 0.0 ),
  label = cms.untracked.string( "slave_0" ),
  appendToDataLabel = cms.string( "" )
)
process.SlaveField20 = cms.ESProducer( "ParametrizedMagneticFieldProducer",
  label = cms.untracked.string( "slave_20" ),
  version = cms.string( "OAE_1103l_071212" ),
  appendToDataLabel = cms.string( "" ),
  parameters = cms.PSet(  BValue = cms.string( "2_0T" ) )
)
process.SlaveField30 = cms.ESProducer( "ParametrizedMagneticFieldProducer",
  label = cms.untracked.string( "slave_30" ),
  version = cms.string( "OAE_1103l_071212" ),
  appendToDataLabel = cms.string( "" ),
  parameters = cms.PSet(  BValue = cms.string( "3_0T" ) )
)
process.SlaveField35 = cms.ESProducer( "ParametrizedMagneticFieldProducer",
  label = cms.untracked.string( "slave_35" ),
  version = cms.string( "OAE_1103l_071212" ),
  appendToDataLabel = cms.string( "" ),
  parameters = cms.PSet(  BValue = cms.string( "3_5T" ) )
)
process.SlaveField38 = cms.ESProducer( "ParametrizedMagneticFieldProducer",
  label = cms.untracked.string( "slave_38" ),
  version = cms.string( "OAE_1103l_071212" ),
  appendToDataLabel = cms.string( "" ),
  parameters = cms.PSet(  BValue = cms.string( "3_8T" ) )
)
process.SlaveField40 = cms.ESProducer( "ParametrizedMagneticFieldProducer",
  label = cms.untracked.string( "slave_40" ),
  version = cms.string( "OAE_1103l_071212" ),
  appendToDataLabel = cms.string( "" ),
  parameters = cms.PSet(  BValue = cms.string( "4_0T" ) )
)
process.SteppingHelixPropagatorAny = cms.ESProducer( "SteppingHelixPropagatorESProducer",
  ComponentName = cms.string( "SteppingHelixPropagatorAny" ),
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
  useTuningForL2Speed = cms.bool( False ),
  useEndcapShiftsInZ = cms.bool( False ),
  endcapShiftInZPos = cms.double( 0.0 ),
  endcapShiftInZNeg = cms.double( 0.0 ),
  appendToDataLabel = cms.string( "" )
)
process.StripCPEfromTrackAngleESProducer = cms.ESProducer( "StripCPEESProducer",
  ComponentName = cms.string( "StripCPEfromTrackAngle" ),
  appendToDataLabel = cms.string( "" ),
  TanDiffusionAngle = cms.double( 0.01 ),
  ThicknessRelativeUncertainty = cms.double( 0.02 ),
  NoiseThreshold = cms.double( 2.3 ),
  MaybeNoiseThreshold = cms.double( 3.5 ),
  UncertaintyScaling = cms.double( 1.42 ),
  MinimumUncertainty = cms.double( 0.01 )
)
process.TrackerDigiGeometryESModule = cms.ESProducer( "TrackerDigiGeometryESModule",
  alignmentsLabel = cms.string( "" ),
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( True ),
  fromDDD = cms.bool( False )
)
process.TrackerGeometricDetESModule = cms.ESProducer( "TrackerGeometricDetESModule",
  fromDDD = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.TransientTrackBuilderESProducer = cms.ESProducer( "TransientTrackBuilderESProducer",
  ComponentName = cms.string( "TransientTrackBuilder" ),
  appendToDataLabel = cms.string( "" )
)
process.VBF0 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "0t" ),
  version = cms.string( "grid_1103l_071212_2t" ),
  overrideMasterSector = cms.bool( True ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_0" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32(  ),
  scalingFactors = cms.vdouble(  ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.VBF20 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "071212_2t" ),
  version = cms.string( "grid_1103l_071212_2t" ),
  overrideMasterSector = cms.bool( True ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_20" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32(  ),
  scalingFactors = cms.vdouble(  ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.VBF30 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "071212_3t" ),
  version = cms.string( "grid_1103l_071212_3t" ),
  overrideMasterSector = cms.bool( True ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_30" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32(  ),
  scalingFactors = cms.vdouble(  ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.VBF35 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "071212_3_5t" ),
  version = cms.string( "grid_1103l_071212_3_5t" ),
  overrideMasterSector = cms.bool( True ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_35" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32(  ),
  scalingFactors = cms.vdouble(  ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.VBF38 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "090322_3_8t" ),
  version = cms.string( "grid_1103l_090322_3_8t" ),
  overrideMasterSector = cms.bool( False ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_38" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32( 14100, 14200, 17600, 17800, 17900, 18100, 18300, 18400, 18600, 23100, 23300, 23400, 23600, 23800, 23900, 24100, 28600, 28800, 28900, 29100, 29300, 29400, 29600, 28609, 28809, 28909, 29109, 29309, 29409, 29609, 28610, 28810, 28910, 29110, 29310, 29410, 29610, 28611, 28811, 28911, 29111, 29311, 29411, 29611 ),
  scalingFactors = cms.vdouble( 1.0, 1.0, 0.994, 1.004, 1.004, 1.005, 1.004, 1.004, 0.994, 0.965, 0.958, 0.958, 0.953, 0.958, 0.958, 0.965, 0.918, 0.924, 0.924, 0.906, 0.924, 0.924, 0.918, 0.991, 0.998, 0.998, 0.978, 0.998, 0.998, 0.991, 0.991, 0.998, 0.998, 0.978, 0.998, 0.998, 0.991, 0.991, 0.998, 0.998, 0.978, 0.998, 0.998, 0.991 ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.VBF40 = cms.ESProducer( "VolumeBasedMagneticFieldESProducer",
  label = cms.untracked.string( "071212_4t" ),
  version = cms.string( "grid_1103l_071212_4t" ),
  overrideMasterSector = cms.bool( True ),
  useParametrizedTrackerField = cms.bool( True ),
  paramLabel = cms.string( "slave_40" ),
  appendToDataLabel = cms.string( "" ),
  scalingVolumes = cms.vint32(  ),
  scalingFactors = cms.vdouble(  ),
  findVolumeTolerance = cms.double( 0.0 ),
  cacheLastVolume = cms.untracked.bool( True )
)
process.XMLFromDBSource = cms.ESProducer( "XMLIdealGeometryESProducer",
  rootDDName = cms.string( "cms:OCMS" ),
  label = cms.string( "Extended" ),
  appendToDataLabel = cms.string( "" )
)
process.ZdcGeometryFromDBEP = cms.ESProducer( "ZdcGeometryFromDBEP",
  appendToDataLabel = cms.string( "" ),
  applyAlignment = cms.bool( False )
)
process.caloDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "CaloDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
process.cosmicsNavigationSchoolESProducer = cms.ESProducer( "NavigationSchoolESProducer",
  ComponentName = cms.string( "CosmicNavigationSchool" ),
  appendToDataLabel = cms.string( "" )
)
process.ecalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "EcalDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.02 ),
  nEta = cms.int32( 300 ),
  nPhi = cms.int32( 360 ),
  includeBadChambers = cms.bool( False )
)
process.ecalSeverityLevel = cms.ESProducer( "EcalSeverityLevelESProducer",
  appendToDataLabel = cms.string( "" ),
  flagMask = cms.vuint32( 1, 114, 896, 4, 49152, 3080 ),
  dbstatusMask = cms.vuint32( 1, 2046, 0, 0, 0, 64512 ),
  timeThresh = cms.double( 2.0 )
)
process.hcalDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HcalDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 70 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
process.hcalRecAlgos = cms.ESProducer( "HcalRecAlgoESProducer",
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
process.hcal_db_producer = cms.ESProducer( "HcalDbProducer",
  appendToDataLabel = cms.string( "" )
)
process.hltESPAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
  ComponentName = cms.string( "hltESPAnalyticalPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  MaxDPhi = cms.double( 1.6 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2EstimatorForRefit = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2EstimatorForRefit" ),
  MaxChi2 = cms.double( 100000.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator" ),
  MaxChi2 = cms.double( 30.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator16 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator16" ),
  MaxChi2 = cms.double( 16.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPChi2MeasurementEstimator9 = cms.ESProducer( "Chi2MeasurementEstimatorESProducer",
  ComponentName = cms.string( "hltESPChi2MeasurementEstimator9" ),
  MaxChi2 = cms.double( 9.0 ),
  nSigma = cms.double( 3.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPCkf3HitTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
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
process.hltESPCkf3HitTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
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
process.hltESPCkfTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
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
process.hltESPCkfTrajectoryBuilderForHI = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPCkfTrajectoryBuilderForHI" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterialForHI" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOppositeForHI" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTrackerForHI" ),
  trajectoryFilterName = cms.string( "hltESPCkfTrajectoryFilterForHI" ),
  maxCand = cms.int32( 5 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( False ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPCkfTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
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
process.hltESPCkfTrajectoryFilterForHI = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPCkfTrajectoryFilterForHI" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minimumNumberOfHits = cms.int32( 6 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( -1 ),
    maxConsecLostHits = cms.int32( 1 ),
    chargeSignificance = cms.double( -1.0 ),
    nSigmaMinPt = cms.double( 5.0 ),
    minPt = cms.double( 11.0 )
  )
)
process.hltESPDummyDetLayerGeometry = cms.ESProducer( "DetLayerGeometryESProducer",
  ComponentName = cms.string( "hltESPDummyDetLayerGeometry" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPESUnpackerWorker = cms.ESProducer( "ESUnpackerWorkerESProducer",
  ComponentName = cms.string( "hltESPESUnpackerWorker" ),
  appendToDataLabel = cms.string( "" ),
  DCCDataUnpacker = cms.PSet(  LookupTable = cms.FileInPath( "EventFilter/ESDigiToRaw/data/ES_lookup_table.dat" ) ),
  RHAlgo = cms.PSet( 
    ESRecoAlgo = cms.int32( 0 ),
    Type = cms.string( "ESRecHitWorker" )
  )
)
process.hltESPEcalRegionCablingESProducer = cms.ESProducer( "EcalRegionCablingESProducer",
  appendToDataLabel = cms.string( "" ),
  esMapping = cms.PSet(  LookupTable = cms.FileInPath( "EventFilter/ESDigiToRaw/data/ES_lookup_table.dat" ) )
)
process.hltESPEcalTrigTowerConstituentsMapBuilder = cms.ESProducer( "EcalTrigTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string( "Geometry/EcalMapping/data/EndCap_TTMap.txt" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPFastSteppingHelixPropagatorAny = cms.ESProducer( "SteppingHelixPropagatorESProducer",
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
process.hltESPFastSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
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
process.hltESPFittingSmootherIT = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPFittingSmootherIT" ),
  Fitter = cms.string( "hltESPTrajectoryFitterRK" ),
  Smoother = cms.string( "hltESPTrajectorySmootherRK" ),
  EstimateCut = cms.double( 10.0 ),
  LogPixelProbabilityCut = cms.double( -16.0 ),
  MinNumberOfHits = cms.int32( 3 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPFittingSmootherRK = cms.ESProducer( "KFFittingSmootherESProducer",
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
process.hltESPGlobalDetLayerGeometry = cms.ESProducer( "GlobalDetLayerGeometryESProducer",
  ComponentName = cms.string( "hltESPGlobalDetLayerGeometry" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPGlobalTrackingGeometryESProducer = cms.ESProducer( "GlobalTrackingGeometryESProducer",
  appendToDataLabel = cms.string( "" )
)
process.hltESPHIPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
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
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPHIPixelLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPHIPixelLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltHISiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPHITTRHBuilderWithoutRefit = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPHITTRHBuilderWithoutRefit" ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "Fake" ),
  Matcher = cms.string( "Fake" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFFittingSmoother = cms.ESProducer( "KFFittingSmootherESProducer",
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
process.hltESPKFFittingSmootherForL2Muon = cms.ESProducer( "KFFittingSmootherESProducer",
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
process.hltESPKFFittingSmootherWithOutliersRejectionAndRK = cms.ESProducer( "KFFittingSmootherESProducer",
  ComponentName = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
  Fitter = cms.string( "hltESPRKFitter" ),
  Smoother = cms.string( "hltESPRKSmoother" ),
  EstimateCut = cms.double( 20.0 ),
  LogPixelProbabilityCut = cms.double( -14.0 ),
  MinNumberOfHits = cms.int32( 3 ),
  RejectTracks = cms.bool( True ),
  BreakTrajWith2ConsecutiveMissing = cms.bool( True ),
  NoInvalidHitsBeginEnd = cms.bool( True ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitter" ),
  Propagator = cms.string( "PropagatorWithMaterial" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectoryFitterForL2Muon = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPKFTrajectoryFitterForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmoother" ),
  Propagator = cms.string( "PropagatorWithMaterial" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmootherForL2Muon = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForL2Muon" ),
  Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFTrajectorySmootherForMuonTrackLoader = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPKFTrajectorySmootherForMuonTrackLoader" ),
  Propagator = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 10.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPKFUpdator = cms.ESProducer( "KFUpdatorESProducer",
  ComponentName = cms.string( "hltESPKFUpdator" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPL3MuKFTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPL3MuKFTrajectoryFitter" ),
  Propagator = cms.string( "hltESPSmartPropagatorAny" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
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
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TID = cms.PSet( 
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
process.hltESPMeasurementTrackerForHI = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltESPMeasurementTrackerForHI" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  StripCPE = cms.string( "StripCPEfromTrackAngle" ),
  HitMatcher = cms.string( "StandardMatcher" ),
  Regional = cms.bool( False ),
  OnDemand = cms.bool( False ),
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
  pixelClusterProducer = cms.string( "hltHISiPixelClusters" ),
  skipClusters = cms.InputTag( "" ),
  stripClusterProducer = cms.string( "hltHISiStripClustersNonRegional" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripRawToDigi' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 2 ),
      maxBad = cms.uint32( 4 )
    ),
    TID = cms.PSet( 
      maxBad = cms.uint32( 4 ),
      maxConsecutiveBad = cms.uint32( 2 )
    ),
    TEC = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 2 ),
      maxBad = cms.uint32( 4 )
    ),
    TIB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 2 ),
      maxBad = cms.uint32( 4 )
    )
  )
)
process.hltESPMixedLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
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
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet( 
    useRingSlector = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    minRing = cms.int32( 1 ),
    maxRing = cms.int32( 1 )
  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPMuTrackJpsiTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
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
process.hltESPMuTrackJpsiTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
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
process.hltESPMuonCkfTrajectoryBuilder = cms.ESProducer( "MuonCkfTrajectoryBuilderESProducer",
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
process.hltESPMuonCkfTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
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
process.hltESPMuonDetLayerGeometryESProducer = cms.ESProducer( "MuonDetLayerGeometryESProducer",
  appendToDataLabel = cms.string( "" )
)
process.hltESPMuonTransientTrackingRecHitBuilder = cms.ESProducer( "MuonTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPMuonTransientTrackingRecHitBuilder" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPPixelCPEGeneric = cms.ESProducer( "PixelCPEGenericESProducer",
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
process.hltESPPixelCPETemplateReco = cms.ESProducer( "PixelCPETemplateRecoESProducer",
  ComponentName = cms.string( "hltESPPixelCPETemplateReco" ),
  DoCosmics = cms.bool( False ),
  LoadTemplatesFromDB = cms.bool( True ),
  speed = cms.int32( -2 ),
  UseClusterSplitter = cms.bool( False ),
  appendToDataLabel = cms.string( "" ),
  TanLorentzAnglePerTesla = cms.double( 0.106 ),
  PixelErrorParametrization = cms.string( "" ),
  Alpha2Order = cms.bool( True ),
  ClusterProbComputationFlag = cms.int32( 0 )
)
process.hltESPPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
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
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPPixelLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPPixelLayerTripletsHITHB = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTripletsHITHB" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3' ),
  BPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPPixelLayerTripletsHITHE = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltESPPixelLayerTripletsHITHE" ),
  layerList = cms.vstring( 'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0027 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 )
  ),
  FPix = cms.PSet( 
    useErrorsFromParam = cms.bool( True ),
    hitErrorRPhi = cms.double( 0.0051 ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltESPPromptTrackCountingESProducer = cms.ESProducer( "PromptTrackCountingESProducer",
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
process.hltESPRKTrajectoryFitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPRKFitter" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPRKTrajectorySmoother = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPRKSmoother" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPGlobalDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPRungeKuttaTrackerPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
  ComponentName = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Mass = cms.double( 0.105 ),
  MaxDPhi = cms.double( 1.6 ),
  useRungeKutta = cms.bool( True ),
  ptMin = cms.double( -1.0 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPSiStripRegionConnectivity = cms.ESProducer( "SiStripRegionConnectivity",
  EtaDivisions = cms.untracked.uint32( 20 ),
  PhiDivisions = cms.untracked.uint32( 20 ),
  EtaMax = cms.untracked.double( 2.5 )
)
process.hltESPSmartPropagator = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "hltESPSteppingHelixPropagatorAlong" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPSmartPropagatorAny = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAny" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterial" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPSmartPropagatorAnyOpposite = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorAnyOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterialOpposite" ),
  MuonPropagator = cms.string( "SteppingHelixPropagatorAny" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPSmartPropagatorOpposite = cms.ESProducer( "SmartPropagatorESProducer",
  ComponentName = cms.string( "hltESPSmartPropagatorOpposite" ),
  PropagationDirection = cms.string( "oppositeToMomentum" ),
  Epsilon = cms.double( 5.0 ),
  TrackerPropagator = cms.string( "PropagatorWithMaterialOpposite" ),
  MuonPropagator = cms.string( "hltESPSteppingHelixPropagatorOpposite" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPSoftLeptonByDistance = cms.ESProducer( "LeptonTaggerByDistanceESProducer",
  appendToDataLabel = cms.string( "" ),
  distance = cms.double( 0.5 )
)
process.hltESPSoftLeptonByPt = cms.ESProducer( "LeptonTaggerByPtESProducer",
  appendToDataLabel = cms.string( "" ),
  ipSign = cms.string( "any" )
)
process.hltESPSteppingHelixPropagatorAlong = cms.ESProducer( "SteppingHelixPropagatorESProducer",
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
process.hltESPSteppingHelixPropagatorOpposite = cms.ESProducer( "SteppingHelixPropagatorESProducer",
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
process.hltESPStraightLinePropagator = cms.ESProducer( "StraightLinePropagatorESProducer",
  ComponentName = cms.string( "hltESPStraightLinePropagator" ),
  PropagationDirection = cms.string( "alongMomentum" ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBWithTrackAngle = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBWithTrackAngle" ),
  StripCPE = cms.string( "StripCPEfromTrackAngle" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderAngleAndTemplate = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderAngleAndTemplate" ),
  StripCPE = cms.string( "StripCPEfromTrackAngle" ),
  PixelCPE = cms.string( "hltESPPixelCPETemplateReco" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderPixelOnly = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderPixelOnly" ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTTRHBuilderWithoutAngle4PixelTriplets = cms.ESProducer( "TkTransientTrackingRecHitBuilderESProducer",
  ComponentName = cms.string( "hltESPTTRHBuilderWithoutAngle4PixelTriplets" ),
  StripCPE = cms.string( "Fake" ),
  PixelCPE = cms.string( "hltESPPixelCPEGeneric" ),
  Matcher = cms.string( "StandardMatcher" ),
  ComputeCoarseLocalPositionFromDisk = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrackCounting3D1st = cms.ESProducer( "TrackCountingESProducer",
  appendToDataLabel = cms.string( "" ),
  nthTrack = cms.int32( 1 ),
  impactParameterType = cms.int32( 0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 5.0 ),
  maximumDistanceToJetAxis = cms.double( 0.07 ),
  trackQualityClass = cms.string( "any" )
)
process.hltESPTrackCounting3D2nd = cms.ESProducer( "TrackCountingESProducer",
  appendToDataLabel = cms.string( "" ),
  nthTrack = cms.int32( 2 ),
  impactParameterType = cms.int32( 0 ),
  deltaR = cms.double( -1.0 ),
  maximumDecayLength = cms.double( 5.0 ),
  maximumDistanceToJetAxis = cms.double( 0.07 ),
  trackQualityClass = cms.string( "any" )
)
process.hltESPTrackerRecoGeometryESProducer = cms.ESProducer( "TrackerRecoGeometryESProducer",
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectoryBuilderIT = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltESPTrajectoryBuilderIT" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator9" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltESPTrajectoryFilterIT" ),
  maxCand = cms.int32( 2 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectoryBuilderL3 = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
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
process.hltESPTrajectoryCleanerBySharedHits = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedHits" ),
  appendToDataLabel = cms.string( "" ),
  fractionShared = cms.double( 0.5 ),
  allowSharedFirstHit = cms.bool( False )
)
process.hltESPTrajectoryCleanerBySharedSeeds = cms.ESProducer( "TrajectoryCleanerESProducer",
  ComponentName = cms.string( "hltESPTrajectoryCleanerBySharedSeeds" ),
  ComponentType = cms.string( "TrajectoryCleanerBySharedSeeds" ),
  appendToDataLabel = cms.string( "" ),
  fractionShared = cms.double( 0.5 ),
  allowSharedFirstHit = cms.bool( True )
)
process.hltESPTrajectoryFilterIT = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFilterIT" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.3 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 100 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 3 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hltESPTrajectoryFilterL3 = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFilterL3" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.5 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 1000000000 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 5 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hltESPTrajectoryFitterRK = cms.ESProducer( "KFTrajectoryFitterESProducer",
  ComponentName = cms.string( "hltESPTrajectoryFitterRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPTrajectorySmootherRK = cms.ESProducer( "KFTrajectorySmootherESProducer",
  ComponentName = cms.string( "hltESPTrajectorySmootherRK" ),
  Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
  Updator = cms.string( "hltESPKFUpdator" ),
  Estimator = cms.string( "hltESPChi2MeasurementEstimator" ),
  RecoGeometry = cms.string( "hltESPDummyDetLayerGeometry" ),
  errorRescaling = cms.double( 100.0 ),
  minHits = cms.int32( 3 ),
  appendToDataLabel = cms.string( "" )
)
process.hltESPbJetRegionalTrajectoryBuilder = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
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
process.hltESPbJetRegionalTrajectoryFilter = cms.ESProducer( "TrajectoryFilterESProducer",
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
process.hltIter1ESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltIter1ESPMeasurementTracker" ),
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
  skipClusters = cms.InputTag( "hltIter1ClustersRefRemoval" ),
  stripClusterProducer = cms.string( "hltIter1SiStripClusters" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TID = cms.PSet( 
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
process.hltIter1ESPPixelLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltIter1ESPPixelLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg' ),
  BPix = cms.PSet( 
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter1ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0027 )
  ),
  FPix = cms.PSet( 
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter1ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0051 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltIter1ESPTrajectoryBuilderIT = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltIter1ESPTrajectoryBuilderIT" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator16" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltIter1ESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltIter1ESPTrajectoryFilterIT" ),
  maxCand = cms.int32( 2 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltIter1ESPTrajectoryFilterIT = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltIter1ESPTrajectoryFilterIT" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.2 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 100 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 3 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hltIter2ESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltIter2ESPMeasurementTracker" ),
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
  skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
  stripClusterProducer = cms.string( "hltIter2SiStripClusters" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TID = cms.PSet( 
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
process.hltIter2ESPPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltIter2ESPPixelLayerPairs" ),
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
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0027 )
  ),
  FPix = cms.PSet( 
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0051 )
  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  ),
  TOB = cms.PSet(  )
)
process.hltIter2ESPTrajectoryBuilderIT = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltIter2ESPTrajectoryBuilderIT" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator16" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltIter2ESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltIter2ESPTrajectoryFilterIT" ),
  maxCand = cms.int32( 2 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltIter2ESPTrajectoryFilterIT = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltIter2ESPTrajectoryFilterIT" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.3 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 1 ),
    maxNumberOfHits = cms.int32( 100 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 3 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hltIter3ESPLayerTriplets = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltIter3ESPLayerTriplets" ),
  layerList = cms.vstring( 'BPix1+BPix2+BPix3',
    'BPix1+BPix2+FPix1_pos',
    'BPix1+BPix2+FPix1_neg',
    'BPix1+FPix1_pos+FPix2_pos',
    'BPix1+FPix1_neg+FPix2_neg',
    'BPix2+FPix1_pos+FPix2_pos',
    'BPix2+FPix1_neg+FPix2_neg',
    'FPix1_pos+FPix2_pos+TEC1_pos',
    'FPix1_neg+FPix2_neg+TEC1_neg',
    'FPix2_pos+TEC2_pos+TEC3_pos',
    'FPix2_neg+TEC2_neg+TEC3_neg',
    'BPix2+BPix3+TIB1',
    'BPix2+BPix3+TIB2',
    'BPix1+BPix3+TIB1',
    'BPix1+BPix3+TIB2',
    'BPix1+BPix2+TIB1',
    'BPix1+BPix2+TIB2' ),
  BPix = cms.PSet( 
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0060 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter3ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0027 )
  ),
  FPix = cms.PSet( 
    HitProducer = cms.string( "hltSiPixelRecHits" ),
    hitErrorRZ = cms.double( 0.0036 ),
    useErrorsFromParam = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
    skipClusters = cms.InputTag( "hltIter3ClustersRefRemoval" ),
    hitErrorRPhi = cms.double( 0.0051 )
  ),
  TEC = cms.PSet( 
    useRingSlector = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    minRing = cms.int32( 1 ),
    maxRing = cms.int32( 1 )
  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ) ),
  TOB = cms.PSet(  )
)
process.hltIter3ESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltIter3ESPMeasurementTracker" ),
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
  skipClusters = cms.InputTag( "hltIter3ClustersRefRemoval" ),
  stripClusterProducer = cms.string( "hltIter3SiStripClusters" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TID = cms.PSet( 
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
process.hltIter3ESPTrajectoryBuilderIT = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltIter3ESPTrajectoryBuilderIT" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator16" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltIter3ESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltIter3ESPTrajectoryFilterIT" ),
  maxCand = cms.int32( 1 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" )
)
process.hltIter3ESPTrajectoryFilterIT = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltIter3ESPTrajectoryFilterIT" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.3 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 0 ),
    maxNumberOfHits = cms.int32( 100 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 3 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hltIter4ESPMeasurementTracker = cms.ESProducer( "MeasurementTrackerESProducer",
  ComponentName = cms.string( "hltIter4ESPMeasurementTracker" ),
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
  skipClusters = cms.InputTag( "hltIter4ClustersRefRemoval" ),
  stripClusterProducer = cms.string( "hltIter4SiStripClusters" ),
  stripLazyGetterProducer = cms.string( "hltSiStripRawToClustersFacility" ),
  appendToDataLabel = cms.string( "" ),
  inactivePixelDetectorLabels = cms.VInputTag(  ),
  inactiveStripDetectorLabels = cms.VInputTag( 'hltSiStripExcludedFEDListProducer' ),
  badStripCuts = cms.PSet( 
    TOB = cms.PSet( 
      maxConsecutiveBad = cms.uint32( 9999 ),
      maxBad = cms.uint32( 9999 )
    ),
    TID = cms.PSet( 
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
process.hltIter4ESPPixelLayerPairs = cms.ESProducer( "SeedingLayersESProducer",
  appendToDataLabel = cms.string( "" ),
  ComponentName = cms.string( "hltIter4ESPPixelLayerPairs" ),
  layerList = cms.vstring( 'TIB1+TIB2' ),
  BPix = cms.PSet(  ),
  FPix = cms.PSet(  ),
  TEC = cms.PSet(  ),
  TID = cms.PSet(  ),
  TIB = cms.PSet(  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ) ),
  TOB = cms.PSet(  )
)
process.hltIter4ESPTrajectoryBuilderIT = cms.ESProducer( "CkfTrajectoryBuilderESProducer",
  ComponentName = cms.string( "hltIter4ESPTrajectoryBuilderIT" ),
  updator = cms.string( "hltESPKFUpdator" ),
  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
  estimator = cms.string( "hltESPChi2MeasurementEstimator16" ),
  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  MeasurementTrackerName = cms.string( "hltIter4ESPMeasurementTracker" ),
  trajectoryFilterName = cms.string( "hltIter4ESPTrajectoryFilterIT" ),
  maxCand = cms.int32( 1 ),
  lostHitPenalty = cms.double( 30.0 ),
  intermediateCleaning = cms.bool( True ),
  alwaysUseInvalidHits = cms.bool( False ),
  appendToDataLabel = cms.string( "" ),
  minNrOfHitsForRebuild = cms.untracked.int32( 4 )
)
process.hltIter4ESPTrajectoryFilterIT = cms.ESProducer( "TrajectoryFilterESProducer",
  ComponentName = cms.string( "hltIter4ESPTrajectoryFilterIT" ),
  appendToDataLabel = cms.string( "" ),
  filterPset = cms.PSet( 
    minPt = cms.double( 0.3 ),
    minHitsMinPt = cms.int32( 3 ),
    ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
    maxLostHits = cms.int32( 0 ),
    maxNumberOfHits = cms.int32( 100 ),
    maxConsecLostHits = cms.int32( 1 ),
    minimumNumberOfHits = cms.int32( 6 ),
    nSigmaMinPt = cms.double( 5.0 ),
    chargeSignificance = cms.double( -1.0 )
  )
)
process.hoDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "HODetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.087 ),
  nEta = cms.int32( 30 ),
  nPhi = cms.int32( 72 ),
  includeBadChambers = cms.bool( False )
)
process.muonDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "MuonDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.125 ),
  nEta = cms.int32( 48 ),
  nPhi = cms.int32( 48 ),
  includeBadChambers = cms.bool( False )
)
process.navigationSchoolESProducer = cms.ESProducer( "NavigationSchoolESProducer",
  ComponentName = cms.string( "SimpleNavigationSchool" ),
  appendToDataLabel = cms.string( "" )
)
process.preshowerDetIdAssociator = cms.ESProducer( "DetIdAssociatorESProducer",
  ComponentName = cms.string( "PreshowerDetIdAssociator" ),
  appendToDataLabel = cms.string( "" ),
  etaBinSize = cms.double( 0.1 ),
  nEta = cms.int32( 60 ),
  nPhi = cms.int32( 30 ),
  includeBadChambers = cms.bool( False )
)
process.siPixelTemplateDBObjectESProducer = cms.ESProducer( "SiPixelTemplateDBObjectESProducer",
  appendToDataLabel = cms.string( "" )
)
process.sistripconn = cms.ESProducer( "SiStripConnectivity" )

process.DQM = cms.Service( "DQM",
)
process.DQMStore = cms.Service( "DQMStore",
)
process.DTDataIntegrityTask = cms.Service( "DTDataIntegrityTask",
    getSCInfo = cms.untracked.bool( True ),
    fedIntegrityFolder = cms.untracked.string( "DT/FEDIntegrity_EvF" ),
    processingMode = cms.untracked.string( "HLT" )
)
process.MessageLogger = cms.Service( "MessageLogger",
    destinations = cms.untracked.vstring( 'warnings',
      'errors',
      'infos',
      'debugs',
      'cout',
      'cerr' ),
    categories = cms.untracked.vstring( 'FwkJob',
      'FwkReport',
      'FwkSummary',
      'Root_NoDictionary' ),
    statistics = cms.untracked.vstring( 'cerr' ),
    cerr = cms.untracked.PSet( 
      INFO = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      noTimeStamps = cms.untracked.bool( False ),
      FwkReport = cms.untracked.PSet( 
        reportEvery = cms.untracked.int32( 1 ),
        limit = cms.untracked.int32( 0 )
      ),
      default = cms.untracked.PSet(  limit = cms.untracked.int32( 10000000 ) ),
      Root_NoDictionary = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      FwkJob = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      FwkSummary = cms.untracked.PSet( 
        reportEvery = cms.untracked.int32( 1 ),
        limit = cms.untracked.int32( 10000000 )
      ),
      threshold = cms.untracked.string( "INFO" ),
    ),
    cout = cms.untracked.PSet( 
      threshold = cms.untracked.string( "ERROR" ),
    ),
    errors = cms.untracked.PSet( 
      threshold = cms.untracked.string( "INFO" ),
      placeholder = cms.untracked.bool( True ),
    ),
    warnings = cms.untracked.PSet( 
      threshold = cms.untracked.string( "INFO" ),
      placeholder = cms.untracked.bool( True ),
    ),
    infos = cms.untracked.PSet( 
      threshold = cms.untracked.string( "INFO" ),
      Root_NoDictionary = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      placeholder = cms.untracked.bool( True ),
    ),
    debugs = cms.untracked.PSet( 
      threshold = cms.untracked.string( "INFO" ),
      placeholder = cms.untracked.bool( True ),
    ),
    fwkJobReports = cms.untracked.vstring( 'FrameworkJobReport' ),
    FrameworkJobReport = cms.untracked.PSet( 
      default = cms.untracked.PSet(  limit = cms.untracked.int32( 0 ) ),
      FwkJob = cms.untracked.PSet(  limit = cms.untracked.int32( 10000000 ) )
    ),
    suppressWarning = cms.untracked.vstring( 'hltL3MuonsOIState',
      'hltPixelVertices3DbbPhi',
      'hltSiPixelDigis',
      'hltPixelTracksForHighMult',
      'hltSiPixelClusters',
      'hltLightPFTracks',
      'hltPixelTracks',
      'hltOnlineBeamSpot',
      'hltL3MuonsOIHit',
      'hltHITPixelTracksHE',
      'hltHITPixelTracksHB',
      'hltL3MuonsIOHit' ),
    threshold = cms.untracked.string( "INFO" ),
)
process.MicroStateService = cms.Service( "MicroStateService",
)
process.ModuleWebRegistry = cms.Service( "ModuleWebRegistry",
)
process.PrescaleService = cms.Service( "PrescaleService",
    lvl1DefaultLabel = cms.untracked.string( "3e33" ),
    lvl1Labels = cms.vstring( '5e33',
      '4e33',
      '3e33',
      '2.5e33',
      'Cosmics',
      'Cosmics + High Random' ),
    prescaleTable = cms.VPSet( 
      cms.PSet(  pathName = cms.string( "HLT_HIZeroBias_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIZeroBiasXOR_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIZeroBiasPixel_SingleTrack_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasBSC_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasBSC_OR_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasHF_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasHf_OR_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasHfOrBSC_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasPixel_SingleTrack_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasZDC_Calo_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasZDC_Scint_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIMinBiasZDCPixel_SingleTrack_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIBptxXOR_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL1Algo_BptxXOR_BSC_OR_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL1SingleMu3_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL1SingleMu5_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL1SingleMu7_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL1DoubleMuOpen_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL2Mu3_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL2Mu5Tight_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL2Mu20_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL2DoubleMu0_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIL2DoubleMu3_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIUpcEcal_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIUpcMu_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIStoppedHSCP35_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIActivityHF_Coincidence3_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIActivityHF_Single3_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIClusterVertexCompatibility_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HICentralityVeto_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLT_HIRandom_v1" ),
        prescales = cms.vuint32( 0, 0, 0, 0, 0, 0 )
      ),
      cms.PSet(  pathName = cms.string( "HLTMONOutput" ),
        prescales = cms.vuint32( 100, 100, 100, 100, 100, 100 )
      )
    )
)
process.UpdaterService = cms.Service( "UpdaterService",
)

process.hltGetRaw = cms.EDAnalyzer( "HLTGetRaw",
    RawDataCollection = cms.InputTag( "source" )
)
process.hltBoolFalse = cms.EDFilter( "HLTBool",
    result = cms.bool( False )
)
process.hltCalibrationEventsFilter = cms.EDFilter( "HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32( 2 )
)
process.hltGtDigis = cms.EDProducer( "L1GlobalTriggerRawToDigi",
    DaqGtInputTag = cms.InputTag( "source" ),
    DaqGtFedId = cms.untracked.int32( 813 ),
    ActiveBoardsMask = cms.uint32( 0xffff ),
    UnpackBxInEvent = cms.int32( 5 ),
    Verbosity = cms.untracked.int32( 0 )
)
process.hltPreDTCalibration = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltDTCalibrationRaw = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "source" ),
    fedList = cms.vuint32( 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780 )
)
process.hltBoolEnd = cms.EDFilter( "HLTBool",
    result = cms.bool( True )
)
process.hltPreEcalCalibration = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltEcalCalibrationRaw = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "source" ),
    fedList = cms.vuint32( 601, 602, 603, 604, 605, 606, 607, 608, 609, 610, 611, 612, 613, 614, 615, 616, 617, 618, 619, 620, 621, 622, 623, 624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654 )
)
process.hltPreHcalCalibration = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHcalCalibTypeFilter = cms.EDFilter( "HLTHcalCalibTypeFilter",
    InputTag = cms.InputTag( "source" ),
    CalibTypes = cms.vint32( 1, 2, 3, 4, 5, 6 )
)
process.hltHcalCalibrationRaw = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "source" ),
    fedList = cms.vuint32( 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 726, 727, 728, 729, 730, 731 )
)
process.hltPreTrackerCalibration = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltLaserAlignmentEventFilter = cms.EDFilter( "LaserAlignmentEventFilter",
    FedInputTag = cms.InputTag( "source" ),
    SIGNAL_Filter = cms.bool( True ),
    SINGLE_CHANNEL_THRESH = cms.uint32( 11 ),
    CHANNEL_COUNT_THRESH = cms.uint32( 8 ),
    FED_IDs = cms.vint32( 260, 261, 262, 263, 264, 265, 266, 267, 269, 270, 273, 274, 277, 278, 281, 282, 284, 285, 288, 289, 292, 293, 294, 295, 300, 301, 304, 305, 308, 309, 310, 311, 316, 317, 324, 325, 329, 330, 331, 332, 339, 340, 341, 342, 349, 350, 351, 352, 164, 165, 172, 173, 177, 178, 179, 180, 187, 188, 189, 190, 197, 198, 199, 200, 204, 205, 208, 209, 212, 213, 214, 215, 220, 221, 224, 225, 228, 229, 230, 231, 236, 237, 238, 239, 240, 241, 242, 243, 245, 246, 249, 250, 253, 254, 257, 258, 478, 476, 477, 482, 484, 480, 481, 474, 459, 460, 461, 463, 485, 487, 488, 489, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 288, 289, 292, 293, 300, 301, 304, 305, 310, 311, 316, 317, 329, 330, 339, 340, 341, 342, 349, 350, 164, 165, 177, 178, 179, 180, 189, 190, 197, 198, 204, 205, 212, 213, 220, 221, 224, 225, 230, 231 ),
    SIGNAL_IDs = cms.vint32( 470389128, 470389384, 470389640, 470389896, 470390152, 470390408, 470390664, 470390920, 470389192, 470389448, 470389704, 470389960, 470390216, 470390472, 470390728, 470390984, 470126984, 470127240, 470127496, 470127752, 470128008, 470128264, 470128520, 470128776, 470127048, 470127304, 470127560, 470127816, 470128072, 470128328, 470128584, 470128840, 436232506, 436232826, 436233146, 436233466, 369174604, 369174812, 369175068, 369175292, 470307468, 470307716, 470308236, 470308748, 470308996, 470045316, 470045580, 470046084, 470046596, 470046860 )
)
process.hltTrackerCalibrationRaw = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "source" ),
    fedList = cms.vuint32( 260, 261, 262, 263, 264, 265, 266, 267, 269, 270, 273, 274, 277, 278, 281, 282, 284, 285, 288, 289, 292, 293, 294, 295, 300, 301, 304, 305, 308, 309, 310, 311, 316, 317, 324, 325, 329, 330, 331, 332, 339, 340, 341, 342, 349, 350, 351, 352, 164, 165, 172, 173, 177, 178, 179, 180, 187, 188, 189, 190, 197, 198, 199, 200, 204, 205, 208, 209, 212, 213, 214, 215, 220, 221, 224, 225, 228, 229, 230, 231, 236, 237, 238, 239, 240, 241, 242, 243, 245, 246, 249, 250, 253, 254, 257, 258, 478, 476, 477, 482, 484, 480, 481, 474, 459, 460, 461, 463, 485, 487, 488, 489, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 288, 289, 292, 293, 300, 301, 304, 305, 310, 311, 316, 317, 329, 330, 339, 340, 341, 342, 349, 350, 164, 165, 177, 178, 179, 180, 189, 190, 197, 198, 204, 205, 212, 213, 220, 221, 224, 225, 230, 231 )
)
process.hltPreLogMonitor = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltLogMonitorFilter = cms.EDFilter( "HLTLogMonitorFilter",
    default_threshold = cms.uint32( 10 ),
    categories = cms.VPSet( 
    )
)
process.hltTriggerType = cms.EDFilter( "HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32( 1 )
)
process.hltGctDigis = cms.EDProducer( "GctRawToDigi",
    inputLabel = cms.InputTag( "source" ),
    gctFedId = cms.untracked.int32( 745 ),
    hltMode = cms.bool( True ),
    numberOfGctSamplesToUnpack = cms.uint32( 1 ),
    numberOfRctSamplesToUnpack = cms.uint32( 1 ),
    unpackSharedRegions = cms.bool( False ),
    unpackerVersion = cms.uint32( 0 )
)
process.hltL1GtObjectMap = cms.EDProducer( "L1GlobalTrigger",
    GmtInputTag = cms.InputTag( "hltGtDigis" ),
    GctInputTag = cms.InputTag( "hltGctDigis" ),
    CastorInputTag = cms.InputTag( "castorL1Digis" ),
    ProduceL1GtDaqRecord = cms.bool( False ),
    ProduceL1GtEvmRecord = cms.bool( False ),
    ProduceL1GtObjectMapRecord = cms.bool( True ),
    WritePsbL1GtDaqRecord = cms.bool( False ),
    ReadTechnicalTriggerRecords = cms.bool( True ),
    EmulateBxInEvent = cms.int32( 1 ),
    AlternativeNrBxBoardDaq = cms.uint32( 0 ),
    AlternativeNrBxBoardEvm = cms.uint32( 0 ),
    BstLengthBytes = cms.int32( -1 ),
    TechnicalTriggersInputTags = cms.VInputTag( 'simBscDigis' ),
    RecordLength = cms.vint32( 3, 0 )
)
process.hltL1extraParticles = cms.EDProducer( "L1ExtraParticlesProd",
    produceMuonParticles = cms.bool( True ),
    muonSource = cms.InputTag( "hltGtDigis" ),
    produceCaloParticles = cms.bool( True ),
    isolatedEmSource = cms.InputTag( 'hltGctDigis','isoEm' ),
    nonIsolatedEmSource = cms.InputTag( 'hltGctDigis','nonIsoEm' ),
    centralJetSource = cms.InputTag( 'hltGctDigis','cenJets' ),
    forwardJetSource = cms.InputTag( 'hltGctDigis','forJets' ),
    tauJetSource = cms.InputTag( 'hltGctDigis','tauJets' ),
    etTotalSource = cms.InputTag( "hltGctDigis" ),
    etHadSource = cms.InputTag( "hltGctDigis" ),
    etMissSource = cms.InputTag( "hltGctDigis" ),
    htMissSource = cms.InputTag( "hltGctDigis" ),
    hfRingEtSumsSource = cms.InputTag( "hltGctDigis" ),
    hfRingBitCountsSource = cms.InputTag( "hltGctDigis" ),
    centralBxOnly = cms.bool( True ),
    ignoreHtMiss = cms.bool( False )
)
process.hltScalersRawToDigi = cms.EDProducer( "ScalersRawToDigi",
    scalersInputTag = cms.InputTag( "source" )
)
process.hltOnlineBeamSpot = cms.EDProducer( "BeamSpotOnlineProducer",
    label = cms.InputTag( "hltScalersRawToDigi" ),
    changeToCMSCoordinates = cms.bool( False ),
    maxRadius = cms.double( 2.0 ),
    maxZ = cms.double( 40.0 ),
    setSigmaZ = cms.double( 0.0 ),
    gtEvmLabel = cms.InputTag( "" )
)
process.hltOfflineBeamSpot = cms.EDProducer( "BeamSpotProducer" )
process.hltL1sHIZeroBias = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BptxPlusANDMinus" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIZeroBias = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1BptxXOR = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BptxXOR" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIZeroBiasXOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIZeroBiasXOR = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BptxPlusANDMinus OR L1_BptxXOR" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIZeroBiasPixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltSiPixelDigis = cms.EDProducer( "SiPixelRawToDigi",
    IncludeErrors = cms.bool( False ),
    UseQualityInfo = cms.bool( False ),
    InputLabel = cms.InputTag( "source" )
)
process.hltHISiPixelClusters = cms.EDProducer( "SiPixelClusterProducer",
    src = cms.InputTag( "hltSiPixelDigis" ),
    maxNumberOfClusters = cms.int32( -1 ),
    payloadType = cms.string( "HLT" ),
    ChannelThreshold = cms.int32( 1000 ),
    SeedThreshold = cms.int32( 1000 ),
    ClusterThreshold = cms.double( 4000.0 ),
    VCaltoElectronGain = cms.int32( 65 ),
    VCaltoElectronOffset = cms.int32( -414 ),
    MissCalibrate = cms.untracked.bool( True ),
    SplitClusters = cms.bool( False )
)
process.hltHISiPixelRecHits = cms.EDProducer( "SiPixelRecHitConverter",
    src = cms.InputTag( "hltHISiPixelClusters" ),
    CPE = cms.string( "hltESPPixelCPEGeneric" )
)
process.hltHIPixelClusterVertices = cms.EDProducer( "HIPixelClusterVtxProducer",
    pixelRecHits = cms.string( "hltHISiPixelRecHits" ),
    minZ = cms.double( -20.0 ),
    maxZ = cms.double( 20.05 ),
    zStep = cms.double( 0.1 )
)
process.hltPixelTracksForHITrackTrigger = cms.EDProducer( "PixelTrackProducer",
    useFilterWithES = cms.bool( False ),
    RegionFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "HITrackingRegionForPrimaryVtxProducer" ),
      RegionPSet = cms.PSet( 
        precise = cms.bool( True ),
        doVariablePtMin = cms.bool( True ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        useFixedError = cms.bool( True ),
        directionYCoord = cms.double( 1.0 ),
        sigmaZVertex = cms.double( 3.0 ),
        fixedError = cms.double( 2.0 ),
        directionXCoord = cms.double( 1.0 ),
        directionZCoord = cms.double( 0.0 ),
        VertexCollection = cms.InputTag( "hltHIPixelClusterVertices" ),
        ptMin = cms.double( 0.7 ),
        siPixelRecHits = cms.InputTag( "hltHISiPixelRecHits" ),
        nSigmaZ = cms.double( 3.0 ),
        useFoundVertices = cms.bool( True ),
        originRadius = cms.double( 0.1 )
      )
    ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitTripletGenerator" ),
      GeneratorPSet = cms.PSet( 
        useBending = cms.bool( True ),
        useFixedPreFiltering = cms.bool( False ),
        maxElement = cms.uint32( 100000 ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRPhitolerance = cms.double( 0.06 ),
        useMultScattering = cms.bool( True ),
        ComponentName = cms.string( "PixelTripletHLTGenerator" ),
        extraHitRZtolerance = cms.double( 0.06 ),
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
      ),
      SeedingLayers = cms.string( "hltESPHIPixelLayerTriplets" )
    ),
    FitterPSet = cms.PSet( 
      ComponentName = cms.string( "PixelFitterByHelixProjections" ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" )
    ),
    FilterPSet = cms.PSet( 
      doVariablePtMin = cms.bool( True ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      chi2 = cms.double( 1000.0 ),
      ComponentName = cms.string( "HIProtoTrackFilter" ),
      ptMin = cms.double( 1.0 ),
      siPixelRecHits = cms.InputTag( "hltHISiPixelRecHits" ),
      tipMax = cms.double( 1.0 )
    ),
    CleanerPSet = cms.PSet(  ComponentName = cms.string( "PixelTrackCleanerBySharedHits" ) )
)
process.hltPixelCandsForHITrackTrigger = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltPixelTracksForHITrackTrigger" ),
    particleType = cms.string( "pi+" )
)
process.hltHISinglePixelTrackFilter = cms.EDFilter( "HLTPixlMBFilt",
    pixlTag = cms.InputTag( "hltPixelCandsForHITrackTrigger" ),
    MinPt = cms.double( 0.0 ),
    MinTrks = cms.uint32( 1 ),
    MinSep = cms.double( 1.0 )
)
process.hltL1sHIMinBiasBSC = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BscMinBiasThreshold1_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasBSC = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasBSCOR = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BscMinBiasOR_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasBSCOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasHF = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HcalHfCoincidencePm_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasHF = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasHfOr = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HcalHfMmOrPpOrPm_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasHfOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasHfOrBSC = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_HcalHfCoincPmORBscMinBiasThresh1_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasHfOrBSC = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreHIMinBiasPixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasZDC = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasZDCCalo = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasZDCCaloPlusOrMinus = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasZDCCaloPlusOrMinus = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasZDCScint = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ZdcScintTightVertex_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasZDCScint = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIMinBiasZDCPixelSingleTrack = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ZdcCaloPlus_ZdcCaloMinus_BptxAND OR L1_ZdcScintTightVertex_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIMinBiasZDCPixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreHIBptxXOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1BptxXORBscMinBiasOR = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_BptxXOR_BscMinBiasOR" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL1AlgoBptxXORBSCOR = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1SingleMu3BptxAND = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu3_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL1SingleMu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1SingleMu5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu5" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL1SingleMu5 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1SingleMu7 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu7" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL1SingleMu7 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1DoubleMuOpenBptxAND = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleMuOpen_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL1DoubleMuOpen = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoubleMuLevel1PathL1OpenFiltered = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "hltL1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMuOpenBptxAND" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
process.hltL1sL1SingleMu3BptxANDwithBSCHF = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleMu3_BptxAND AND ( L1_BscMinBiasOR_BptxPlusORMinus OR L1_HcalHfMmOrPpOrPm OR L1_BscMinBiasThreshold1_BptxAND)" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL2Mu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIL1SingleMu3withBSCHFL1Filtered = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "hltL1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1SingleMu3BptxANDwithBSCHF" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 1 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( False ),
    SelectQualities = cms.vint32(  )
)
process.hltMuonDTDigis = cms.EDProducer( "DTUnpackingModule",
    dataType = cms.string( "DDU" ),
    fedbyType = cms.bool( False ),
    inputLabel = cms.InputTag( "source" ),
    useStandardFEDid = cms.bool( True ),
    dqmOnly = cms.bool( False ),
    rosParameters = cms.PSet(  ),
    readOutParameters = cms.PSet( 
      debug = cms.untracked.bool( False ),
      rosParameters = cms.PSet( 
        writeSC = cms.untracked.bool( True ),
        readingDDU = cms.untracked.bool( True ),
        performDataIntegrityMonitor = cms.untracked.bool( False ),
        readDDUIDfromDDU = cms.untracked.bool( True ),
        debug = cms.untracked.bool( False ),
        localDAQ = cms.untracked.bool( False )
      ),
      localDAQ = cms.untracked.bool( False ),
      performDataIntegrityMonitor = cms.untracked.bool( False )
    )
)
process.hltDt1DRecHits = cms.EDProducer( "DTRecHitProducer",
    dtDigiLabel = cms.InputTag( "hltMuonDTDigis" ),
    recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
    recAlgoConfig = cms.PSet( 
      tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
      minTime = cms.double( -3.0 ),
      stepTwoFromDigi = cms.bool( False ),
      doVdriftCorr = cms.bool( False ),
      debug = cms.untracked.bool( False ),
      maxTime = cms.double( 420.0 ),
      tTrigModeConfig = cms.PSet( 
        vPropWire = cms.double( 24.4 ),
        doTOFCorrection = cms.bool( True ),
        tofCorrType = cms.int32( 0 ),
        wirePropCorrType = cms.int32( 0 ),
        tTrigLabel = cms.string( "" ),
        doWirePropCorrection = cms.bool( True ),
        doT0Correction = cms.bool( True ),
        debug = cms.untracked.bool( False )
      )
    )
)
process.hltDt4DSegments = cms.EDProducer( "DTRecSegment4DProducer",
    debug = cms.untracked.bool( False ),
    recHits1DLabel = cms.InputTag( "hltDt1DRecHits" ),
    recHits2DLabel = cms.InputTag( "dt2DSegments" ),
    Reco4DAlgoName = cms.string( "DTCombinatorialPatternReco4D" ),
    Reco4DAlgoConfig = cms.PSet( 
      segmCleanerMode = cms.int32( 2 ),
      Reco2DAlgoName = cms.string( "DTCombinatorialPatternReco" ),
      recAlgoConfig = cms.PSet( 
        tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
        minTime = cms.double( -3.0 ),
        stepTwoFromDigi = cms.bool( False ),
        doVdriftCorr = cms.bool( False ),
        debug = cms.untracked.bool( False ),
        maxTime = cms.double( 420.0 ),
        tTrigModeConfig = cms.PSet( 
          vPropWire = cms.double( 24.4 ),
          doTOFCorrection = cms.bool( True ),
          tofCorrType = cms.int32( 0 ),
          wirePropCorrType = cms.int32( 0 ),
          tTrigLabel = cms.string( "" ),
          doWirePropCorrection = cms.bool( True ),
          doT0Correction = cms.bool( True ),
          debug = cms.untracked.bool( False )
        )
      ),
      nSharedHitsMax = cms.int32( 2 ),
      hit_afterT0_resolution = cms.double( 0.03 ),
      Reco2DAlgoConfig = cms.PSet( 
        segmCleanerMode = cms.int32( 2 ),
        recAlgoConfig = cms.PSet( 
          tTrigMode = cms.string( "DTTTrigSyncFromDB" ),
          minTime = cms.double( -3.0 ),
          stepTwoFromDigi = cms.bool( False ),
          doVdriftCorr = cms.bool( False ),
          debug = cms.untracked.bool( False ),
          maxTime = cms.double( 420.0 ),
          tTrigModeConfig = cms.PSet( 
            vPropWire = cms.double( 24.4 ),
            doTOFCorrection = cms.bool( True ),
            tofCorrType = cms.int32( 0 ),
            wirePropCorrType = cms.int32( 0 ),
            tTrigLabel = cms.string( "" ),
            doWirePropCorrection = cms.bool( True ),
            doT0Correction = cms.bool( True ),
            debug = cms.untracked.bool( False )
          )
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
        performT0SegCorrection = cms.bool( False ),
        perform_delta_rejecting = cms.bool( False )
      ),
      performT0_vdriftSegCorrection = cms.bool( False ),
      debug = cms.untracked.bool( False ),
      recAlgo = cms.string( "DTLinearDriftFromDBAlgo" ),
      nUnSharedHitsMin = cms.int32( 2 ),
      AllDTRecHits = cms.bool( True ),
      performT0SegCorrection = cms.bool( False ),
      perform_delta_rejecting = cms.bool( False )
    )
)
process.hltMuonCSCDigis = cms.EDProducer( "CSCDCCUnpacker",
    InputObjects = cms.InputTag( "source" ),
    UseExaminer = cms.bool( True ),
    ExaminerMask = cms.uint32( 0x1febf3f6 ),
    UseSelectiveUnpacking = cms.bool( True ),
    ErrorMask = cms.uint32( 0x0 ),
    UnpackStatusDigis = cms.bool( False ),
    UseFormatStatus = cms.bool( True ),
    PrintEventNumber = cms.untracked.bool( False )
)
process.hltCsc2DRecHits = cms.EDProducer( "CSCRecHitDProducer",
    CSCUseCalibrations = cms.bool( True ),
    CSCUseStaticPedestals = cms.bool( False ),
    CSCUseTimingCorrections = cms.bool( True ),
    stripDigiTag = cms.InputTag( 'hltMuonCSCDigis','MuonCSCStripDigi' ),
    wireDigiTag = cms.InputTag( 'hltMuonCSCDigis','MuonCSCWireDigi' ),
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
process.hltCscSegments = cms.EDProducer( "CSCSegmentProducer",
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
process.hltMuonRPCDigis = cms.EDProducer( "RPCUnpackingModule",
    InputLabel = cms.InputTag( "source" ),
    doSynchro = cms.bool( False )
)
process.hltRpcRecHits = cms.EDProducer( "RPCRecHitProducer",
    rpcDigiLabel = cms.InputTag( "hltMuonRPCDigis" ),
    recAlgo = cms.string( "RPCRecHitStandardAlgo" ),
    maskSource = cms.string( "File" ),
    maskvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCMaskVec.dat" ),
    deadSource = cms.string( "File" ),
    deadvecfile = cms.FileInPath( "RecoLocalMuon/RPCRecHit/data/RPCDeadVec.dat" ),
    recAlgoConfig = cms.PSet(  )
)
process.hltL2MuonSeeds = cms.EDProducer( "L2MuonSeedGenerator",
    InputObjects = cms.InputTag( "hltL1extraParticles" ),
    GMTReadoutCollection = cms.InputTag( "hltGtDigis" ),
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
process.hltL2Muons = cms.EDProducer( "L2MuonProducer",
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
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      MuonUpdatorAtVertexParameters = cms.PSet( 
        MaxChi2 = cms.double( 1000000.0 ),
        BeamSpotPosition = cms.vdouble( 0.0, 0.0, 0.0 ),
        Propagator = cms.string( "hltESPFastSteppingHelixPropagatorOpposite" ),
        BeamSpotPositionErrors = cms.vdouble( 0.1, 0.1, 5.3 )
      ),
      VertexConstraint = cms.bool( True )
    )
)
process.hltL2MuonCandidates = cms.EDProducer( "L2MuonCandidateProducer",
    InputObjects = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' )
)
process.hltHIL2Mu3L2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHIL1SingleMu3withBSCHFL1Filtered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 3.0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True ),
    AbsEtaBins = cms.vdouble( 5.0 ),
    MinNstations = cms.vint32( 0 ),
    MinNhits = cms.vint32( 0 )
)
process.hltPreHIL2Mu5Tight = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIL2Mu5TightL2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHIL1SingleMu3withBSCHFL1Filtered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 3.0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 5.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True ),
    AbsEtaBins = cms.vdouble( 5.0 ),
    MinNstations = cms.vint32( 0 ),
    MinNhits = cms.vint32( 0 )
)
process.hltPreHIL2Mu20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIL2Mu20L2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHIL1SingleMu3withBSCHFL1Filtered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 1 ),
    MaxEta = cms.double( 3.0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 20.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True ),
    AbsEtaBins = cms.vdouble( 5.0 ),
    MinNstations = cms.vint32( 0 ),
    MinNhits = cms.vint32( 0 )
)
process.hltL1sL1DoubleMuOpenBptxANDwithBSCHF = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_DoubleMuOpen_BptxAND AND ( L1_BscMinBiasOR_BptxPlusORMinus OR L1_HcalHfMmOrPpOrPm OR L1_BscMinBiasThreshold1_BptxAND)" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIL2DoubleMu0 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoubleMuLevel1PathL1OpenWithBSCHFFiltered = cms.EDFilter( "HLTMuonL1Filter",
    CandTag = cms.InputTag( "hltL1extraParticles" ),
    PreviousCandTag = cms.InputTag( "hltL1sL1DoubleMuOpenBptxANDwithBSCHF" ),
    MaxEta = cms.double( 2.5 ),
    MinPt = cms.double( 0.0 ),
    MinN = cms.int32( 2 ),
    ExcludeSingleSegmentCSC = cms.bool( False ),
    CSCTFtag = cms.InputTag( "unused" ),
    saveTags = cms.bool( True ),
    SelectQualities = cms.vint32(  )
)
process.hltHIL2DoubleMu0L2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHIDoubleMuLevel1PathL1OpenWithBSCHFFiltered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 3.0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 0.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True ),
    AbsEtaBins = cms.vdouble( 5.0 ),
    MinNstations = cms.vint32( 0 ),
    MinNhits = cms.vint32( 0 )
)
process.hltPreHIL2DoubleMu3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIL2DoubleMu3L2Filtered = cms.EDFilter( "HLTMuonL2PreFilter",
    BeamSpotTag = cms.InputTag( "hltOnlineBeamSpot" ),
    CandTag = cms.InputTag( "hltL2MuonCandidates" ),
    PreviousCandTag = cms.InputTag( "hltHIDoubleMuLevel1PathL1OpenWithBSCHFFiltered" ),
    SeedMapTag = cms.InputTag( "hltL2Muons" ),
    MinN = cms.int32( 2 ),
    MaxEta = cms.double( 3.0 ),
    MaxDr = cms.double( 9999.0 ),
    MaxDz = cms.double( 9999.0 ),
    MinPt = cms.double( 3.0 ),
    NSigmaPt = cms.double( 0.0 ),
    saveTags = cms.bool( True ),
    AbsEtaBins = cms.vdouble( 5.0 ),
    MinNstations = cms.vint32( 0 ),
    MinNhits = cms.vint32( 0 )
)
process.hltL1sHIUpcEcal = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "NOT L1_ETT60 AND L1_DoubleEG5_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUpcEcal = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUpcMu = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "NOT L1_ETT60 AND L1_SingleMu3_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUpcMu = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sL1SingleEG5BptxAND = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG5_BptxAND" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIPhoton15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltEcalRawToRecHitFacility = cms.EDProducer( "EcalRawToRecHitFacility",
    sourceTag = cms.InputTag( "source" ),
    workerName = cms.string( "" )
)
process.hltEcalRegionalRestFEDs = cms.EDProducer( "EcalRawToRecHitRoI",
    sourceTag = cms.InputTag( "hltEcalRawToRecHitFacility" ),
    type = cms.string( "all" ),
    doES = cms.bool( False ),
    sourceTag_es = cms.InputTag( "NotNeededoESfalse" ),
    MuJobPSet = cms.PSet(  ),
    JetJobPSet = cms.VPSet( 
    ),
    EmJobPSet = cms.VPSet( 
    ),
    CandJobPSet = cms.VPSet( 
    )
)
process.hltEcalRecHitAll = cms.EDProducer( "EcalRawToRecHitProducer",
    lazyGetterTag = cms.InputTag( "hltEcalRawToRecHitFacility" ),
    sourceTag = cms.InputTag( "hltEcalRegionalRestFEDs" ),
    splitOutput = cms.bool( True ),
    EBrechitCollection = cms.string( "EcalRecHitsEB" ),
    EErechitCollection = cms.string( "EcalRecHitsEE" ),
    rechitCollection = cms.string( "NotNeededsplitOutputTrue" ),
    cleaningConfig = cms.PSet( 
      tightenCrack_e1_double = cms.double( 2.0 ),
      tightenCrack_e6e2_double = cms.double( 3.0 ),
      e4e1Threshold_endcap = cms.double( 0.3 ),
      tightenCrack_e4e1_single = cms.double( 3.0 ),
      cThreshold_barrel = cms.double( 4.0 ),
      e4e1Threshold_barrel = cms.double( 0.08 ),
      tightenCrack_e1_single = cms.double( 2.0 ),
      e4e1_b_barrel = cms.double( -0.024 ),
      e4e1_a_barrel = cms.double( 0.04 ),
      ignoreOutOfTimeThresh = cms.double( 1000000.0 ),
      cThreshold_endcap = cms.double( 15.0 ),
      e4e1_b_endcap = cms.double( -0.0125 ),
      e4e1_a_endcap = cms.double( 0.02 ),
      e6e2thresh = cms.double( 0.04 ),
      cThreshold_double = cms.double( 10.0 ),
      swissCrossThreshold = cms.double( 0.95 ),
      recHitThreshold = cms.double( 4.0 ),
      useieta85 = cms.bool( True )
    )
)
process.hltHcalDigis = cms.EDProducer( "HcalRawToDigi",
    InputLabel = cms.InputTag( "source" ),
    UnpackCalib = cms.untracked.bool( True ),
    UnpackZDC = cms.untracked.bool( True ),
    firstSample = cms.int32( 0 ),
    lastSample = cms.int32( 9 ),
    FilterDataQuality = cms.bool( True )
)
process.hltHbhereco = cms.EDProducer( "HcalHitReconstructor",
    correctForTimeslew = cms.bool( True ),
    correctForPhaseContainment = cms.bool( True ),
    correctionPhaseNS = cms.double( 13.0 ),
    digiLabel = cms.InputTag( "hltHcalDigis" ),
    correctTiming = cms.bool( False ),
    setNoiseFlags = cms.bool( False ),
    setHSCPFlags = cms.bool( False ),
    setSaturationFlags = cms.bool( False ),
    setTimingTrustFlags = cms.bool( False ),
    setPulseShapeFlags = cms.bool( False ),
    dropZSmarkedPassed = cms.bool( True ),
    firstAuxTS = cms.int32( 4 ),
    firstSample = cms.int32( 4 ),
    samplesToAdd = cms.int32( 4 ),
    tsFromDB = cms.bool( True ),
    useLeakCorrection = cms.bool( False ),
    Subdetector = cms.string( "HBHE" ),
    setTimingShapedCutsFlags = cms.bool( False ),
    digistat = cms.PSet(  ),
    HFInWindowStat = cms.PSet(  ),
    S9S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 24 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      isS8S1 = cms.bool( False ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    S8S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      shortEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      flagsToSkip = cms.int32( 16 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      longEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      long_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      isS8S1 = cms.bool( True ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    PETstat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_R_29 = cms.vdouble( 0.8 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 0 ),
      short_R = cms.vdouble( 0.8 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      long_R_29 = cms.vdouble( 0.8 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_R = cms.vdouble( 0.98 ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    saturationParameters = cms.PSet(  maxADCvalue = cms.int32( 127 ) ),
    timingshapedcutsParameters = cms.PSet( 
      ignorelowest = cms.bool( True ),
      win_offset = cms.double( 0.0 ),
      ignorehighest = cms.bool( False ),
      win_gain = cms.double( 1.0 ),
      tfilterEnvelope = cms.vdouble( 4.0, 12.04, 13.0, 10.56, 23.5, 8.82, 37.0, 7.38, 56.0, 6.3, 81.0, 5.64, 114.5, 5.44, 175.5, 5.38, 350.5, 5.14 )
    ),
    flagParameters = cms.PSet( 
      nominalPedestal = cms.double( 3.0 ),
      hitMultiplicityThreshold = cms.int32( 17 ),
      hitEnergyMinimum = cms.double( 1.0 ),
      pulseShapeParameterSets = cms.VPSet( 
        cms.PSet(  pulseShapeParameters = cms.vdouble( 0.0, 100.0, -50.0, 0.0, -15.0, 0.15 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 100.0, 2000.0, -50.0, 0.0, -5.0, 0.05 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 2000.0, 1000000.0, -50.0, 0.0, 95.0, 0.0 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 0.0 )        )
      )
    ),
    hscpParameters = cms.PSet( 
      slopeMax = cms.double( -0.6 ),
      r1Max = cms.double( 1.0 ),
      r1Min = cms.double( 0.15 ),
      TimingEnergyThreshold = cms.double( 30.0 ),
      slopeMin = cms.double( -1.5 ),
      outerMin = cms.double( 0.0 ),
      outerMax = cms.double( 0.1 ),
      fracLeaderMin = cms.double( 0.4 ),
      r2Min = cms.double( 0.1 ),
      r2Max = cms.double( 0.5 ),
      fracLeaderMax = cms.double( 0.7 )
    ),
    pulseShapeParameters = cms.PSet(  ),
    hfTimingTrustParameters = cms.PSet( 
      hfTimingTrustLevel2 = cms.int32( 4 ),
      hfTimingTrustLevel1 = cms.int32( 1 )
    ),
    firstAuxOffset = cms.int32( 0 )
)
process.hltHfreco = cms.EDProducer( "HcalHitReconstructor",
    correctForTimeslew = cms.bool( False ),
    correctForPhaseContainment = cms.bool( False ),
    correctionPhaseNS = cms.double( 13.0 ),
    digiLabel = cms.InputTag( "hltHcalDigis" ),
    correctTiming = cms.bool( False ),
    setNoiseFlags = cms.bool( True ),
    setHSCPFlags = cms.bool( False ),
    setSaturationFlags = cms.bool( False ),
    setTimingTrustFlags = cms.bool( False ),
    setPulseShapeFlags = cms.bool( False ),
    dropZSmarkedPassed = cms.bool( True ),
    firstAuxTS = cms.int32( 1 ),
    firstSample = cms.int32( 2 ),
    samplesToAdd = cms.int32( 2 ),
    tsFromDB = cms.bool( True ),
    useLeakCorrection = cms.bool( False ),
    Subdetector = cms.string( "HF" ),
    setTimingShapedCutsFlags = cms.bool( False ),
    digistat = cms.PSet( 
      HFdigiflagFirstSample = cms.int32( 1 ),
      HFdigiflagMinEthreshold = cms.double( 40.0 ),
      HFdigiflagSamplesToAdd = cms.int32( 3 ),
      HFdigiflagCoef0 = cms.double( 0.93 ),
      HFdigiflagCoef2 = cms.double( -0.012667 ),
      HFdigiflagCoef1 = cms.double( -0.38275 ),
      HFdigiflagExpectedPeak = cms.int32( 2 )
    ),
    HFInWindowStat = cms.PSet( 
      hflongEthresh = cms.double( 40.0 ),
      hflongMinWindowTime = cms.vdouble( -10.0 ),
      hfshortEthresh = cms.double( 40.0 ),
      hflongMaxWindowTime = cms.vdouble( 10.0 ),
      hfshortMaxWindowTime = cms.vdouble( 10.0 ),
      hfshortMinWindowTime = cms.vdouble( -12.0 )
    ),
    S9S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 24 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      isS8S1 = cms.bool( False ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    S8S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      shortEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      flagsToSkip = cms.int32( 16 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      longEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      long_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      isS8S1 = cms.bool( True ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    PETstat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_R_29 = cms.vdouble( 0.8 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 0 ),
      short_R = cms.vdouble( 0.8 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      long_R_29 = cms.vdouble( 0.8 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_R = cms.vdouble( 0.98 ),
      HcalAcceptSeverityLevel = cms.int32( 9 )
    ),
    saturationParameters = cms.PSet(  maxADCvalue = cms.int32( 127 ) ),
    timingshapedcutsParameters = cms.PSet( 
      ignorelowest = cms.bool( True ),
      win_offset = cms.double( 0.0 ),
      ignorehighest = cms.bool( False ),
      win_gain = cms.double( 1.0 ),
      tfilterEnvelope = cms.vdouble( 4.0, 12.04, 13.0, 10.56, 23.5, 8.82, 37.0, 7.38, 56.0, 6.3, 81.0, 5.64, 114.5, 5.44, 175.5, 5.38, 350.5, 5.14 )
    ),
    flagParameters = cms.PSet( 
      nominalPedestal = cms.double( 3.0 ),
      hitMultiplicityThreshold = cms.int32( 17 ),
      hitEnergyMinimum = cms.double( 1.0 ),
      pulseShapeParameterSets = cms.VPSet( 
        cms.PSet(  pulseShapeParameters = cms.vdouble( 0.0, 100.0, -50.0, 0.0, -15.0, 0.15 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 100.0, 2000.0, -50.0, 0.0, -5.0, 0.05 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 2000.0, 1000000.0, -50.0, 0.0, 95.0, 0.0 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 0.0 )        )
      )
    ),
    hscpParameters = cms.PSet( 
      slopeMax = cms.double( -0.6 ),
      r1Max = cms.double( 1.0 ),
      r1Min = cms.double( 0.15 ),
      TimingEnergyThreshold = cms.double( 30.0 ),
      slopeMin = cms.double( -1.5 ),
      outerMin = cms.double( 0.0 ),
      outerMax = cms.double( 0.1 ),
      fracLeaderMin = cms.double( 0.4 ),
      r2Min = cms.double( 0.1 ),
      r2Max = cms.double( 0.5 ),
      fracLeaderMax = cms.double( 0.7 )
    ),
    pulseShapeParameters = cms.PSet(  ),
    hfTimingTrustParameters = cms.PSet( 
      hfTimingTrustLevel2 = cms.int32( 4 ),
      hfTimingTrustLevel1 = cms.int32( 1 )
    ),
    firstAuxOffset = cms.int32( 0 )
)
process.hltHoreco = cms.EDProducer( "HcalHitReconstructor",
    correctForTimeslew = cms.bool( True ),
    correctForPhaseContainment = cms.bool( True ),
    correctionPhaseNS = cms.double( 13.0 ),
    digiLabel = cms.InputTag( "hltHcalDigis" ),
    correctTiming = cms.bool( False ),
    setNoiseFlags = cms.bool( False ),
    setHSCPFlags = cms.bool( False ),
    setSaturationFlags = cms.bool( False ),
    setTimingTrustFlags = cms.bool( False ),
    setPulseShapeFlags = cms.bool( False ),
    dropZSmarkedPassed = cms.bool( True ),
    firstAuxTS = cms.int32( 4 ),
    firstSample = cms.int32( 4 ),
    samplesToAdd = cms.int32( 4 ),
    tsFromDB = cms.bool( True ),
    useLeakCorrection = cms.bool( False ),
    Subdetector = cms.string( "HO" ),
    setTimingShapedCutsFlags = cms.bool( False ),
    digistat = cms.PSet(  ),
    HFInWindowStat = cms.PSet(  ),
    S9S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      HcalAcceptSeverityLevel = cms.int32( 9 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 24 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_optimumSlope = cms.vdouble( -99999.0, 0.0164905, 0.0238698, 0.0321383, 0.041296, 0.0513428, 0.0622789, 0.0741041, 0.0868186, 0.100422, 0.135313, 0.136289, 0.0589927 ),
      isS8S1 = cms.bool( False )
    ),
    S8S1stat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      HcalAcceptSeverityLevel = cms.int32( 9 ),
      shortEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      flagsToSkip = cms.int32( 16 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      longEnergyParams = cms.vdouble( 40.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 ),
      long_optimumSlope = cms.vdouble( 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ),
      isS8S1 = cms.bool( True )
    ),
    PETstat = cms.PSet( 
      longETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_R_29 = cms.vdouble( 0.8 ),
      HcalAcceptSeverityLevel = cms.int32( 9 ),
      shortEnergyParams = cms.vdouble( 35.1773, 35.37, 35.7933, 36.4472, 37.3317, 38.4468, 39.7925, 41.3688, 43.1757, 45.2132, 47.4813, 49.98, 52.7093 ),
      flagsToSkip = cms.int32( 0 ),
      long_R_29 = cms.vdouble( 0.8 ),
      shortETParams = cms.vdouble( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ),
      short_R = cms.vdouble( 0.8 ),
      longEnergyParams = cms.vdouble( 43.5, 45.7, 48.32, 51.36, 54.82, 58.7, 63.0, 67.72, 72.86, 78.42, 84.4, 90.8, 97.62 ),
      long_R = cms.vdouble( 0.98 )
    ),
    saturationParameters = cms.PSet(  maxADCvalue = cms.int32( 127 ) ),
    timingshapedcutsParameters = cms.PSet( 
      ignorelowest = cms.bool( True ),
      win_offset = cms.double( 0.0 ),
      ignorehighest = cms.bool( False ),
      win_gain = cms.double( 1.0 ),
      tfilterEnvelope = cms.vdouble( 4.0, 12.04, 13.0, 10.56, 23.5, 8.82, 37.0, 7.38, 56.0, 6.3, 81.0, 5.64, 114.5, 5.44, 175.5, 5.38, 350.5, 5.14 )
    ),
    flagParameters = cms.PSet( 
      nominalPedestal = cms.double( 3.0 ),
      hitMultiplicityThreshold = cms.int32( 17 ),
      hitEnergyMinimum = cms.double( 1.0 ),
      pulseShapeParameterSets = cms.VPSet( 
        cms.PSet(  pulseShapeParameters = cms.vdouble( 0.0, 100.0, -50.0, 0.0, -15.0, 0.15 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 100.0, 2000.0, -50.0, 0.0, -5.0, 0.05 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( 2000.0, 1000000.0, -50.0, 0.0, 95.0, 0.0 )        ),
        cms.PSet(  pulseShapeParameters = cms.vdouble( -1000000.0, 1000000.0, 45.0, 0.1, 1000000.0, 0.0 )        )
      )
    ),
    hscpParameters = cms.PSet( 
      slopeMax = cms.double( -0.6 ),
      r1Max = cms.double( 1.0 ),
      r1Min = cms.double( 0.15 ),
      TimingEnergyThreshold = cms.double( 30.0 ),
      slopeMin = cms.double( -1.5 ),
      outerMin = cms.double( 0.0 ),
      outerMax = cms.double( 0.1 ),
      fracLeaderMin = cms.double( 0.4 ),
      r2Min = cms.double( 0.1 ),
      r2Max = cms.double( 0.5 ),
      fracLeaderMax = cms.double( 0.7 )
    ),
    pulseShapeParameters = cms.PSet(  ),
    hfTimingTrustParameters = cms.PSet( 
      hfTimingTrustLevel2 = cms.int32( 4 ),
      hfTimingTrustLevel1 = cms.int32( 1 )
    ),
    firstAuxOffset = cms.int32( 0 )
)
process.hltTowerMakerForAll = cms.EDProducer( "CaloTowersCreator",
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
process.hltIslandBasicClustersHI = cms.EDProducer( "IslandClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    barrelHitProducer = cms.string( "hltEcalRecHitAll" ),
    endcapHitProducer = cms.string( "hltEcalRecHitAll" ),
    barrelHitCollection = cms.string( "EcalRecHitsEB" ),
    endcapHitCollection = cms.string( "EcalRecHitsEE" ),
    barrelClusterCollection = cms.string( "islandBarrelBasicClustersHI" ),
    endcapClusterCollection = cms.string( "islandEndcapBasicClustersHI" ),
    IslandBarrelSeedThr = cms.double( 0.5 ),
    IslandEndcapSeedThr = cms.double( 0.18 ),
    clustershapecollectionEB = cms.string( "islandBarrelShape" ),
    clustershapecollectionEE = cms.string( "islandEndcapShape" ),
    barrelShapeAssociation = cms.string( "islandBarrelShapeAssoc" ),
    endcapShapeAssociation = cms.string( "islandEndcapShapeAssoc" ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    )
)
process.hltIslandSuperClustersHI = cms.EDProducer( "SuperClusterProducer",
    VerbosityLevel = cms.string( "ERROR" ),
    endcapClusterProducer = cms.string( "hltIslandBasicClustersHI" ),
    barrelClusterProducer = cms.string( "hltIslandBasicClustersHI" ),
    endcapClusterCollection = cms.string( "islandEndcapBasicClustersHI" ),
    barrelClusterCollection = cms.string( "islandBarrelBasicClustersHI" ),
    endcapSuperclusterCollection = cms.string( "islandEndcapSuperClustersHI" ),
    barrelSuperclusterCollection = cms.string( "islandBarrelSuperClustersHI" ),
    doBarrel = cms.bool( True ),
    doEndcaps = cms.bool( True ),
    barrelEtaSearchRoad = cms.double( 0.06 ),
    barrelPhiSearchRoad = cms.double( 0.8 ),
    endcapEtaSearchRoad = cms.double( 0.14 ),
    endcapPhiSearchRoad = cms.double( 0.6 ),
    seedTransverseEnergyThreshold = cms.double( 1.0 )
)
process.hltCorrectedIslandBarrelSuperClustersHI = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    rawSuperClusterProducer = cms.InputTag( 'hltIslandSuperClustersHI','islandBarrelSuperClustersHI' ),
    superClusterAlgo = cms.string( "Island" ),
    applyEnergyCorrection = cms.bool( False ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 0.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet(  ),
    isl_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 0.0 ),
      fEtEtaVec = cms.vdouble( 0.0 ),
      brLinearHighThr = cms.double( 0.0 ),
      fBremVec = cms.vdouble( 0.0 )
    ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet(  )
)
process.hltCorrectedIslandEndcapSuperClustersHI = cms.EDProducer( "EgammaSCCorrectionMaker",
    VerbosityLevel = cms.string( "ERROR" ),
    recHitProducer = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' ),
    rawSuperClusterProducer = cms.InputTag( 'hltIslandSuperClustersHI','islandEndcapSuperClustersHI' ),
    superClusterAlgo = cms.string( "Island" ),
    applyEnergyCorrection = cms.bool( False ),
    applyCrackCorrection = cms.bool( False ),
    sigmaElectronicNoise = cms.double( 0.15 ),
    etThresh = cms.double( 0.0 ),
    corectedSuperClusterCollection = cms.string( "" ),
    hyb_fCorrPset = cms.PSet(  ),
    isl_fCorrPset = cms.PSet( 
      brLinearLowThr = cms.double( 0.0 ),
      fEtEtaVec = cms.vdouble( 0.0 ),
      brLinearHighThr = cms.double( 0.0 ),
      fBremVec = cms.vdouble( 0.0 )
    ),
    dyn_fCorrPset = cms.PSet(  ),
    fix_fCorrPset = cms.PSet(  )
)
process.hltCleanedCorrectedIslandBarrelSuperClustersHI = cms.EDProducer( "HiSpikeCleaner",
    recHitProducerBarrel = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
    recHitProducerEndcap = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' ),
    originalSuperClusterProducer = cms.InputTag( "hltCorrectedIslandBarrelSuperClustersHI" ),
    TimingCut = cms.untracked.double( 9999999.0 ),
    etCut = cms.double( 8.0 ),
    outputColl = cms.string( "" )
)
process.hltRecoHIEcalWithCleaningCandidate = cms.EDProducer( "EgammaHLTRecoEcalCandidateProducers",
    scHybridBarrelProducer = cms.InputTag( "hltCleanedCorrectedIslandBarrelSuperClustersHI" ),
    scIslandEndcapProducer = cms.InputTag( "hltCorrectedIslandEndcapSuperClustersHI" ),
    recoEcalCandidateCollection = cms.string( "" )
)
process.hltHIPhoton15 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIPhoton20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIPhoton20 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIPhoton30 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIPhoton30 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 30.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIPhoton40 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIPhoton40 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 40.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIDoublePhoton10and15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoublePhoton1015Filter1 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 10.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 2 )
)
process.hltHIDoublePhoton1015Filter2 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIDoublePhoton15and20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoublePhoton1520Filter1 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 2 )
)
process.hltHIDoublePhoton1520Filter2 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 1 )
)
process.hltPreHIDoublePhoton10 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoublePhoton10 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 10.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 2 )
)
process.hltPreHIDoublePhoton15 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoublePhoton15 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 15.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 2 )
)
process.hltPreHIDoublePhoton20 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIDoublePhoton20 = cms.EDFilter( "HLT1Photon",
    inputTag = cms.InputTag( "hltRecoHIEcalWithCleaningCandidate" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 20.0 ),
    MaxEta = cms.double( 2.0 ),
    MinN = cms.int32( 2 )
)
process.hltBPTXAntiCoincidence = cms.EDFilter( "HLTLevel1Activity",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    daqPartitions = cms.uint32( 1 ),
    ignoreL1Mask = cms.bool( True ),
    invert = cms.bool( True ),
    bunchCrossings = cms.vint32( 0, 1, -1 ),
    physicsLoBits = cms.uint64( 0x1 ),
    physicsHiBits = cms.uint64( 0x0 ),
    technicalBits = cms.uint64( 0x11 )
)
process.hltL1sL1SingleJet20UNotBptxOR = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleJet20U_NotBptxOR" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( False )
)
process.hltPreHIStoppedHSCP35 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltStoppedHSCPHpdFilter = cms.EDFilter( "HLTHPDFilter",
    inputTag = cms.InputTag( "hltHbhereco" ),
    energy = cms.double( -99.0 ),
    hpdSpikeEnergy = cms.double( 10.0 ),
    hpdSpikeIsolationEnergy = cms.double( 1.0 ),
    rbxSpikeEnergy = cms.double( 50.0 ),
    rbxSpikeUnbalance = cms.double( 0.2 )
)
process.hltStoppedHSCPTowerMakerForAll = cms.EDProducer( "CaloTowersCreator",
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
    hoInput = cms.InputTag( "" ),
    hfInput = cms.InputTag( "" ),
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
process.hltStoppedHSCPIterativeCone5CaloJets = cms.EDProducer( "FastjetJetProducer",
    UseOnlyVertexTracks = cms.bool( False ),
    UseOnlyOnePV = cms.bool( False ),
    DzTrVtxMax = cms.double( 0.0 ),
    DxyTrVtxMax = cms.double( 0.0 ),
    MinVtxNdof = cms.int32( 5 ),
    MaxVtxZ = cms.double( 15.0 ),
    jetAlgorithm = cms.string( "IterativeCone" ),
    rParam = cms.double( 0.5 ),
    src = cms.InputTag( "hltStoppedHSCPTowerMakerForAll" ),
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
    sumRecHits = cms.bool( False ),
    doAreaDiskApprox = cms.bool( False ),
    voronoiRfact = cms.double( -9.0 ),
    useDeterministicSeed = cms.bool( False ),
    minSeed = cms.uint32( 0 ),
    Rho_EtaMax = cms.double( 4.4 ),
    puPtMin = cms.double( 10.0 )
)
process.hltStoppedHSCP1CaloJetEnergy35 = cms.EDFilter( "HLT1CaloJetEnergy",
    inputTag = cms.InputTag( "hltStoppedHSCPIterativeCone5CaloJets" ),
    saveTags = cms.bool( True ),
    MinE = cms.double( 35.0 ),
    MaxEta = cms.double( 3.0 ),
    MinN = cms.int32( 1 )
)
process.hltL1sL1GlobalDecision = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1GlobalDecision" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIActivityHFCoincidence3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHcalSimpleRecHitFilterCoincidence = cms.EDFilter( "HLTHcalSimpleRecHitFilter",
    threshold = cms.double( 3.0 ),
    minNHitsNeg = cms.int32( 1 ),
    minNHitsPos = cms.int32( 1 ),
    doCoincidence = cms.bool( True ),
    HFRecHitCollection = cms.InputTag( "hltHfreco" ),
    maskedChannels = cms.vint32(  )
)
process.hltPreHIActivityHFSingle3 = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHcalSimpleRecHitFilter = cms.EDFilter( "HLTHcalSimpleRecHitFilter",
    threshold = cms.double( 3.0 ),
    minNHitsNeg = cms.int32( 1 ),
    minNHitsPos = cms.int32( 1 ),
    doCoincidence = cms.bool( False ),
    HFRecHitCollection = cms.InputTag( "hltHfreco" ),
    maskedChannels = cms.vint32(  )
)
process.hltPreHIClusterVertexCompatibility = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHIPixelClusterShapeFilter = cms.EDFilter( "HLTPixelClusterShapeFilter",
    inputTag = cms.InputTag( "hltHISiPixelRecHits" ),
    saveTags = cms.bool( False ),
    minZ = cms.double( -20.0 ),
    maxZ = cms.double( 20.05 ),
    zStep = cms.double( 0.2 ),
    nhitsTrunc = cms.int32( 150 ),
    clusterTrunc = cms.double( 2.0 ),
    clusterPars = cms.vdouble(  )
)
process.hltPreHICentralityVeto = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPixelActivityFilterCentralityVeto = cms.EDFilter( "HLTPixelActivityFilter",
    inputTag = cms.InputTag( "hltHISiPixelClusters" ),
    saveTags = cms.bool( False ),
    minClusters = cms.uint32( 3 ),
    maxClusters = cms.uint32( 1000 )
)
process.hltBPTXCoincidence = cms.EDFilter( "HLTLevel1Activity",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    daqPartitions = cms.uint32( 1 ),
    ignoreL1Mask = cms.bool( True ),
    invert = cms.bool( False ),
    bunchCrossings = cms.vint32( 0, -1, 1 ),
    physicsLoBits = cms.uint64( 0x1 ),
    physicsHiBits = cms.uint64( 0x40000 ),
    technicalBits = cms.uint64( 0x1 )
)
process.hltL1sETT100 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_ETT100" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIFullTrack12L1Central = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHICaloTowerFilter4 = cms.EDFilter( "HLTCaloTowerFilter",
    inputTag = cms.InputTag( "hltTowerMakerForAll" ),
    saveTags = cms.bool( False ),
    MinPt = cms.double( 4.0 ),
    MaxEta = cms.double( 2.4 ),
    MinN = cms.uint32( 1 )
)
process.hltHIPixelClusterVerticesForHITrackTrigger = cms.EDProducer( "HIPixelClusterVtxProducer",
    pixelRecHits = cms.string( "hltHISiPixelRecHits" ),
    minZ = cms.double( -15.0 ),
    maxZ = cms.double( 15.0 ),
    zStep = cms.double( 1.0 )
)
process.hltHIPixel3ProtoTracks = cms.EDProducer( "PixelTrackProducer",
    useFilterWithES = cms.bool( False ),
    RegionFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "HITrackingRegionForPrimaryVtxProducer" ),
      RegionPSet = cms.PSet( 
        precise = cms.bool( True ),
        originRadius = cms.double( 0.1 ),
        ptMin = cms.double( 1.0 ),
        directionXCoord = cms.double( 1.0 ),
        directionZCoord = cms.double( 0.0 ),
        directionYCoord = cms.double( 1.0 ),
        useFoundVertices = cms.bool( True ),
        doVariablePtMin = cms.bool( True ),
        nSigmaZ = cms.double( 3.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        useFixedError = cms.bool( True ),
        fixedError = cms.double( 1.0 ),
        sigmaZVertex = cms.double( 3.0 ),
        siPixelRecHits = cms.InputTag( "hltHISiPixelRecHits" ),
        VertexCollection = cms.InputTag( "hltHIPixelClusterVerticesForHITrackTrigger" )
      )
    ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitTripletGenerator" ),
      SeedingLayers = cms.string( "hltESPHIPixelLayerTriplets" ),
      GeneratorPSet = cms.PSet( 
        useBending = cms.bool( True ),
        useFixedPreFiltering = cms.bool( False ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        useMultScattering = cms.bool( True ),
        ComponentName = cms.string( "PixelTripletHLTGenerator" ),
        extraHitRZtolerance = cms.double( 0.037 ),
        maxElement = cms.uint32( 100000 ),
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
      )
    ),
    FitterPSet = cms.PSet( 
      ComponentName = cms.string( "PixelFitterByHelixProjections" ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderWithoutAngle4PixelTriplets" )
    ),
    FilterPSet = cms.PSet( 
      chi2 = cms.double( 1000.0 ),
      ComponentName = cms.string( "HIProtoTrackFilter" ),
      ptMin = cms.double( 1.0 ),
      tipMax = cms.double( 1.0 ),
      doVariablePtMin = cms.bool( True ),
      beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      siPixelRecHits = cms.InputTag( "hltHISiPixelRecHits" )
    ),
    CleanerPSet = cms.PSet(  ComponentName = cms.string( "PixelTrackCleanerBySharedHits" ) )
)
process.hltHIPixelMedianVertex = cms.EDProducer( "HIPixelMedianVtxProducer",
    TrackCollection = cms.InputTag( "hltHIPixel3ProtoTracks" ),
    PtMin = cms.double( 0.075 ),
    PeakFindThreshold = cms.uint32( 100 ),
    PeakFindMaxZ = cms.double( 15.0 ),
    PeakFindBinsPerCm = cms.int32( 10 ),
    FitThreshold = cms.int32( 5 ),
    FitMaxZ = cms.double( 0.1 ),
    FitBinsPerCm = cms.int32( 500 )
)
process.hltHISelectedProtoTracks = cms.EDFilter( "HIProtoTrackSelection",
    src = cms.InputTag( "hltHIPixel3ProtoTracks" ),
    VertexCollection = cms.InputTag( "hltHIPixelMedianVertex" ),
    beamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
    ptMin = cms.double( 0.0 ),
    nSigmaZ = cms.double( 5.0 ),
    minZCut = cms.double( 0.2 ),
    maxD0Significance = cms.double( 5.0 )
)
process.hltHIPixelAdaptiveVertex = cms.EDProducer( "PrimaryVertexProducer",
    TrackLabel = cms.InputTag( "hltHISelectedProtoTracks" ),
    beamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
    algorithm = cms.string( "AdaptiveVertexFitter" ),
    label = cms.string( "" ),
    minNdof = cms.double( 0.0 ),
    useBeamConstraint = cms.bool( False ),
    maxDistanceToBeam = cms.double( 2.0 ),
    TkFilterParameters = cms.PSet( 
      maxD0Significance = cms.double( 3.0 ),
      minPt = cms.double( 0.0 ),
      maxNormalizedChi2 = cms.double( 5.0 ),
      minSiliconLayersWithHits = cms.int32( 0 ),
      minPixelLayersWithHits = cms.int32( 2 ),
      trackQuality = cms.string( "any" ),
      numTracksThreshold = cms.int32( 2 ),
      algorithm = cms.string( "filterWithThreshold" )
    ),
    PVSelParameters = cms.PSet(  maxDistanceToBeam = cms.double( 0.1 ) ),
    TkClusParameters = cms.PSet( 
      algorithm = cms.string( "gap" ),
      TkGapClusParameters = cms.PSet(  zSeparation = cms.double( 1.0 ) )
    ),
    vertexCollections = cms.VPSet( 
      cms.PSet(  maxDistanceToBeam = cms.double( 2.0 ),
        useBeamConstraint = cms.bool( False ),
        minNdof = cms.double( 0.0 ),
        algorithm = cms.string( "AdaptiveVertexFitter" ),
        label = cms.string( "" )
      )
    )
)
process.hltHIBestAdaptiveVertex = cms.EDFilter( "HIBestVertexSelection",
    src = cms.InputTag( "hltHIPixelAdaptiveVertex" ),
    maxNumber = cms.uint32( 1 )
)
process.hltHISelectedVertex = cms.EDProducer( "HIBestVertexProducer",
    medianVertexCollection = cms.InputTag( "hltHIPixelMedianVertex" ),
    adaptiveVertexCollection = cms.InputTag( "hltHIBestAdaptiveVertex" )
)
process.hltHIPixel3PrimTracks = cms.EDProducer( "PixelTrackProducer",
    useFilterWithES = cms.bool( True ),
    RegionFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "GlobalTrackingRegionWithVerticesProducer" ),
      RegionPSet = cms.PSet( 
        precise = cms.bool( True ),
        originRadius = cms.double( 0.05 ),
        ptMin = cms.double( 0.9 ),
        useFoundVertices = cms.bool( True ),
        nSigmaZ = cms.double( 2.0 ),
        beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
        useFixedError = cms.bool( True ),
        fixedError = cms.double( 0.2 ),
        sigmaZVertex = cms.double( 3.0 ),
        VertexCollection = cms.InputTag( "hltHISelectedVertex" )
      )
    ),
    OrderedHitsFactoryPSet = cms.PSet( 
      ComponentName = cms.string( "StandardHitTripletGenerator" ),
      SeedingLayers = cms.string( "hltESPHIPixelLayerTriplets" ),
      GeneratorPSet = cms.PSet( 
        useBending = cms.bool( True ),
        useFixedPreFiltering = cms.bool( False ),
        phiPreFiltering = cms.double( 0.3 ),
        extraHitRPhitolerance = cms.double( 0.032 ),
        useMultScattering = cms.bool( True ),
        ComponentName = cms.string( "PixelTripletHLTGenerator" ),
        extraHitRZtolerance = cms.double( 0.037 ),
        maxElement = cms.uint32( 100000 ),
        SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) )
      )
    ),
    FitterPSet = cms.PSet( 
      ComponentName = cms.string( "PixelFitterByHelixProjections" ),
      TTRHBuilder = cms.string( "hltESPTTRHBuilderWithoutAngle4PixelTriplets" )
    ),
    FilterPSet = cms.PSet( 
      chi2 = cms.double( 1000.0 ),
      ComponentName = cms.string( "HIPixelTrackFilter" ),
      ptMin = cms.double( 0.9 ),
      tipMax = cms.double( 0.0 ),
      useClusterShape = cms.bool( False ),
      VertexCollection = cms.InputTag( "hltHISelectedVertex" ),
      nSigmaTipMaxTolerance = cms.double( 6.0 ),
      nSigmaLipMaxTolerance = cms.double( 0.0 ),
      lipMax = cms.double( 0.3 )
    ),
    CleanerPSet = cms.PSet(  ComponentName = cms.string( "TrackCleaner" ) )
)
process.hltHIPixelTrackCandsForHITrackTrigger = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltHIPixel3PrimTracks" ),
    particleType = cms.string( "pi+" )
)
process.hltHISinglePixelTrackFilter10 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltHISelectedVertex" ),
    trackCollection = cms.InputTag( "hltHIPixelTrackCandsForHITrackTrigger" ),
    MinPt = cms.double( 10.0 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 15.0 ),
    MinTrks = cms.int32( 1 ),
    MinSep = cms.double( 0.2 )
)
process.hltSiStripRawToDigi = cms.EDProducer( "SiStripRawToDigiModule",
    ProductLabel = cms.InputTag( "source" ),
    AppendedBytes = cms.int32( 0 ),
    UseDaqRegister = cms.bool( False ),
    UseFedKey = cms.bool( False ),
    UnpackBadChannels = cms.bool( False ),
    MarkModulesOnMissingFeds = cms.bool( True ),
    TriggerFedId = cms.int32( 0 ),
    UnpackCommonModeValues = cms.bool( False ),
    DoAllCorruptBufferChecks = cms.bool( False ),
    ErrorThreshold = cms.uint32( 7174 )
)
process.hltSiStripZeroSuppression = cms.EDProducer( "SiStripZeroSuppression",
    fixCM = cms.bool( False ),
    produceCalculatedBaseline = cms.bool( False ),
    doAPVRestore = cms.bool( True ),
    produceBaselinePoints = cms.bool( False ),
    useCMMeanMap = cms.bool( False ),
    RawDigiProducersList = cms.VInputTag( 'hltSiStripRawToDigi:VirginRaw','hltSiStripRawToDigi:ProcessedRaw','hltSiStripRawToDigi:ScopeMode' ),
    storeInZScollBadAPV = cms.bool( True ),
    mergeCollections = cms.bool( False ),
    Algorithms = cms.PSet( 
      CutToAvoidSignal = cms.double( 2.0 ),
      PedestalSubtractionFedMode = cms.bool( False ),
      Fraction = cms.double( 0.2 ),
      minStripsToFit = cms.uint32( 4 ),
      consecThreshold = cms.uint32( 5 ),
      hitStripThreshold = cms.uint32( 40 ),
      Deviation = cms.uint32( 25 ),
      CommonModeNoiseSubtractionMode = cms.string( "IteratedMedian" ),
      TruncateInSuppressor = cms.bool( True ),
      restoreThreshold = cms.double( 0.5 ),
      APVInspectMode = cms.string( "BaselineFollower" ),
      ForceNoRestore = cms.bool( False ),
      useRealMeanCM = cms.bool( False ),
      DeltaCMThreshold = cms.uint32( 20 ),
      nSigmaNoiseDerTh = cms.uint32( 4 ),
      nSaturatedStrip = cms.uint32( 2 ),
      SiStripFedZeroSuppressionMode = cms.uint32( 4 ),
      APVRestoreMode = cms.string( "BaselineFollower" ),
      distortionThreshold = cms.uint32( 20 ),
      Iterations = cms.int32( 3 ),
      nSmooth = cms.uint32( 9 ),
      SelfSelectRestoreAlgo = cms.bool( False ),
      doAPVRestore = cms.bool( True ),
      useCMMeanMap = cms.bool( False ),
      ApplyBaselineCleaner = cms.bool( True ),
      MeanCM = cms.int32( 0 )
    ),
    storeCM = cms.bool( False ),
    produceRawDigis = cms.bool( True ),
    NumberOfModuleSkippedBeforeForcingRestore = cms.uint32( 20 )
)
process.hltHISiStripClustersNonRegional = cms.EDProducer( "SiStripClusterizer",
    DigiProducersList = cms.VInputTag( 'hltSiStripZeroSuppression:VirginRaw','hltSiStripZeroSuppression:ProcessedRaw','hltSiStripZeroSuppression:ScopeMode' ),
    Clusterizer = cms.PSet( 
      ChannelThreshold = cms.double( 2.0 ),
      MaxSequentialBad = cms.uint32( 1 ),
      Algorithm = cms.string( "ThreeThresholdAlgorithm" ),
      MaxSequentialHoles = cms.uint32( 0 ),
      MaxAdjacentBad = cms.uint32( 0 ),
      QualityLabel = cms.string( "" ),
      SeedThreshold = cms.double( 3.0 ),
      ClusterThreshold = cms.double( 5.0 )
    )
)
process.hltHIPixelTrackSeeds = cms.EDProducer( "SeedGeneratorFromProtoTracksEDProducer",
    InputCollection = cms.InputTag( "hltHIPixel3PrimTracks" ),
    InputVertexCollection = cms.InputTag( "" ),
    originHalfLength = cms.double( 1.0E9 ),
    originRadius = cms.double( 1.0E9 ),
    useProtoTrackKinematics = cms.bool( False ),
    useEventsWithNoVertex = cms.bool( True ),
    TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" )
)
process.hltHIPrimTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    src = cms.InputTag( "hltHIPixelTrackSeeds" ),
    TrajectoryBuilder = cms.string( "hltESPCkfTrajectoryBuilderForHI" ),
    TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedSeeds" ),
    NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    RedundantSeedCleaner = cms.string( "none" ),
    useHitsSplitting = cms.bool( True ),
    doSeedingRegionRebuilding = cms.bool( False ),
    TransientInitialStateEstimatorParameters = cms.PSet( 
      propagatorAlongTISE = cms.string( "PropagatorWithMaterialForHI" ),
      propagatorOppositeTISE = cms.string( "PropagatorWithMaterialOppositeForHI" ),
      numberMeasurementsForFit = cms.int32( 4 )
    ),
    cleanTrajectoryAfterInOut = cms.bool( True ),
    maxNSeeds = cms.uint32( 100000 )
)
process.hltHIGlobalPrimTracks = cms.EDProducer( "TrackProducer",
    TrajectoryInEvent = cms.bool( True ),
    useHitsSplitting = cms.bool( True ),
    clusterRemovalInfo = cms.InputTag( "" ),
    alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    Fitter = cms.string( "hltESPKFFittingSmootherWithOutliersRejectionAndRK" ),
    Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" ),
    src = cms.InputTag( "hltHIPrimTrackCandidates" ),
    beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    TTRHBuilder = cms.string( "hltESPTTRHBuilderAngleAndTemplate" ),
    AlgorithmName = cms.string( "PropagatorWithMaterialOppositeForHI" ),
    NavigationSchool = cms.string( "" )
)
process.hltHIGoodLooseTracks = cms.EDProducer( "AnalyticalTrackSelector",
    src = cms.InputTag( "hltHIGlobalPrimTracks" ),
    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    useVertices = cms.bool( True ),
    useVtxError = cms.bool( True ),
    vertices = cms.InputTag( "hltHISelectedVertex" ),
    copyTrajectories = cms.untracked.bool( True ),
    vtxNumber = cms.int32( -1 ),
    vertexCut = cms.string( "" ),
    chi2n_par = cms.double( 9999.0 ),
    chi2n_no1Dmod_par = cms.double( 0.2 ),
    applyAdaptedPVCuts = cms.bool( True ),
    max_d0 = cms.double( 100.0 ),
    max_z0 = cms.double( 100.0 ),
    nSigmaZ = cms.double( 4.0 ),
    minNumberLayers = cms.uint32( 6 ),
    minNumber3DLayers = cms.uint32( 0 ),
    maxNumberLostLayers = cms.uint32( 999 ),
    max_relpterr = cms.double( 0.1 ),
    min_nhits = cms.uint32( 10 ),
    applyAbsCutsIfNoPV = cms.bool( False ),
    keepAllTracks = cms.bool( False ),
    qualityBit = cms.string( "loose" ),
    max_d0NoPV = cms.double( 0.2 ),
    max_z0NoPV = cms.double( 15.0 ),
    res_par = cms.vdouble( 99999.0, 99999.0 ),
    d0_par1 = cms.vdouble( 9999.0, 0.0 ),
    dz_par1 = cms.vdouble( 9999.0, 0.0 ),
    d0_par2 = cms.vdouble( 8.0, 0.0 ),
    dz_par2 = cms.vdouble( 8.0, 0.0 )
)
process.hltHIFullTrackCandsForHITrackTrigger = cms.EDProducer( "ConcreteChargedCandidateProducer",
    src = cms.InputTag( "hltHIGoodLooseTracks" ),
    particleType = cms.string( "pi+" )
)
process.hltHISingleFullTrackFilter12 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltHISelectedVertex" ),
    trackCollection = cms.InputTag( "hltHIFullTrackCandsForHITrackTrigger" ),
    MinPt = cms.double( 12.0 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 15.0 ),
    MinTrks = cms.int32( 1 ),
    MinSep = cms.double( 0.2 )
)
process.hltPreHIFullTrack14L1Central = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHISingleFullTrackFilter14 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltHISelectedVertex" ),
    trackCollection = cms.InputTag( "hltHIFullTrackCandsForHITrackTrigger" ),
    MinPt = cms.double( 14.0 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 15.0 ),
    MinTrks = cms.int32( 1 ),
    MinSep = cms.double( 0.2 )
)
process.hltPreHIFullTrack20L1Central = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHISingleFullTrackFilter20 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltHISelectedVertex" ),
    trackCollection = cms.InputTag( "hltHIFullTrackCandsForHITrackTrigger" ),
    MinPt = cms.double( 20.0 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 15.0 ),
    MinTrks = cms.int32( 1 ),
    MinSep = cms.double( 0.2 )
)
process.hltPreHIFullTrack25L1Central = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltHISingleFullTrackFilter25 = cms.EDFilter( "HLTSingleVertexPixelTrackFilter",
    vertexCollection = cms.InputTag( "hltHISelectedVertex" ),
    trackCollection = cms.InputTag( "hltHIFullTrackCandsForHITrackTrigger" ),
    MinPt = cms.double( 25.0 ),
    MaxPt = cms.double( 10000.0 ),
    MaxEta = cms.double( 2.4 ),
    MaxVz = cms.double( 15.0 ),
    MinTrks = cms.int32( 1 ),
    MinSep = cms.double( 0.2 )
)
process.hltRandomEventsFilter = cms.EDFilter( "HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32( 3 )
)
process.hltPreHIRandom = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuMu = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_BscMinBiasThreshold2_BptxAND) AND L1_SingleMuOpen" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuMuPixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuEG2 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_BscMinBiasThreshold2_BptxAND) AND L1_SingleEG2" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuEG2PixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_BscMinBiasThreshold2_BptxAND) AND L1_SingleEG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuEG5PixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuHcalHfMu = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_HcalHfCoincidencePm_BptxAND ) AND L1_SingleMuOpen" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuHcalHfMuPixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuHcalHfEG2 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_HcalHfCoincidencePm_BptxAND) AND L1_SingleEG2" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuHcalHfEG2PixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltL1sHIUPCNeuHcalHfEG5 = cms.EDFilter( "HLTLevel1GTSeed",
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1TechTriggerSeeding = cms.bool( False ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "(L1_ZdcCaloPlus_BptxAND OR L1_ZdcCaloMinus_BptxAND) AND (NOT L1_HcalHfCoincidencePm_BptxAND) AND L1_SingleEG5" ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    saveTags = cms.bool( True )
)
process.hltPreHIUPCNeuHcalHfEG5PixelSingleTrack = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltFEDSelector = cms.EDProducer( "EvFFEDSelector",
    inputTag = cms.InputTag( "source" ),
    fedList = cms.vuint32( 1023 )
)
process.hltTriggerSummaryAOD = cms.EDProducer( "TriggerSummaryProducerAOD",
    processName = cms.string( "@" )
)
process.hltTriggerSummaryRAW = cms.EDProducer( "TriggerSummaryProducerRAW",
    processName = cms.string( "@" )
)
process.hltL1GtTrigReport = cms.EDAnalyzer( "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "hltGtDigis" )
)
process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
)
process.hltPreALCAP0Output = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreALCAPHISYMOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreCalibrationOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltDQML1Scalers = cms.EDAnalyzer( "L1Scalers",
    l1GtData = cms.InputTag( "hltGtDigis" ),
    fedRawData = cms.InputTag( "source" ),
    HFRecHitCollection = cms.InputTag( "hltHfreco" ),
    maskedChannels = cms.untracked.vint32( 8137, 8141, 8146, 8149, 8150, 8153 )
)
process.hltDQML1SeedLogicScalers = cms.EDAnalyzer( "HLTSeedL1LogicScalers",
    l1BeforeMask = cms.bool( False ),
    processname = cms.string( "HLT" ),
    L1GtDaqReadoutRecordInputTag = cms.InputTag( "hltGtDigis" ),
    L1GtRecordInputTag = cms.InputTag( "unused" ),
    DQMFolder = cms.untracked.string( "HLT/HLTSeedL1LogicScalers_EvF" ),
    monitorPaths = cms.untracked.vstring( 'HLT_L1MuOpen',
      'HLT_L1Mu',
      'HLT_Mu3',
      'HLT_L1SingleForJet',
      'HLT_SingleLooseIsoTau20',
      'HLT_MinBiasEcal' )
)
process.hltDQMHLTScalers = cms.EDAnalyzer( "HLTScalers",
    processname = cms.string( "HLT" ),
    triggerResults = cms.InputTag( 'TriggerResults','','HLT' )
)
process.hltPreDQMOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreDQMOutputSmart = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'HLT_DTCalibration_v1 / 10',
      'HLT_EcalCalibration_v2 / 10',
      'HLT_HcalCalibration_v2',
      'HLT_TrackerCalibration_v2',
      'HLT_LogMonitor_v1',
      'HLT_HIFullTrack12_L1Central_v1',
      'HLT_HIFullTrack14_L1Central_v1',
      'HLT_HIFullTrack20_L1Central_v1',
      'HLT_HIFullTrack25_L1Central_v1',
      'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1' ),
    hltResults = cms.InputTag( "TriggerResults" ),
    l1tResults = cms.InputTag( "hltGtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( False )
)
process.hltPreEcalCalibrationOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreExpressOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreExpressOutputSmart = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(  ),
    hltResults = cms.InputTag( "TriggerResults" ),
    l1tResults = cms.InputTag( "hltGtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( False )
)
process.hltPreHLTDQMOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreHLTDQMOutputSmart = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'HLT_LogMonitor_v1',
      'HLT_HIFullTrack12_L1Central_v1',
      'HLT_HIFullTrack14_L1Central_v1',
      'HLT_HIFullTrack20_L1Central_v1',
      'HLT_HIFullTrack25_L1Central_v1',
      'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1' ),
    hltResults = cms.InputTag( "TriggerResults" ),
    l1tResults = cms.InputTag( "hltGtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( False )
)
process.hltPreHLTDQMResultsOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreHLTDQMResultsOutputSmart = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'HLT_* AND NOT HLT_*Calibration*' ),
    hltResults = cms.InputTag( "TriggerResults" ),
    l1tResults = cms.InputTag( "hltGtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( False )
)
process.hltPreHLTMONOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreHLTMONOutputSmart = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring( 'HLT_HcalCalibration_v2',
      'HLT_LogMonitor_v1',
      'HLT_HIFullTrack12_L1Central_v1',
      'HLT_HIFullTrack14_L1Central_v1',
      'HLT_HIFullTrack20_L1Central_v1',
      'HLT_HIFullTrack25_L1Central_v1',
      'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
      'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1' ),
    hltResults = cms.InputTag( "TriggerResults" ),
    l1tResults = cms.InputTag( "hltGtDigis" ),
    l1tIgnoreMask = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True ),
    l1techIgnorePrescales = cms.bool( False )
)
process.hltPreNanoDSTOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreRPCMONOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltPreTrackerCalibrationOutput = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)

process.hltOutputA = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputA.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
  'HLT_HIActivityHF_Single3_v1',
  'HLT_HIBptxXOR_v1',
  'HLT_HICentralityVeto_v1',
  'HLT_HIClusterVertexCompatibility_v1',
  'HLT_HIDoublePhoton10_v1',
  'HLT_HIDoublePhoton15_v1',
  'HLT_HIDoublePhoton20_v1',
  'HLT_HIFullTrack12_L1Central_v1',
  'HLT_HIFullTrack14_L1Central_v1',
  'HLT_HIFullTrack20_L1Central_v1',
  'HLT_HIFullTrack25_L1Central_v1',
  'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
  'HLT_HIL1DoubleMuOpen_v1',
  'HLT_HIL1SingleMu3_v1',
  'HLT_HIL1SingleMu5_v1',
  'HLT_HIL1SingleMu7_v1',
  'HLT_HIL2DoubleMu0_v1',
  'HLT_HIL2DoubleMu3_v1',
  'HLT_HIL2Mu20_v1',
  'HLT_HIL2Mu3_v1',
  'HLT_HIL2Mu5Tight_v1',
  'HLT_HIMinBiasBSC_OR_v1',
  'HLT_HIMinBiasBSC_v1',
  'HLT_HIMinBiasHF_v1',
  'HLT_HIMinBiasHfOrBSC_v1',
  'HLT_HIMinBiasHf_OR_v1',
  'HLT_HIMinBiasPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
  'HLT_HIMinBiasZDC_Calo_v1',
  'HLT_HIMinBiasZDC_Scint_v1',
  'HLT_HIPhoton10_Photon15_v1',
  'HLT_HIPhoton15_Photon20_v1',
  'HLT_HIRandom_v1',
  'HLT_HISinglePhoton15_v1',
  'HLT_HISinglePhoton20_v1',
  'HLT_HISinglePhoton30_v1',
  'HLT_HISinglePhoton40_v1',
  'HLT_HIStoppedHSCP35_v1',
  'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
  'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
  'HLT_HIUpcEcal_v1',
  'HLT_HIUpcMu_v1',
  'HLT_HIZeroBiasPixel_SingleTrack_v1',
  'HLT_HIZeroBiasXOR_v1',
  'HLT_HIZeroBias_v1',
  'HLT_LogMonitor_v1' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltL1GtObjectMap_*_*',
      'keep FEDRawDataCollection_rawDataCollector_*_*',
      'keep FEDRawDataCollection_source_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputCalibration = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputCalibration.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_DTCalibration_v1',
  'HLT_EcalCalibration_v2',
  'HLT_HcalCalibration_v2' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltDTCalibrationRaw_*_*',
      'keep *_hltEcalCalibrationRaw_*_*',
      'keep *_hltHcalCalibrationRaw_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputDQM = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputDQM.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_DTCalibration_v1',
  'HLT_EcalCalibration_v2',
  'HLT_HIActivityHF_Coincidence3_v1',
  'HLT_HIActivityHF_Single3_v1',
  'HLT_HIBptxXOR_v1',
  'HLT_HICentralityVeto_v1',
  'HLT_HIClusterVertexCompatibility_v1',
  'HLT_HIDoublePhoton10_v1',
  'HLT_HIDoublePhoton15_v1',
  'HLT_HIDoublePhoton20_v1',
  'HLT_HIFullTrack12_L1Central_v1',
  'HLT_HIFullTrack14_L1Central_v1',
  'HLT_HIFullTrack20_L1Central_v1',
  'HLT_HIFullTrack25_L1Central_v1',
  'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
  'HLT_HIL1DoubleMuOpen_v1',
  'HLT_HIL1SingleMu3_v1',
  'HLT_HIL1SingleMu5_v1',
  'HLT_HIL1SingleMu7_v1',
  'HLT_HIL2DoubleMu0_v1',
  'HLT_HIL2DoubleMu3_v1',
  'HLT_HIL2Mu20_v1',
  'HLT_HIL2Mu3_v1',
  'HLT_HIL2Mu5Tight_v1',
  'HLT_HIMinBiasBSC_OR_v1',
  'HLT_HIMinBiasBSC_v1',
  'HLT_HIMinBiasHF_v1',
  'HLT_HIMinBiasHfOrBSC_v1',
  'HLT_HIMinBiasHf_OR_v1',
  'HLT_HIMinBiasPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
  'HLT_HIMinBiasZDC_Calo_v1',
  'HLT_HIMinBiasZDC_Scint_v1',
  'HLT_HIPhoton10_Photon15_v1',
  'HLT_HIPhoton15_Photon20_v1',
  'HLT_HIRandom_v1',
  'HLT_HISinglePhoton15_v1',
  'HLT_HISinglePhoton20_v1',
  'HLT_HISinglePhoton30_v1',
  'HLT_HISinglePhoton40_v1',
  'HLT_HIStoppedHSCP35_v1',
  'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
  'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
  'HLT_HIUpcEcal_v1',
  'HLT_HIUpcMu_v1',
  'HLT_HIZeroBiasPixel_SingleTrack_v1',
  'HLT_HIZeroBiasXOR_v1',
  'HLT_HIZeroBias_v1',
  'HLT_HcalCalibration_v2',
  'HLT_LogMonitor_v1',
  'HLT_TrackerCalibration_v2' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltDt4DSegments_*_*',
      'keep *_hltL1GtObjectMap_*_*',
      'keep FEDRawDataCollection_rawDataCollector_*_*',
      'keep FEDRawDataCollection_source_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEventWithRefs_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputEcalCalibration = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputEcalCalibration.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_EcalCalibration_v2' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltEcalCalibrationRaw_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputExpress = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputExpress.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_HIL2DoubleMu3_v1',
  'HLT_HIMinBiasHfOrBSC_v1' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltL1GtObjectMap_*_*',
      'keep FEDRawDataCollection_rawDataCollector_*_*',
      'keep FEDRawDataCollection_source_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputHLTDQM = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputHLTDQM.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
  'HLT_HIActivityHF_Single3_v1',
  'HLT_HIBptxXOR_v1',
  'HLT_HICentralityVeto_v1',
  'HLT_HIClusterVertexCompatibility_v1',
  'HLT_HIDoublePhoton10_v1',
  'HLT_HIDoublePhoton15_v1',
  'HLT_HIDoublePhoton20_v1',
  'HLT_HIFullTrack12_L1Central_v1',
  'HLT_HIFullTrack14_L1Central_v1',
  'HLT_HIFullTrack20_L1Central_v1',
  'HLT_HIFullTrack25_L1Central_v1',
  'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
  'HLT_HIL1DoubleMuOpen_v1',
  'HLT_HIL1SingleMu3_v1',
  'HLT_HIL1SingleMu5_v1',
  'HLT_HIL1SingleMu7_v1',
  'HLT_HIL2DoubleMu0_v1',
  'HLT_HIL2DoubleMu3_v1',
  'HLT_HIL2Mu20_v1',
  'HLT_HIL2Mu3_v1',
  'HLT_HIL2Mu5Tight_v1',
  'HLT_HIMinBiasBSC_OR_v1',
  'HLT_HIMinBiasBSC_v1',
  'HLT_HIMinBiasHF_v1',
  'HLT_HIMinBiasHfOrBSC_v1',
  'HLT_HIMinBiasHf_OR_v1',
  'HLT_HIMinBiasPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
  'HLT_HIMinBiasZDC_Calo_v1',
  'HLT_HIMinBiasZDC_Scint_v1',
  'HLT_HIPhoton10_Photon15_v1',
  'HLT_HIPhoton15_Photon20_v1',
  'HLT_HIRandom_v1',
  'HLT_HISinglePhoton15_v1',
  'HLT_HISinglePhoton20_v1',
  'HLT_HISinglePhoton30_v1',
  'HLT_HISinglePhoton40_v1',
  'HLT_HIStoppedHSCP35_v1',
  'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
  'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
  'HLT_HIUpcEcal_v1',
  'HLT_HIUpcMu_v1',
  'HLT_HIZeroBiasPixel_SingleTrack_v1',
  'HLT_HIZeroBiasXOR_v1',
  'HLT_HIZeroBias_v1',
  'HLT_LogMonitor_v1' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltAlCaEtaRecHitsFilter_*_*',
      'keep *_hltAlCaPhiSymStream_*_*',
      'keep *_hltAlCaPi0RecHitsFilter_*_*',
      'keep *_hltAntiKT5CaloJets_*_*',
      'keep *_hltAntiKT5PFJetsForTaus_*_*',
      'keep *_hltAntiKT5PFJets_*_*',
      'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonMu5L3_*_*',
      'keep *_hltCaloJetCorrectedRegional_*_*',
      'keep *_hltCaloJetCorrected_*_*',
      'keep *_hltConvPFTauTightIsoTrackPt20Isolation_*_*',
      'keep *_hltConvPFTauTightIsoTrackPt20_*_*',
      'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*',
      'keep *_hltConvPFTausTightIsoTrackFinding_*_*',
      'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*',
      'keep *_hltConvPFTausTightIsoTrackPt5_*_*',
      'keep *_hltConvPFTausTightIso_*_*',
      'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*',
      'keep *_hltConvPFTausTrackFinding_*_*',
      'keep *_hltConvPFTaus_*_*',
      'keep *_hltCsc2DRecHits_*_*',
      'keep *_hltCscSegments_*_*',
      'keep *_hltDoublePFTauTightIso25Track5_*_*',
      'keep *_hltDoublePFTauTightIso25Track_*_*',
      'keep *_hltDoublePFTauTightIso40Track5_*_*',
      'keep *_hltDoublePFTauTightIso40Track_*_*',
      'keep *_hltDoublePFTauTightIso45Track5_*_*',
      'keep *_hltDoublePFTauTightIso45Track_*_*',
      'keep *_hltDt4DSegments_*_*',
      'keep *_hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter_*_*',
      'keep *_hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter_*_*',
      'keep *_hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter_*_*',
      'keep *_hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter_*_*',
      'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*',
      'keep *_hltFilterDoubleIsoPFTau25Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterDoubleIsoPFTau40Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau25Trk5_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau40Trk5_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*',
      'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET45_*_*',
      'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*',
      'keep *_hltHITIPTCorrectorHB_*_*',
      'keep *_hltHITIPTCorrectorHE_*_*',
      'keep *_hltIsolPixelTrackProdHB_*_*',
      'keep *_hltIsolPixelTrackProdHE_*_*',
      'keep *_hltL1HLTDoubleIsoPFTau25Trk5JetsMatch_*_*',
      'keep *_hltL1HLTDoubleIsoPFTau40Trk5JetsMatch_*_*',
      'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*',
      'keep *_hltL1IsoElectronTrackIsol_*_*',
      'keep *_hltL1NonIsoElectronTrackIsol_*_*',
      'keep *_hltL1extraParticles_*_*',
      'keep *_hltL1sDoubleIsoTau25Trk5_*_*',
      'keep *_hltL1sDoubleIsoTau40Trk5_*_*',
      'keep *_hltL1sDoubleTauJet40Eta2p17orDoubleJet52Central_*_*',
      'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*',
      'keep *_hltL1sL1Jet52ETM30_*_*',
      'keep *_hltL1sL1SingleEG12_*_*',
      'keep *_hltL1sL1SingleEG15_*_*',
      'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*',
      'keep *_hltL1sL1SingleJet52ETM30_*_*',
      'keep *_hltL1sL1SingleMu10_*_*',
      'keep *_hltL1sL1SingleMu14Eta2p1_*_*',
      'keep *_hltL1sSingleIsoTau35Trk20MET45_*_*',
      'keep *_hltL2MuonCandidates_*_*',
      'keep *_hltL2MuonIsolations_*_*',
      'keep *_hltL2MuonSeeds_*_*',
      'keep *_hltL2Muons_*_*',
      'keep *_hltL2TauJets_*_*',
      'keep *_hltL3MuonCandidates_*_*',
      'keep *_hltL3MuonIsolations_*_*',
      'keep *_hltL3MuonsIOHit_*_*',
      'keep *_hltL3MuonsLinksCombination_*_*',
      'keep *_hltL3MuonsOIHit_*_*',
      'keep *_hltL3MuonsOIState_*_*',
      'keep *_hltL3Muons_*_*',
      'keep *_hltL3TkFromL2OICombination_*_*',
      'keep *_hltL3TkTracksFromL2IOHit_*_*',
      'keep *_hltL3TkTracksFromL2OIHit_*_*',
      'keep *_hltL3TkTracksFromL2OIState_*_*',
      'keep *_hltL3TkTracksFromL2_*_*',
      'keep *_hltL3TrackCandidateFromL2IOHit_*_*',
      'keep *_hltL3TrackCandidateFromL2OIHit_*_*',
      'keep *_hltL3TrackCandidateFromL2OIState_*_*',
      'keep *_hltL3TrajSeedIOHit_*_*',
      'keep *_hltL3TrajSeedOIHit_*_*',
      'keep *_hltL3TrajSeedOIState_*_*',
      'keep *_hltL3TrajectorySeed_*_*',
      'keep *_hltMet_*_*',
      'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*',
      'keep *_hltMuTrackJpsiCtfTrackCands_*_*',
      'keep *_hltMuTrackJpsiPixelTrackCands_*_*',
      'keep *_hltMuonCSCDigis_*_*',
      'keep *_hltMuonDTDigis_*_*',
      'keep *_hltMuonRPCDigis_*_*',
      'keep *_hltOfflineBeamSpot_*_*',
      'keep *_hltOnlineBeamSpot_*_*',
      'keep *_hltOverlapFilterIsoEle15IsoPFTau15_*_*',
      'keep *_hltOverlapFilterIsoEle18LooseIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoEle18TightIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15IsoPFTau15_*_*',
      'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*',
      'keep *_hltPFTau15TrackLooseIso_*_*',
      'keep *_hltPFTau15Track_*_*',
      'keep *_hltPFTau15_*_*',
      'keep *_hltPFTau20TrackLooseIso_*_*',
      'keep *_hltPFTau20Track_*_*',
      'keep *_hltPFTau20_*_*',
      'keep *_hltPFTauJetTracksAssociator_*_*',
      'keep *_hltPFTauMediumIso20TrackMediumIso_*_*',
      'keep *_hltPFTauMediumIso20Track_*_*',
      'keep *_hltPFTauMediumIso20_*_*',
      'keep *_hltPFTauMediumIso35Track_*_*',
      'keep *_hltPFTauMediumIso35_*_*',
      'keep *_hltPFTauTagInfo_*_*',
      'keep *_hltPFTauTightIso20TrackTightIso_*_*',
      'keep *_hltPFTauTightIso20Track_*_*',
      'keep *_hltPFTauTightIso20_*_*',
      'keep *_hltPFTauTightIso35Track_*_*',
      'keep *_hltPFTauTightIso35_*_*',
      'keep *_hltParticleFlowForTaus_*_*',
      'keep *_hltParticleFlow_*_*',
      'keep *_hltPixelMatchElectronsL1Iso_*_*',
      'keep *_hltPixelMatchElectronsL1NonIso_*_*',
      'keep *_hltRpcRecHits_*_*',
      'keep *_hltSiStripRawToClustersFacility_*_*',
      'keep *_hltSingleMu15L3Filtered15_*_*',
      'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*',
      'keep *_hltSingleMuIsoL3IsoFiltered15_*_*',
      'keep *_hltTowerMakerForMuons_*_*',
      'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*',
      'keep TrackingRecHitsOwned_hltL3Muons_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep recoCaloMETs_hltMet_*_*',
      'keep recoElectronSeeds_hltL1IsoStartUpElectronPixelSeeds_*_*',
      'keep recoElectronSeeds_hltL1NonIsoStartUpElectronPixelSeeds_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*',
      'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*',
      'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*',
      'keep triggerTriggerEventWithRefs_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputHLTDQMResults = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputHLTDQMResults.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLTriggerFinalPath' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep L1GlobalTriggerReadoutRecord_hltGtDigis_*_*',
      'keep LumiScalerss_hltScalersRawToDigi_*_*',
      'keep edmTriggerResults_*_*_*' )
)
process.hltOutputHLTMON = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputHLTMON.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_HIActivityHF_Coincidence3_v1',
  'HLT_HIActivityHF_Single3_v1',
  'HLT_HIBptxXOR_v1',
  'HLT_HICentralityVeto_v1',
  'HLT_HIClusterVertexCompatibility_v1',
  'HLT_HIDoublePhoton10_v1',
  'HLT_HIDoublePhoton15_v1',
  'HLT_HIDoublePhoton20_v1',
  'HLT_HIFullTrack12_L1Central_v1',
  'HLT_HIFullTrack14_L1Central_v1',
  'HLT_HIFullTrack20_L1Central_v1',
  'HLT_HIFullTrack25_L1Central_v1',
  'HLT_HIL1Algo_BptxXOR_BSC_OR_v1',
  'HLT_HIL1DoubleMuOpen_v1',
  'HLT_HIL1SingleMu3_v1',
  'HLT_HIL1SingleMu5_v1',
  'HLT_HIL1SingleMu7_v1',
  'HLT_HIL2DoubleMu0_v1',
  'HLT_HIL2DoubleMu3_v1',
  'HLT_HIL2Mu20_v1',
  'HLT_HIL2Mu3_v1',
  'HLT_HIL2Mu5Tight_v1',
  'HLT_HIMinBiasBSC_OR_v1',
  'HLT_HIMinBiasBSC_v1',
  'HLT_HIMinBiasHF_v1',
  'HLT_HIMinBiasHfOrBSC_v1',
  'HLT_HIMinBiasHf_OR_v1',
  'HLT_HIMinBiasPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDCPixel_SingleTrack_v1',
  'HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1',
  'HLT_HIMinBiasZDC_Calo_v1',
  'HLT_HIMinBiasZDC_Scint_v1',
  'HLT_HIPhoton10_Photon15_v1',
  'HLT_HIPhoton15_Photon20_v1',
  'HLT_HIRandom_v1',
  'HLT_HISinglePhoton15_v1',
  'HLT_HISinglePhoton20_v1',
  'HLT_HISinglePhoton30_v1',
  'HLT_HISinglePhoton40_v1',
  'HLT_HIStoppedHSCP35_v1',
  'HLT_HIUPCNeuEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1',
  'HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1',
  'HLT_HIUPCNeuMuPixel_SingleTrack_v1',
  'HLT_HIUpcEcal_v1',
  'HLT_HIUpcMu_v1',
  'HLT_HIZeroBiasPixel_SingleTrack_v1',
  'HLT_HIZeroBiasXOR_v1',
  'HLT_HIZeroBias_v1',
  'HLT_HcalCalibration_v2',
  'HLT_LogMonitor_v1' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltAlCaEtaRecHitsFilter_*_*',
      'keep *_hltAlCaPi0RecHitsFilter_*_*',
      'keep *_hltAntiKT5CaloJets_*_*',
      'keep *_hltAntiKT5PFJetsForTaus_*_*',
      'keep *_hltAntiKT5PFJets_*_*',
      'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*',
      'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*',
      'keep *_hltBSoftMuonMu5L3_*_*',
      'keep *_hltCaloJetCorrectedRegional_*_*',
      'keep *_hltCaloJetCorrected_*_*',
      'keep *_hltConvPFTauTightIsoTrackPt20Isolation_*_*',
      'keep *_hltConvPFTauTightIsoTrackPt20_*_*',
      'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*',
      'keep *_hltConvPFTausTightIsoTrackFinding_*_*',
      'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*',
      'keep *_hltConvPFTausTightIsoTrackPt5_*_*',
      'keep *_hltConvPFTausTightIso_*_*',
      'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*',
      'keep *_hltConvPFTausTrackFinding_*_*',
      'keep *_hltConvPFTaus_*_*',
      'keep *_hltCscSegments_*_*',
      'keep *_hltDoublePFTauTightIso25Track5_*_*',
      'keep *_hltDoublePFTauTightIso25Track_*_*',
      'keep *_hltDoublePFTauTightIso40Track5_*_*',
      'keep *_hltDoublePFTauTightIso40Track_*_*',
      'keep *_hltDoublePFTauTightIso45Track5_*_*',
      'keep *_hltDoublePFTauTightIso45Track_*_*',
      'keep *_hltDt4DSegments_*_*',
      'keep *_hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTDphiFilter_*_*',
      'keep *_hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter_*_*',
      'keep *_hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter_*_*',
      'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*',
      'keep *_hltFilterDoubleIsoPFTau25Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterDoubleIsoPFTau40Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau25Trk5_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau40Trk5_*_*',
      'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*',
      'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET45_*_*',
      'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*',
      'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*',
      'keep *_hltL1HLTDoubleIsoPFTau25Trk5JetsMatch_*_*',
      'keep *_hltL1HLTDoubleIsoPFTau40Trk5JetsMatch_*_*',
      'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*',
      'keep *_hltL1IsoElectronTrackIsol_*_*',
      'keep *_hltL1NonIsoElectronTrackIsol_*_*',
      'keep *_hltL1extraParticles_*_*',
      'keep *_hltL1sDoubleIsoTau25Trk5_*_*',
      'keep *_hltL1sDoubleIsoTau40Trk5_*_*',
      'keep *_hltL1sDoubleTauJet40Eta2p17orDoubleJet52Central_*_*',
      'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*',
      'keep *_hltL1sL1Jet52ETM30_*_*',
      'keep *_hltL1sL1SingleEG12_*_*',
      'keep *_hltL1sL1SingleEG15_*_*',
      'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*',
      'keep *_hltL1sL1SingleJet52ETM30_*_*',
      'keep *_hltL1sL1SingleMu10_*_*',
      'keep *_hltL1sL1SingleMu14Eta2p1_*_*',
      'keep *_hltL1sSingleIsoTau35Trk20MET45_*_*',
      'keep *_hltL2MuonCandidatesNoVtx_*_*',
      'keep *_hltL2Muons_*_*',
      'keep *_hltL2TauJets_*_*',
      'keep *_hltL3MuonCandidatesNoVtx_*_*',
      'keep *_hltL3Muons_*_*',
      'keep *_hltL3TkTracksFromL2_*_*',
      'keep *_hltMet_*_*',
      'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*',
      'keep *_hltMuonCSCDigis_*_*',
      'keep *_hltMuonDTDigis_*_*',
      'keep *_hltMuonRPCDigis_*_*',
      'keep *_hltOnlineBeamSpot_*_*',
      'keep *_hltOverlapFilterIsoEle15IsoPFTau15_*_*',
      'keep *_hltOverlapFilterIsoEle18LooseIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoEle18TightIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15IsoPFTau15_*_*',
      'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*',
      'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*',
      'keep *_hltPFTau15TrackLooseIso_*_*',
      'keep *_hltPFTau15Track_*_*',
      'keep *_hltPFTau15_*_*',
      'keep *_hltPFTau20TrackLooseIso_*_*',
      'keep *_hltPFTau20Track_*_*',
      'keep *_hltPFTau20_*_*',
      'keep *_hltPFTauJetTracksAssociator_*_*',
      'keep *_hltPFTauMediumIso20TrackMediumIso_*_*',
      'keep *_hltPFTauMediumIso20Track_*_*',
      'keep *_hltPFTauMediumIso20_*_*',
      'keep *_hltPFTauMediumIso35Track_*_*',
      'keep *_hltPFTauMediumIso35_*_*',
      'keep *_hltPFTauTagInfo_*_*',
      'keep *_hltPFTauTightIso20TrackTightIso_*_*',
      'keep *_hltPFTauTightIso20Track_*_*',
      'keep *_hltPFTauTightIso20_*_*',
      'keep *_hltPFTauTightIso35Track_*_*',
      'keep *_hltPFTauTightIso35_*_*',
      'keep *_hltParticleFlowForTaus_*_*',
      'keep *_hltParticleFlow_*_*',
      'keep *_hltPixelMatchElectronsL1Iso_*_*',
      'keep *_hltPixelMatchElectronsL1NonIso_*_*',
      'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*',
      'keep *_hltRpcRecHits_*_*',
      'keep *_hltSingleMu15L3Filtered15_*_*',
      'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*',
      'keep *_hltSingleMuIsoL3IsoFiltered15_*_*',
      'keep FEDRawDataCollection_rawDataCollector_*_*',
      'keep FEDRawDataCollection_source_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*',
      'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*',
      'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*',
      'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*',
      'keep triggerTriggerEventWithRefs_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)
process.hltOutputTrackerCalibration = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputTrackerCalibration.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'HLT_TrackerCalibration_v2' ) ),
    outputCommands = cms.untracked.vstring( 'drop *_hlt*_*_*',
      'keep *_hltTrackerCalibrationRaw_*_*',
      'keep edmTriggerResults_*_*_*',
      'keep triggerTriggerEvent_*_*_*' )
)

process.HLTBeginSequenceCalibration = cms.Sequence( process.hltCalibrationEventsFilter + process.hltGtDigis )
process.HLTEndSequence = cms.Sequence( process.hltBoolEnd )
process.HLTL1UnpackerSequence = cms.Sequence( process.hltGtDigis + process.hltGctDigis + process.hltL1GtObjectMap + process.hltL1extraParticles )
process.HLTBeamSpot = cms.Sequence( process.hltScalersRawToDigi + process.hltOnlineBeamSpot + process.hltOfflineBeamSpot )
process.HLTBeginSequence = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence + process.HLTBeamSpot )
process.HLTDoHILocalPixelSequence = cms.Sequence( process.hltSiPixelDigis + process.hltHISiPixelClusters + process.hltHISiPixelRecHits )
process.HLTPixelTrackingForHITrackTrigger = cms.Sequence( process.hltHIPixelClusterVertices + process.hltPixelTracksForHITrackTrigger + process.hltPixelCandsForHITrackTrigger )
process.HLTMuonLocalRecoSequence = cms.Sequence( process.hltMuonDTDigis + process.hltDt1DRecHits + process.hltDt4DSegments + process.hltMuonCSCDigis + process.hltCsc2DRecHits + process.hltCscSegments + process.hltMuonRPCDigis + process.hltRpcRecHits )
process.HLTL2muonrecoNocandSequence = cms.Sequence( process.HLTMuonLocalRecoSequence + process.hltL2MuonSeeds + process.hltL2Muons )
process.HLTL2muonrecoSequence = cms.Sequence( process.HLTL2muonrecoNocandSequence + process.hltL2MuonCandidates )
process.HLTDoLocalHcalSequence = cms.Sequence( process.hltHcalDigis + process.hltHbhereco + process.hltHfreco + process.hltHoreco )
process.HLTDoCaloSequence = cms.Sequence( process.hltEcalRawToRecHitFacility + process.hltEcalRegionalRestFEDs + process.hltEcalRecHitAll + process.HLTDoLocalHcalSequence + process.hltTowerMakerForAll )
process.HLTDoHIEcalClusWithCleaningSequence = cms.Sequence( process.hltIslandBasicClustersHI + process.hltIslandSuperClustersHI + process.hltCorrectedIslandBarrelSuperClustersHI + process.hltCorrectedIslandEndcapSuperClustersHI + process.hltCleanedCorrectedIslandBarrelSuperClustersHI + process.hltRecoHIEcalWithCleaningCandidate )
process.HLTBeginSequenceAntiBPTX = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence + process.hltBPTXAntiCoincidence + process.HLTBeamSpot )
process.HLTBeginSequenceBPTX = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence + process.hltBPTXCoincidence + process.HLTBeamSpot )
process.HLTPixelSeedingForHITrackTrigger = cms.Sequence( process.hltHIPixelClusterVerticesForHITrackTrigger + process.hltHIPixel3ProtoTracks + process.hltHIPixelMedianVertex + process.hltHISelectedProtoTracks + process.hltHIPixelAdaptiveVertex + process.hltHIBestAdaptiveVertex + process.hltHISelectedVertex + process.hltHIPixel3PrimTracks + process.hltHIPixelTrackCandsForHITrackTrigger )
process.HLTDoHILocalStripSequenceNonRegional = cms.Sequence( process.hltSiStripRawToDigi + process.hltSiStripZeroSuppression + process.hltHISiStripClustersNonRegional )
process.HLTFullTrackingForHITrackTrigger = cms.Sequence( process.hltHIPixelTrackSeeds + process.hltHIPrimTrackCandidates + process.hltHIGlobalPrimTracks + process.hltHIGoodLooseTracks + process.hltHIFullTrackCandsForHITrackTrigger )

process.HLTriggerFirstPath = cms.Path( process.hltGetRaw + process.hltBoolFalse )
process.HLT_DTCalibration_v1 = cms.Path( process.HLTBeginSequenceCalibration + process.hltPreDTCalibration + process.hltDTCalibrationRaw + process.HLTEndSequence )
process.HLT_EcalCalibration_v2 = cms.Path( process.HLTBeginSequenceCalibration + process.hltPreEcalCalibration + process.hltEcalCalibrationRaw + process.HLTEndSequence )
process.HLT_HcalCalibration_v2 = cms.Path( process.HLTBeginSequenceCalibration + process.hltPreHcalCalibration + process.hltHcalCalibTypeFilter + process.hltHcalCalibrationRaw + process.HLTEndSequence )
process.HLT_TrackerCalibration_v2 = cms.Path( process.HLTBeginSequenceCalibration + process.hltPreTrackerCalibration + process.hltLaserAlignmentEventFilter + process.hltTrackerCalibrationRaw + process.HLTEndSequence )
process.HLT_LogMonitor_v1 = cms.Path( process.hltGtDigis + process.hltPreLogMonitor + process.hltLogMonitorFilter + process.HLTEndSequence )
process.HLT_HIZeroBias_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIZeroBias + process.hltPreHIZeroBias + process.HLTEndSequence )
process.HLT_HIZeroBiasXOR_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1BptxXOR + process.hltPreHIZeroBiasXOR + process.HLTEndSequence )
process.HLT_HIZeroBiasPixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIZeroBiasXOR + process.hltPreHIZeroBiasPixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIMinBiasBSC_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasBSC + process.hltPreHIMinBiasBSC + process.HLTEndSequence )
process.HLT_HIMinBiasBSC_OR_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasBSCOR + process.hltPreHIMinBiasBSCOR + process.HLTEndSequence )
process.HLT_HIMinBiasHF_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasHF + process.hltPreHIMinBiasHF + process.HLTEndSequence )
process.HLT_HIMinBiasHf_OR_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasHfOr + process.hltPreHIMinBiasHfOR + process.HLTEndSequence )
process.HLT_HIMinBiasHfOrBSC_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasHfOrBSC + process.hltPreHIMinBiasHfOrBSC + process.HLTEndSequence )
process.HLT_HIMinBiasPixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasHfOrBSC + process.hltPreHIMinBiasPixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIMinBiasZDC_Calo_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasZDC + process.hltPreHIMinBiasZDCCalo + process.HLTEndSequence )
process.HLT_HIMinBiasZDC_Calo_PlusOrMinus_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasZDCCaloPlusOrMinus + process.hltPreHIMinBiasZDCCaloPlusOrMinus + process.HLTEndSequence )
process.HLT_HIMinBiasZDC_Scint_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasZDCScint + process.hltPreHIMinBiasZDCScint + process.HLTEndSequence )
process.HLT_HIMinBiasZDCPixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasZDCPixelSingleTrack + process.hltPreHIMinBiasZDCPixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIBptxXOR_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1BptxXOR + process.hltPreHIBptxXOR + process.HLTEndSequence )
process.HLT_HIL1Algo_BptxXOR_BSC_OR_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1BptxXORBscMinBiasOR + process.hltPreHIL1AlgoBptxXORBSCOR + process.HLTEndSequence )
process.HLT_HIL1SingleMu3_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu3BptxAND + process.hltPreHIL1SingleMu3 + process.HLTEndSequence )
process.HLT_HIL1SingleMu5_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu5 + process.hltPreHIL1SingleMu5 + process.HLTEndSequence )
process.HLT_HIL1SingleMu7_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu7 + process.hltPreHIL1SingleMu7 + process.HLTEndSequence )
process.HLT_HIL1DoubleMuOpen_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1DoubleMuOpenBptxAND + process.hltPreHIL1DoubleMuOpen + process.hltHIDoubleMuLevel1PathL1OpenFiltered + process.HLTEndSequence )
process.HLT_HIL2Mu3_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu3BptxANDwithBSCHF + process.hltPreHIL2Mu3 + process.hltHIL1SingleMu3withBSCHFL1Filtered + process.HLTL2muonrecoSequence + process.hltHIL2Mu3L2Filtered + process.HLTEndSequence )
process.HLT_HIL2Mu5Tight_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu3BptxANDwithBSCHF + process.hltPreHIL2Mu5Tight + process.hltHIL1SingleMu3withBSCHFL1Filtered + process.HLTL2muonrecoSequence + process.hltHIL2Mu5TightL2Filtered + process.HLTEndSequence )
process.HLT_HIL2Mu20_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleMu3BptxANDwithBSCHF + process.hltPreHIL2Mu20 + process.hltHIL1SingleMu3withBSCHFL1Filtered + process.HLTL2muonrecoSequence + process.hltHIL2Mu20L2Filtered + process.HLTEndSequence )
process.HLT_HIL2DoubleMu0_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1DoubleMuOpenBptxANDwithBSCHF + process.hltPreHIL2DoubleMu0 + process.hltHIDoubleMuLevel1PathL1OpenWithBSCHFFiltered + process.HLTL2muonrecoSequence + process.hltHIL2DoubleMu0L2Filtered + process.HLTEndSequence )
process.HLT_HIL2DoubleMu3_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1DoubleMuOpenBptxANDwithBSCHF + process.hltPreHIL2DoubleMu3 + process.hltHIDoubleMuLevel1PathL1OpenWithBSCHFFiltered + process.HLTL2muonrecoSequence + process.hltHIL2DoubleMu3L2Filtered + process.HLTEndSequence )
process.HLT_HIUpcEcal_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUpcEcal + process.hltPreHIUpcEcal + process.HLTEndSequence )
process.HLT_HIUpcMu_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUpcMu + process.hltPreHIUpcMu + process.HLTEndSequence )
process.HLT_HISinglePhoton15_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIPhoton15 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIPhoton15 + process.HLTEndSequence )
process.HLT_HISinglePhoton20_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIPhoton20 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIPhoton20 + process.HLTEndSequence )
process.HLT_HISinglePhoton30_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIPhoton30 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIPhoton30 + process.HLTEndSequence )
process.HLT_HISinglePhoton40_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIPhoton40 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIPhoton40 + process.HLTEndSequence )
process.HLT_HIPhoton10_Photon15_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIDoublePhoton10and15 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIDoublePhoton1015Filter1 + process.hltHIDoublePhoton1015Filter2 + process.HLTEndSequence )
process.HLT_HIPhoton15_Photon20_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIDoublePhoton15and20 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIDoublePhoton1520Filter1 + process.hltHIDoublePhoton1520Filter2 + process.HLTEndSequence )
process.HLT_HIDoublePhoton10_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIDoublePhoton10 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIDoublePhoton10 + process.HLTEndSequence )
process.HLT_HIDoublePhoton15_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIDoublePhoton15 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIDoublePhoton15 + process.HLTEndSequence )
process.HLT_HIDoublePhoton20_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1SingleEG5BptxAND + process.hltPreHIDoublePhoton20 + process.HLTDoCaloSequence + process.HLTDoHIEcalClusWithCleaningSequence + process.hltHIDoublePhoton20 + process.HLTEndSequence )
process.HLT_HIStoppedHSCP35_v1 = cms.Path( process.HLTBeginSequenceAntiBPTX + process.hltL1sL1SingleJet20UNotBptxOR + process.hltPreHIStoppedHSCP35 + process.hltHcalDigis + process.hltHbhereco + process.hltStoppedHSCPHpdFilter + process.hltStoppedHSCPTowerMakerForAll + process.hltStoppedHSCPIterativeCone5CaloJets + process.hltStoppedHSCP1CaloJetEnergy35 + process.HLTEndSequence )
process.HLT_HIActivityHF_Coincidence3_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1GlobalDecision + process.hltPreHIActivityHFCoincidence3 + process.hltHcalDigis + process.hltHfreco + process.hltHcalSimpleRecHitFilterCoincidence + process.HLTEndSequence )
process.HLT_HIActivityHF_Single3_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1GlobalDecision + process.hltPreHIActivityHFSingle3 + process.hltHcalDigis + process.hltHfreco + process.hltHcalSimpleRecHitFilter + process.HLTEndSequence )
process.HLT_HIClusterVertexCompatibility_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sL1GlobalDecision + process.hltPreHIClusterVertexCompatibility + process.HLTDoHILocalPixelSequence + process.hltHIPixelClusterShapeFilter + process.HLTEndSequence )
process.HLT_HICentralityVeto_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIMinBiasHfOrBSC + process.hltPreHICentralityVeto + process.HLTDoHILocalPixelSequence + process.hltPixelActivityFilterCentralityVeto + process.HLTEndSequence )
process.HLT_HIFullTrack12_L1Central_v1 = cms.Path( process.HLTBeginSequenceBPTX + process.hltL1sETT100 + process.hltPreHIFullTrack12L1Central + process.HLTDoCaloSequence + process.hltHICaloTowerFilter4 + process.HLTDoHILocalPixelSequence + process.HLTPixelSeedingForHITrackTrigger + process.hltHISinglePixelTrackFilter10 + process.HLTDoHILocalStripSequenceNonRegional + process.HLTFullTrackingForHITrackTrigger + process.hltHISingleFullTrackFilter12 + process.HLTEndSequence )
process.HLT_HIFullTrack14_L1Central_v1 = cms.Path( process.HLTBeginSequenceBPTX + process.hltL1sETT100 + process.hltPreHIFullTrack14L1Central + process.HLTDoCaloSequence + process.hltHICaloTowerFilter4 + process.HLTDoHILocalPixelSequence + process.HLTPixelSeedingForHITrackTrigger + process.hltHISinglePixelTrackFilter10 + process.HLTDoHILocalStripSequenceNonRegional + process.HLTFullTrackingForHITrackTrigger + process.hltHISingleFullTrackFilter14 + process.HLTEndSequence )
process.HLT_HIFullTrack20_L1Central_v1 = cms.Path( process.HLTBeginSequenceBPTX + process.hltL1sETT100 + process.hltPreHIFullTrack20L1Central + process.HLTDoCaloSequence + process.hltHICaloTowerFilter4 + process.HLTDoHILocalPixelSequence + process.HLTPixelSeedingForHITrackTrigger + process.hltHISinglePixelTrackFilter10 + process.HLTDoHILocalStripSequenceNonRegional + process.HLTFullTrackingForHITrackTrigger + process.hltHISingleFullTrackFilter20 + process.HLTEndSequence )
process.HLT_HIFullTrack25_L1Central_v1 = cms.Path( process.HLTBeginSequenceBPTX + process.hltL1sETT100 + process.hltPreHIFullTrack25L1Central + process.HLTDoCaloSequence + process.hltHICaloTowerFilter4 + process.HLTDoHILocalPixelSequence + process.HLTPixelSeedingForHITrackTrigger + process.hltHISinglePixelTrackFilter10 + process.HLTDoHILocalStripSequenceNonRegional + process.HLTFullTrackingForHITrackTrigger + process.hltHISingleFullTrackFilter25 + process.HLTEndSequence )
process.HLT_HIRandom_v1 = cms.Path( process.hltRandomEventsFilter + process.HLTL1UnpackerSequence + process.hltPreHIRandom + process.HLTEndSequence )
process.HLT_HIUPCNeuMuPixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuMu + process.hltPreHIUPCNeuMuPixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIUPCNeuEG2Pixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuEG2 + process.hltPreHIUPCNeuEG2PixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIUPCNeuEG5Pixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuEG5 + process.hltPreHIUPCNeuEG5PixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIUPCNeuHcalHfMuPixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuHcalHfMu + process.hltPreHIUPCNeuHcalHfMuPixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIUPCNeuHcalHfEG2Pixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuHcalHfEG2 + process.hltPreHIUPCNeuHcalHfEG2PixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLT_HIUPCNeuHcalHfEG5Pixel_SingleTrack_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sHIUPCNeuHcalHfEG5 + process.hltPreHIUPCNeuHcalHfEG5PixelSingleTrack + process.HLTDoHILocalPixelSequence + process.HLTPixelTrackingForHITrackTrigger + process.hltHISinglePixelTrackFilter + process.HLTEndSequence )
process.HLTriggerFinalPath = cms.Path( process.hltGtDigis + process.hltScalersRawToDigi + process.hltFEDSelector + process.hltTriggerSummaryAOD + process.hltTriggerSummaryRAW )
process.HLTAnalyzerEndpath = cms.EndPath( process.hltL1GtTrigReport + process.hltTrigReport )
process.AOutput = cms.EndPath( process.hltOutputA )
process.ALCAP0Output = cms.EndPath( process.hltPreALCAP0Output )
process.ALCAPHISYMOutput = cms.EndPath( process.hltPreALCAPHISYMOutput )
process.CalibrationOutput = cms.EndPath( process.hltPreCalibrationOutput + process.hltOutputCalibration )
process.DQMOutput = cms.EndPath( process.hltDQML1Scalers + process.hltDQML1SeedLogicScalers + process.hltDQMHLTScalers + process.hltPreDQMOutput + process.hltPreDQMOutputSmart + process.hltOutputDQM )
process.EcalCalibrationOutput = cms.EndPath( process.hltPreEcalCalibrationOutput + process.hltOutputEcalCalibration )
process.ExpressOutput = cms.EndPath( process.hltPreExpressOutput + process.hltPreExpressOutputSmart + process.hltOutputExpress )
process.HLTDQMOutput = cms.EndPath( process.hltPreHLTDQMOutput + process.hltPreHLTDQMOutputSmart + process.hltOutputHLTDQM )
process.HLTDQMResultsOutput = cms.EndPath( process.hltPreHLTDQMResultsOutput + process.hltPreHLTDQMResultsOutputSmart + process.hltOutputHLTDQMResults )
process.HLTMONOutput = cms.EndPath( process.hltPreHLTMONOutput + process.hltPreHLTMONOutputSmart + process.hltOutputHLTMON )
process.NanoDSTOutput = cms.EndPath( process.hltPreNanoDSTOutput )
process.RPCMONOutput = cms.EndPath( process.hltPreRPCMONOutput )
process.TrackerCalibrationOutput = cms.EndPath( process.hltPreTrackerCalibrationOutput + process.hltOutputTrackerCalibration )


# override the process name
process.setName_('HLTHIon')

# adapt HLT modules to the correct process name
if 'hltTrigReport' in process.__dict__:
    process.hltTrigReport.HLTriggerResults       = cms.InputTag( 'TriggerResults', '', 'HLTHIon' )

if 'hltPreExpressSmart' in process.__dict__:
    process.hltPreExpressSmart.TriggerResultsTag = cms.InputTag( 'TriggerResults', '', 'HLTHIon' )

if 'hltPreHLTMONSmart' in process.__dict__:
    process.hltPreHLTMONSmart.TriggerResultsTag  = cms.InputTag( 'TriggerResults', '', 'HLTHIon' )

if 'hltPreDQMSmart' in process.__dict__:
    process.hltPreDQMSmart.TriggerResultsTag     = cms.InputTag( 'TriggerResults', '', 'HLTHIon' )

if 'hltDQMHLTScalers' in process.__dict__:
    process.hltDQMHLTScalers.triggerResults      = cms.InputTag( 'TriggerResults', '', 'HLTHIon' )
    process.hltDQMHLTScalers.processname         = 'HLTHIon'

if 'hltDQML1SeedLogicScalers' in process.__dict__:
    process.hltDQML1SeedLogicScalers.processname = 'HLTHIon'

# remove the HLT prescales
if 'PrescaleService' in process.__dict__:
    process.PrescaleService.lvl1DefaultLabel = cms.untracked.string( '0' )
    process.PrescaleService.lvl1Labels       = cms.vstring( '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
    process.PrescaleService.prescaleTable    = cms.VPSet( )

# Disable HF Noise filters in HIon menu
if 'hltHfreco' in process.__dict__:
    process.hltHfreco.setNoiseFlags = cms.bool( False )

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

# enable the TrigReport and TimeReport
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True )
)

# override the GlobalTag, connection string and pfnPrefix
if 'GlobalTag' in process.__dict__:
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = autoCond['hltonline']

# override the L1 menu
if 'GlobalTag' in process.__dict__:
    process.GlobalTag.toGet.append(
        cms.PSet(
            record  = cms.string( 'L1GtTriggerMenuRcd' ),
            tag     = cms.string( 'L1GtTriggerMenu_L1Menu_CollisionsHeavyIons2010_v2_mc' ),
            label   = cms.untracked.string( '' ),
            connect = cms.untracked.string( 'frontier://FrontierProd/CMS_COND_31X_L1T' )
        )
    )

if 'MessageLogger' in process.__dict__:
    process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
    process.MessageLogger.categories.append('L1GtTrigReport')
    process.MessageLogger.categories.append('HLTrigReport')

# version specific customizations
import os
cmsswVersion = os.environ['CMSSW_VERSION']

# from CMSSW_4_4_0 / CMSSW_5_0_0_pre1: change input label for V00-04-17-00 RecoVertex/BeamSpotProducer
if cmsswVersion > "CMSSW_4_4":
    if 'hltOnlineBeamSpot' in process.__dict__:
        process.hltOnlineBeamSpot.src = process.hltOnlineBeamSpot.label
        del process.hltOnlineBeamSpot.label

# from CMSSW_4_4_0_pre8: update HF configuration for V00-09-18 RecoLocalCalo/HcalRecProducers
if cmsswVersion > "CMSSW_4_4":
    if 'hltHfreco' in process.__dict__:
        process.hltHfreco.digiTimeFromDB = cms.bool( False )
        process.hltHfreco.digistat.HFdigiflagCoef = cms.vdouble(
            process.hltHfreco.digistat.HFdigiflagCoef0.value(),
            process.hltHfreco.digistat.HFdigiflagCoef1.value(),
            process.hltHfreco.digistat.HFdigiflagCoef2.value()
        )
        del process.hltHfreco.digistat.HFdigiflagCoef0
        del process.hltHfreco.digistat.HFdigiflagCoef1
        del process.hltHfreco.digistat.HFdigiflagCoef2

# from CMSSW_4_4_0_pre6: updated configuration for the HybridClusterProducer's and EgammaHLTHybridClusterProducer's
if cmsswVersion > "CMSSW_4_4":
    if 'hltHybridSuperClustersActivity' in process.__dict__:
        process.hltHybridSuperClustersActivity.xi               = cms.double( 0.0 )
        process.hltHybridSuperClustersActivity.useEtForXi       = cms.bool( False )
    if 'hltHybridSuperClustersL1Isolated' in process.__dict__:
        process.hltHybridSuperClustersL1Isolated.xi             = cms.double( 0.0 )
        process.hltHybridSuperClustersL1Isolated.useEtForXi     = cms.bool( False )
    if 'hltHybridSuperClustersL1NonIsolated' in process.__dict__:
        process.hltHybridSuperClustersL1NonIsolated.xi          = cms.double( 0.0 )
        process.hltHybridSuperClustersL1NonIsolated.useEtForXi  = cms.bool( False )

# from CMSSW_4_4_0_pre5: updated configuration for the PFRecoTauDiscriminationByIsolation producers
if cmsswVersion > "CMSSW_4_4":
    if 'hltPFTauTightIsoIsolationDiscriminator' in process.__dict__:
        process.hltPFTauTightIsoIsolationDiscriminator.qualityCuts.primaryVertexSrc = process.hltPFTauTightIsoIsolationDiscriminator.PVProducer
        process.hltPFTauTightIsoIsolationDiscriminator.qualityCuts.pvFindingAlgo    = cms.string('highestPtInEvent')
        del process.hltPFTauTightIsoIsolationDiscriminator.PVProducer
    if 'hltPFTauLooseIsolationDiscriminator' in process.__dict__:
        process.hltPFTauLooseIsolationDiscriminator.qualityCuts.primaryVertexSrc = process.hltPFTauLooseIsolationDiscriminator.PVProducer
        process.hltPFTauLooseIsolationDiscriminator.qualityCuts.pvFindingAlgo    = cms.string('highestPtInEvent')
        del process.hltPFTauLooseIsolationDiscriminator.PVProducer

# from CMSSW_4_4_0_pre5: updated configuration for the EcalSeverityLevelESProducer
if cmsswVersion > "CMSSW_4_4":
    process.ecalSeverityLevel = cms.ESProducer("EcalSeverityLevelESProducer",
        appendToDataLabel = cms.string(''),
        dbstatusMask=cms.PSet(
            kGood        = cms.vuint32(0),
            kProblematic = cms.vuint32(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
            kRecovered   = cms.vuint32(),
            kTime        = cms.vuint32(),
            kWeird       = cms.vuint32(),
            kBad         = cms.vuint32(11, 12, 13, 14, 15, 16)
        ),
        flagMask = cms.PSet (
            kGood        = cms.vstring('kGood'),
            kProblematic = cms.vstring('kPoorReco', 'kPoorCalib', 'kNoisy', 'kSaturated'),
            kRecovered   = cms.vstring('kLeadingEdgeRecovered', 'kTowerRecovered'),
            kTime        = cms.vstring('kOutOfTime'),
            kWeird       = cms.vstring('kWeird', 'kDiWeird'),
            kBad         = cms.vstring('kFaultyHardware', 'kDead', 'kKilled')
        ),
        timeThresh = cms.double(2.0)
    )

# from CMSSW_4_4_0_pre3: additional ESProducer in cfg files
if cmsswVersion > "CMSSW_4_4":
    process.hltSiPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
        ListOfRecordToMerge = cms.VPSet(
            cms.PSet( record = cms.string("SiPixelQualityFromDbRcd"),
                      tag    = cms.string("")
                    ),
            cms.PSet( record = cms.string("SiPixelDetVOffRcd"),
                      tag    = cms.string("")
                    )
        )
    )

# from CMSSW_4_3_0_pre6: additional ESProducer in cfg files
if cmsswVersion > "CMSSW_4_3":
    process.hltESPStripLorentzAngleDep = cms.ESProducer("SiStripLorentzAngleDepESProducer",
        LatencyRecord = cms.PSet(
            record = cms.string('SiStripLatencyRcd'),
            label = cms.untracked.string('')
        ),
        LorentzAngleDeconvMode = cms.PSet(
            record = cms.string('SiStripLorentzAngleRcd'),
            label = cms.untracked.string('deconvolution')
        ),
        LorentzAnglePeakMode = cms.PSet(
            record = cms.string('SiStripLorentzAngleRcd'),
            label = cms.untracked.string('peak')
        )
    )

# from CMSSW_4_3_0_pre6: ECAL severity flags migration
if cmsswVersion > "CMSSW_4_3":
  import HLTrigger.Configuration.Tools.updateEcalSeverityFlags
  HLTrigger.Configuration.Tools.updateEcalSeverityFlags.update( process.__dict__ )

