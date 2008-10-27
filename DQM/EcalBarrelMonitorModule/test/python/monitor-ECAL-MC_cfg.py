import FWCore.ParameterSet.Config as cms

process = cms.Process("ECALDQM")

import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.EcalMapping.EcalMapping_cfi")

process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

process.load("DQM.EcalBarrelMonitorModule.EcalBarrelMonitorModule_cfi")

process.load("DQM.EcalBarrelMonitorTasks.EcalBarrelMonitorTasks_cfi")

process.load("DQM.EcalBarrelMonitorTasks.mergeRuns_cff")

process.load("DQM.EcalBarrelMonitorClient.EcalBarrelMonitorClient_cfi")

process.load("DQM.EcalEndcapMonitorModule.EcalEndcapMonitorModule_cfi")

process.load("DQM.EcalEndcapMonitorTasks.EcalEndcapMonitorTasks_cfi")

process.load("DQM.EcalEndcapMonitorTasks.mergeRuns_cff")

process.load("DQM.EcalEndcapMonitorClient.EcalEndcapMonitorClient_cfi")

process.load("RecoEcal.EgammaClusterProducers.ecalClusteringSequence_cff")

process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")

process.load("DQMServices.Core.DQM_cfg")

process.preScaler = cms.EDFilter("Prescaler",
    prescaleFactor = cms.int32(1)
)

process.ecalPrescaler = cms.EDFilter("EcalMonitorPrescaler",
    EcalRawDataCollection = cms.InputTag("ecalEBunpacker"),
    laserPrescaleFactor = cms.untracked.int32(1),
    ledPrescaleFactor = cms.untracked.int32(1),
    pedestalPrescaleFactor = cms.untracked.int32(1),
    testpulsePrescaleFactor = cms.untracked.int32(1)
)

process.dqmInfoEB = cms.EDFilter("DQMEventInfo",
    subSystemFolder = cms.untracked.string('EcalBarrel')
)

process.dqmInfoEE = cms.EDFilter("DQMEventInfo",
    subSystemFolder = cms.untracked.string('EcalEndcap')
)

process.dqmSaver = cms.EDFilter("DQMFileSaver",
    dirName = cms.untracked.string('.'),
    convention = cms.untracked.string('Online')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(150)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/users/dellaric/data/relval_zee.root')
)

process.EcalTrivialConditionRetriever = cms.ESSource("EcalTrivialConditionRetriever",
    adcToGeVEBConstant = cms.untracked.double(0.035),
    adcToGeVEEConstant = cms.untracked.double(0.06),
    pedWeights = cms.untracked.vdouble(0.333, 0.333, 0.333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    amplWeights = cms.untracked.vdouble(-0.333, -0.333, -0.333, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
    jittWeights = cms.untracked.vdouble(0.041, 0.041, 0.041, 0.0, 1.325, -0.05, -0.504, -0.502, -0.390, 0.0)
)

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING'),
        noLineBreaks = cms.untracked.bool(True),
        noTimeStamps = cms.untracked.bool(True),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalBarrelMonitorModule = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalEndcapMonitorModule = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    categories = cms.untracked.vstring('EcalBarrelMonitorModule',
                                       'EcalEndcapMonitorModule'),
    destinations = cms.untracked.vstring('cout')
)

process.ecalDataSequence = cms.Sequence(process.preScaler*process.ecalUncalibHit*process.ecalRecHit)

process.ecalBarrelMonitorSequence = cms.Sequence(process.ecalBarrelMonitorModule*process.dqmInfoEB*process.ecalBarrelMonitorClient)

process.ecalEndcapMonitorSequence = cms.Sequence(process.ecalEndcapMonitorModule*process.dqmInfoEE*process.ecalEndcapMonitorClient)

process.p = cms.Path(process.ecalDataSequence*process.ecalBarrelMonitorSequence*process.ecalEndcapMonitorSequence*process.dqmSaver)
process.q = cms.Path(process.ecalDataSequence*~process.ecalPrescaler*process.hybridSuperClusters*process.correctedHybridSuperClusters*process.multi5x5BasicClusters*process.multi5x5SuperClusters)
process.r = cms.EndPath(process.ecalBarrelDefaultTasksSequence*process.ecalBarrelClusterTask*process.ecalEndcapDefaultTasksSequence*process.ecalEndcapClusterTask)

process.ecalUncalibHit.MinAmplBarrel = 12.
process.ecalUncalibHit.MinAmplEndcap = 16.
process.ecalUncalibHit.EBdigiCollection = 'simEcalDigis:ebDigis'
process.ecalUncalibHit.EEdigiCollection = 'simEcalDigis:eeDigis'

process.ecalRecHit.EBuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEB'
process.ecalRecHit.EEuncalibRecHitCollection = 'ecalUncalibHit:EcalUncalibRecHitsEE'

process.ecalBarrelMonitorModule.mergeRuns = True
process.ecalBarrelMonitorModule.EBDigiCollection = 'simEcalDigis:ebDigis'
process.ecalBarrelMonitorModule.runType = 3
process.ecalBarrelMonitorModule.EcalTrigPrimDigiCollection = 'ecalTriggerPrimitiveDigis'

process.ecalBarrelOccupancyTask.EBDigiCollection = 'simEcalDigis:ebDigis'
process.ecalBarrelOccupancyTask.EcalTrigPrimDigiCollection = 'ecalTriggerPrimitiveDigis'

process.ecalBarrelPedestalOnlineTask.EBDigiCollection = 'simEcalDigis:ebDigis'

process.ecalBarrelMonitorClient.maskFile = '../data/maskfile-EB.dat'
process.ecalBarrelMonitorClient.mergeRuns = True
process.ecalBarrelMonitorClient.location = 'H4'
process.ecalBarrelMonitorClient.baseHtmlDir = 'HTML_EB'
process.ecalBarrelMonitorClient.enabledClients = ['Integrity', 'Occupancy', 'PedestalOnline', 'Timing', 'Cluster', 'Summary']

process.ecalEndcapMonitorModule.mergeRuns = True
process.ecalEndcapMonitorModule.EEDigiCollection = 'simEcalDigis:eeDigis'
process.ecalEndcapMonitorModule.runType = 3
process.ecalEndcapMonitorModule.EcalTrigPrimDigiCollection = 'ecalTriggerPrimitiveDigis'

process.ecalEndcapOccupancyTask.EEDigiCollection = 'simEcalDigis:eeDigis'
process.ecalEndcapOccupancyTask.EcalTrigPrimDigiCollection = 'ecalTriggerPrimitiveDigis'

process.ecalEndcapPedestalOnlineTask.EEDigiCollection = 'simEcalDigis:eeDigis'

process.ecalEndcapMonitorClient.maskFile = '../data/maskfile-EE.dat'
process.ecalEndcapMonitorClient.mergeRuns = True
process.ecalEndcapMonitorClient.location = 'H4'
process.ecalEndcapMonitorClient.baseHtmlDir = 'HTML_EE'
process.ecalEndcapMonitorClient.enabledClients = ['Integrity', 'Occupancy', 'PedestalOnline', 'Timing', 'Cluster', 'Summary']

process.DQM.collectorHost = ''

