import FWCore.ParameterSet.Config as cms

process = cms.Process("ECALDQM")

process.load("EventFilter.EcalTBRawToDigi.EcalUnpackerData_cfi")

import RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi
process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi.ecalFixedAlphaBetaFitUncalibRecHit.clone()

#import RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi
#process.ecalUncalibHit = RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi.ecalWeightUncalibRecHit.clone()

process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.EcalMapping.EcalMapping_cfi")

process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

process.load("DQM.EcalBarrelMonitorModule.EcalBarrelMonitorModule_cfi")

process.load("DQM.EcalBarrelMonitorTasks.EcalBarrelMonitorTasks_cfi")

process.load("DQM.EcalBarrelMonitorClient.EcalBarrelMonitorClient_cfi")

process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")

process.load("DQMServices.Core.DQM_cfg")

process.preScaler = cms.EDFilter("Prescaler",
    prescaleFactor = cms.int32(1)
)

process.dqmInfoEB = cms.EDFilter("DQMEventInfo",
    subSystemFolder = cms.untracked.string('EcalBarrel')
)

process.dqmQTestEB = cms.EDFilter("QualityTester",
    reportThreshold = cms.untracked.string('red'),
    prescaleFactor = cms.untracked.int32(1),
    qtList = cms.untracked.FileInPath('DQM/EcalBarrelMonitorModule/test/data/EcalBarrelQualityTests.xml'),
    getQualityTestsFromFile = cms.untracked.bool(True)
)

process.dqmSaverEB = cms.EDFilter("DQMFileSaver",
    dirName = cms.untracked.string('.'),
    fileName = cms.untracked.string('EcalBarrel'),
    convention = cms.untracked.string('Online'),
# special fix for H4 cosmic runs
    saveByRun = cms.untracked.int32(2)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(150)
)
process.source = cms.Source("DAQEcalTBInputService",
#---
#    fileNames = cms.untracked.vstring('/store/users/dellaric/data/H4-000006100-SM5-COSMIC-STD'),
#    runNumber = cms.untracked.uint32(6100),
#    isBinary = cms.untracked.bool(False)
#--
#    fileNames = cms.untracked.vstring('/store/users/dellaric/data/H4-000008203-SM5-PEDESTAL-STD'),
#    runNumber = cms.untracked.uint32(8203),
#    isBinary = cms.untracked.bool(False)
#---
    fileNames = cms.untracked.vstring('/store/users/dellaric/data/H4-000008205-SM5-TEST_PULSE-MGPA'),
    runNumber = cms.untracked.uint32(8205),
    isBinary = cms.untracked.bool(False)
#---
#    fileNames = cms.untracked.vstring('/store/users/dellaric/data/H4-000008206-SM5-PEDESTAL-STD'),
#    runNumber = cms.untracked.uint32(8206),
#    isBinary = cms.untracked.bool(False)
#---
#    fileNames = cms.untracked.vstring('/store/users/dellaric/data/H4-000015188-SM5-PEDESTAL-STD'),
#    runNumber = cms.untracked.uint32(15188),
#    isBinary = cms.untracked.bool(True)
#---
)

process.EcalTrivialConditionRetriever = cms.ESSource("EcalTrivialConditionRetriever",
    adcToGeVEBConstant = cms.untracked.double(0.00875), ## 0.035
    pedWeights = cms.untracked.vdouble(0.333, 0.333, 0.333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    amplWeights = cms.untracked.vdouble(-0.333, -0.333, -0.333, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
    jittWeights = cms.untracked.vdouble(0.041, 0.041, 0.041, 0.0, 1.325, -0.05, -0.504, -0.502, -0.39, 0.0)
)

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        EcalTBRawToDigiDccBlockSize = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevTowerId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiNumTowerBlocks = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        threshold = cms.untracked.string('WARNING'),
        noLineBreaks = cms.untracked.bool(True),
        EcalRawToDigiDevSRP = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiChId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiTpg = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDev = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevMemGain = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiGainSwitch = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevMemBlock = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevDccBlockSize = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevTriggerType = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevTpg = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevNumTowerBlocks = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiTowerId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalBarrelMonitorModule = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiTowerSize = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevGainSwitch = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBInputService = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalTBRawToDigiGainZero = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevTowerSize = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevMemChId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        noTimeStamps = cms.untracked.bool(True),
        EcalRawToDigiDevChId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevGainZero = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigiTriggerType = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalRawToDigiDevTCC = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalDCCHeaderRuntypeDecoder = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalRawToDigiDevMemTowerId = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        ),
        EcalTBRawToDigi = cms.untracked.PSet(
            limit = cms.untracked.int32(1000)
        )
    ),
    categories = cms.untracked.vstring('EcalTBInputService', 
        'EcalTBRawToDigi', 
        'EcalTBRawToDigiTriggerType', 
        'EcalTBRawToDigiTpg', 
        'EcalTBRawToDigiNumTowerBlocks', 
        'EcalTBRawToDigiTowerId', 
        'EcalTBRawToDigiTowerSize', 
        'EcalTBRawToDigiChId', 
        'EcalTBRawToDigiGainZero', 
        'EcalTBRawToDigiGainSwitch', 
        'EcalTBRawToDigiDccBlockSize', 
        'EcalRawToDigiDev', 
        'EcalRawToDigiDevTriggerType', 
        'EcalRawToDigiDevTpg', 
        'EcalRawToDigiDevNumTowerBlocks', 
        'EcalRawToDigiDevTowerId', 
        'EcalRawToDigiDevTowerSize', 
        'EcalRawToDigiDevChId', 
        'EcalRawToDigiDevGainZero', 
        'EcalRawToDigiDevGainSwitch', 
        'EcalRawToDigiDevDccBlockSize', 
        'EcalRawToDigiDevMemBlock', 
        'EcalRawToDigiDevMemTowerId', 
        'EcalRawToDigiDevMemChId', 
        'EcalRawToDigiDevMemGain', 
        'EcalRawToDigiDevTCC', 
        'EcalRawToDigiDevSRP', 
        'EcalDCCHeaderRuntypeDecoder', 
        'EcalBarrelMonitorModule'),
    destinations = cms.untracked.vstring('cout')
)

process.ecalDataSequence = cms.Sequence(process.preScaler*process.ecalEBunpacker*process.ecalUncalibHit*process.ecalRecHit)
process.ecalBarrelMonitorSequence = cms.Sequence(process.ecalBarrelMonitorModule*process.dqmInfoEB*process.ecalBarrelMonitorClient*process.dqmQTestEB*process.dqmSaverEB)

process.p = cms.Path(process.ecalDataSequence*process.ecalBarrelMonitorSequence)
process.q = cms.EndPath(process.ecalBarrelCosmicTasksSequence)

process.ecalUncalibHit.MinAmplBarrel = 12.
process.ecalUncalibHit.MinAmplEndcap = 16.
process.ecalUncalibHit.EBdigiCollection = cms.InputTag("ecalEBunpacker","ebDigis")
process.ecalUncalibHit.EEdigiCollection = cms.InputTag("ecalEBunpacker","eeDigis")

process.ecalRecHit.EBuncalibRecHitCollection = cms.InputTag("ecalUncalibHit","EcalUncalibRecHitsEB")
process.ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("ecalUncalibHit","EcalUncalibRecHitsEE")

process.ecalBarrelMonitorClient.maskFile = ''
process.ecalBarrelMonitorClient.location = 'H4'
process.ecalBarrelMonitorClient.baseHtmlDir = '.'
process.ecalBarrelMonitorClient.superModules = [19]

