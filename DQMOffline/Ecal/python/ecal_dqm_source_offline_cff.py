import FWCore.ParameterSet.Config as cms

from DQM.EcalBarrelMonitorModule.EcalBarrelMonitorModule_cfi import *
from DQM.EcalBarrelMonitorTasks.EcalBarrelMonitorTasks_cfi import *

dqmInfoEB = cms.EDFilter("DQMEventInfo",
    subSystemFolder = cms.untracked.string('EcalBarrel')
)

ecal_dqm_source_offline1 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelIntegrityTask)
ecal_dqm_source_offline2 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask)
ecal_dqm_source_offline3 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelCosmicTask)
ecal_dqm_source_offline4 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelCosmicTask)
ecal_dqm_source_offline5 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelPedestalOnlineTask*ecalBarrelCosmicTask)
ecal_dqm_source_offline6 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelPedestalOnlineTask*ecalBarrelTriggerTowerTask*ecalBarrelCosmicTask)
ecal_dqm_source_offline7 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelPedestalOnlineTask*ecalBarrelTriggerTowerTask*ecalBarrelCosmicTask*ecalBarrelClusterTask)
ecal_dqm_source_offline9 = cms.Sequence(ecalBarrelMonitorModule*dqmInfoEB*ecalBarrelOccupancyTask*ecalBarrelIntegrityTask*ecalBarrelStatusFlagsTask*ecalBarrelLaserTask*ecalBarrelPedestalOnlineTask*ecalBarrelPedestalTask*ecalBarrelTestPulseTask*ecalBarrelTriggerTowerTask*ecalBarrelTimingTask*ecalBarrelCosmicTask*ecalBarrelClusterTask)

ecal_dqm_source_offline = cms.Sequence(ecal_dqm_source_offline1)

ecalBarrelMonitorModule.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelMonitorModule.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
ecalBarrelMonitorModule.EcalTrigPrimDigiCollection = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
ecalBarrelMonitorModule.verbose = False

ecalBarrelCosmicTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelCosmicTask.EcalUncalibratedRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEB")

ecalBarrelIntegrityTask.EBDetIdCollection0 = cms.InputTag("ecalDigis","EcalIntegrityDCCSizeErrors")
ecalBarrelIntegrityTask.EBDetIdCollection1 = cms.InputTag("ecalDigis","EcalIntegrityGainErrors")
ecalBarrelIntegrityTask.EBDetIdCollection2 = cms.InputTag("ecalDigis","EcalIntegrityChIdErrors")
ecalBarrelIntegrityTask.EBDetIdCollection3 = cms.InputTag("ecalDigis","EcalIntegrityGainSwitchErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection1 = cms.InputTag("ecalDigis","EcalIntegrityTTIdErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection2 = cms.InputTag("ecalDigis","EcalIntegrityBlockSizeErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection3 = cms.InputTag("ecalDigis","EcalIntegrityMemTtIdErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection4 = cms.InputTag("ecalDigis","EcalIntegrityMemBlockSizeErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection5 = cms.InputTag("ecalDigis","EcalIntegrityMemChIdErrors")
ecalBarrelIntegrityTask.EcalElectronicsIdCollection6 = cms.InputTag("ecalDigis","EcalIntegrityMemGainErrors")

ecalBarrelLaserTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelLaserTask.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
ecalBarrelLaserTask.EcalPnDiodeDigiCollection = cms.InputTag("ecalDigis")
ecalBarrelLaserTask.EcalUncalibratedRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEB")

ecalBarrelOccupancyTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelOccupancyTask.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
ecalBarrelOccupancyTask.EcalPnDiodeDigiCollection = cms.InputTag("ecalDigis")
ecalBarrelOccupancyTask.EcalTrigPrimDigiCollection = cms.InputTag("ecalDigis","EcalTriggerPrimitives")

ecalBarrelPedestalOnlineTask.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")

ecalBarrelPedestalTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelPedestalTask.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
ecalBarrelPedestalTask.EcalPnDiodeDigiCollection = cms.InputTag("ecalDigis")

ecalBarrelStatusFlagsTask.EcalRawDataCollection = cms.InputTag("ecalDigis")

ecalBarrelTestPulseTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelTestPulseTask.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
ecalBarrelTestPulseTask.EcalPnDiodeDigiCollection = cms.InputTag("ecalDigis")
ecalBarrelTestPulseTask.EcalUncalibratedRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEB")

ecalBarrelTimingTask.EcalRawDataCollection = cms.InputTag("ecalDigis")
ecalBarrelTimingTask.EcalUncalibratedRecHitCollection = cms.InputTag("ecalFixedAlphaBetaFitUncalibRecHit","EcalUncalibRecHitsEB")

ecalBarrelTriggerTowerTask.EcalTrigPrimDigiCollectionReal = cms.InputTag("ecalDigis","EcalTriggerPrimitives")
ecalBarrelTriggerTowerTask.EcalTrigPrimDigiCollectionEmul = 'valEcalTriggerPrimitiveDigis'

ecalBarrelClusterTask.BasicClusterCollection = cms.InputTag("cosmicBasicClusters","CosmicBarrelBasicClusters")

