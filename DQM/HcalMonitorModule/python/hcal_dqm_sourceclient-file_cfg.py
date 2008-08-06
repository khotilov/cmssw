import FWCore.ParameterSet.Config as cms

process = cms.Process("HCALDQM")
#
# BEGIN DQM Online Environment ###########################
#
process.load("DQMServices.Core.DQM_cfg")

# use include file for dqmEnv dqmSaver
process.load("DQMServices.Components.DQMEnvironment_cfi")

#
# END ################################################
#
process.load("DQM.HcalMonitorModule.HcalMonitorModule_cfi")

process.load("DQM.HcalMonitorClient.HcalMonitorClient_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hbhe_cfi")

process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_ho_cfi")

process.load("RecoLocalCalo.HcalRecProducers.HcalSimpleReconstructor_hf_cfi")

# will add switch to select histograms to be saved soon
# Quality Tester #### 
process.qTester = cms.EDFilter("QualityTester",
    prescaleFactor = cms.untracked.int32(1),
    qtList = cms.untracked.FileInPath('DQM/HcalMonitorClient/data/hcal_qualitytest_config.xml'),
    getQualityTestsFromFile = cms.untracked.bool(True)
)

process.prefer("GlobalTag")
process.Timing = cms.Service("Timing")

#path p = {hcalDigis, horeco, hfreco, hbhereco, hcalMonitor}
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound', 
        'TooManyProducts', 
        'TooFewProducts')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
#To run on HCAL local runs:  (you need to have access to /bigspool or modify the path to the datafile)	
#source = HcalTBSource {
#	untracked vstring fileNames = {'file:/afs/cern.ch/user/s/stjohn/scratch0/USC_044790.root'}
#untracked vstring streams = { 
#HBHE a, b, c:
#'HCAL_DCC700','HCAL_DCC701','HCAL_DCC702','HCAL_DCC703','HCAL_DCC704','HCAL_DCC705',
#'HCAL_DCC706','HCAL_DCC707','HCAL_DCC708','HCAL_DCC709','HCAL_DCC710','HCAL_DCC711',
#'HCAL_DCC712','HCAL_DCC713','HCAL_DCC714','HCAL_DCC715','HCAL_DCC716','HCAL_DCC717',
#HF:
#'HCAL_DCC718','HCAL_DCC719','HCAL_DCC720','HCAL_DCC721','HCAL_DCC722','HCAL_DCC723',
#HO:
#'HCAL_DCC724','HCAL_DCC725','HCAL_DCC726','HCAL_DCC727','HCAL_DCC728','HCAL_DCC729',
#'HCAL_DCC730','HCAL_DCC731',
#'HCAL_Trigger','HCAL_SlowData'
#}
#}	
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/CRUZET3/Cosmics/RAW/v1/000/051/199/EA1F908F-AD4E-DD11-8235-000423D6A6F4.root')
)

process.p = cms.Path(process.hcalDigis*process.horeco*process.hfreco*process.hbhereco*process.hcalMonitor*process.hcalClient*process.dqmEnv*process.dqmSaver)
process.DQM.collectorHost = 'myhost'
process.DQM.collectorPort = 9092
process.dqmSaver.convention = 'Online'
#replace dqmSaver.dirName          = "."
process.dqmSaver.producer = 'DQM'
process.dqmEnv.subSystemFolder = 'Hcal'
# optionally change fileSaving  conditions
# replace dqmSaver.prescaleLS =   -1
# replace dqmSaver.prescaleTime = -1 # in minutes
# replace dqmSaver.prescaleEvt =  -1
#> For Hcal local run files, replace dqmSaver.saveByRun = 2 
process.dqmSaver.saveByRun = 1
process.hcalMonitor.DigiOccThresh = -999999999 ##Temporary measure while DigiOcc is reworked.

process.hcalMonitor.debug = False
process.hcalMonitor.PedestalsPerChannel = False
process.hcalMonitor.DataFormatMonitor = True
process.hcalMonitor.DigiMonitor = False
process.hcalMonitor.PedestalMonitor = False
process.hcalMonitor.LEDMonitor = False
process.hcalMonitor.TrigPrimMonitor = True
process.hcalMonitor.HotCellMonitor = False
process.hcalMonitor.DeadCellMonitor = False
process.hcalMonitor.RecHitMonitor = False
process.hcalMonitor.MTCCMonitor = False
process.hcalMonitor.CaloTowerMonitor = False
process.hcalMonitor.HcalAnalysis = False
process.hcalClient.plotPedRAW = True
process.hcalClient.DoPerChanTests = False
process.hcalClient.SummaryClient = True
process.hcalClient.DataFormatClient = True
process.hcalClient.DigiClient = False
process.hcalClient.LEDClient = False
process.hcalClient.PedestalClient = False
process.hcalClient.TrigPrimClient = True
process.hcalClient.RecHitClient = False
process.hcalClient.HotCellClient = False
process.hcalClient.DeadCellClient = False
process.hcalClient.CaloTowerClient = False
process.GlobalTag.globaltag = 'CRUZET3_V2P::All'
process.GlobalTag.connect = 'frontier://Frontier/CMS_COND_20X_GLOBALTAG' ##Frontier/CMS_COND_20X_GLOBALTAG"


