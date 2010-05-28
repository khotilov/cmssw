import FWCore.ParameterSet.Config as cms

process = cms.Process("TopDQM")
process.load("DQM.Physics.topSingleLeptonDQM_cfi")
process.load("DQM.Physics.topDiLeptonOfflineDQM_cfi")

process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.DQM.collectorHost = ''

process.dqmSaver.workflow = cms.untracked.string('/Physics/TopSingleLeptonDQM/DataSet')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource"
    ,fileNames = cms.untracked.vstring(
     # '/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/1A19B479-BA3A-DF11-8E43-0017A4770410.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/1A2CED78-BA3A-DF11-98CD-0017A4771010.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/3AE61B7A-BA3A-DF11-BA4C-0017A477040C.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/3CBA7F7C-BA3A-DF11-9ECE-0017A4770C14.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/443CAD79-BA3A-DF11-9F90-0017A4770818.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/4C91A47A-BA3A-DF11-B3D2-0017A4771004.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/5225C429-BB3A-DF11-AD90-0017A4770020.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/62BC7102-BB3A-DF11-8D7C-0017A4771028.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/64FCA77B-BA3A-DF11-8514-0017A477042C.root'
     #,'/store/mc/Spring10/TTbar/GEN-SIM-RECO/MC_3XY_V25_S09_preproduction-v2/0106/7AE57478-BA3A-DF11-BA3C-0017A4771034.root'



 '/store/data/Run2010A/Mu/RECO/v2/000/136/089/D802058C-6E67-DF11-B721-000423D98634.root',
 '/store/data/Run2010A/Mu/RECO/v2/000/136/087/F63A8D4C-6167-DF11-BD3A-000423D944F8.root',
 '/store/data/Run2010A/Mu/RECO/v2/000/136/087/7AE9EF4A-7D67-DF11-97EC-0030487CD7CA.root'



      
    )
)

## load jet corrections
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
process.prefer("ak5CaloL2L3")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TopSingleLeptonDQM'   )
process.MessageLogger.cerr.TopSingleLeptonDQM    = cms.untracked.PSet(limit = cms.untracked.int32(1))
process.MessageLogger.categories.append('TopDiLeptonOfflineDQM')
process.MessageLogger.cerr.TopDiLeptonOfflineDQM = cms.untracked.PSet(limit = cms.untracked.int32(1))

#process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(#process.content *
                     process.topSingleMuonMediumDQM +
                     process.topSingleElectronMediumDQM +
                     process.topSingleMuonLooseDQM +
                     process.topSingleElectronLooseDQM +
                     process.topDiLeptonOfflineDQM +
                     process.topSingleLeptonDQM +
                     process.dqmSaver
                     )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
