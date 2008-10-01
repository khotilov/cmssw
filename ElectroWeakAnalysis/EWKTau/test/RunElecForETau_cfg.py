import FWCore.ParameterSet.Config as cms

process = cms.Process("Comm")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("DQMOffline.Trigger.Tau.HLTTauDQMOffline_cff")
process.TauRefProducer.Electrons.ElectronCollection = "allLayer0ElecForETau"
process.TauRefProducer.Electrons.IdCollection = "elecIdForETauCutBasedRobust"
process.TauRefProducer.Electrons.doElectrons = False
process.TauRefProducer.Electrons.doElecFromZ = True
process.TauRefProducer.Electrons.doTrackIso = False
process.TauRefProducer.Electrons.doID = True
process.TauRefProducer.Electrons.ElecEtFromZcut = 15.
process.TauRefProducer.Electrons.ptMin = 15.

process.load("ElectroWeakAnalysis.EWKTau.elecForETauMod.elecForETauPatProducer_cff")
process.load("ElectroWeakAnalysis.EWKTau.tauForETauMod.tauForETauPFPatProducer_cff")
process.load("ElectroWeakAnalysis.EWKTau.analyzerForETau_cfi")
process.load("ElectroWeakAnalysis.EWKTau.metForETauPatConfig_cff")

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeETau_test1.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(6205)         
)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/0C3B40D7-F87D-DD11-A9FB-000423D998BA.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/3A5455F3-F87D-DD11-AEF4-000423D94534.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/420F04EB-F87D-DD11-95CF-001617E30D0A.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/42A0DDBC-F87D-DD11-BBB3-000423D98800.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/4E209DDA-F87D-DD11-965A-000423D98834.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/5201CDC0-F87D-DD11-A0DC-001617C3B710.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/702355EE-F87D-DD11-8502-000423D6B48C.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/7AAE47C1-F87D-DD11-9F0D-000423D94AA8.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/84FCF4D4-F87D-DD11-9373-000423D8FA38.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/869D29B7-F87D-DD11-8387-001617C3B5D8.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/86E7A6C2-F87D-DD11-88FE-000423D951D4.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/90A7FBDF-F87D-DD11-B0FE-000423D99A8E.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/BA8F86D9-F87D-DD11-B1BE-000423D987FC.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/D28EBFDE-F87D-DD11-8F47-000423D6B42C.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/E20313D4-F87D-DD11-A3B8-001617DC1F70.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0001/FE2C34C8-F87D-DD11-B0BF-000423D94534.root',
        '/store/relval/CMSSW_2_1_7/RelValZEE/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V7_v1/0002/2AA81C2A-437E-DD11-8C3C-000423DD2F34.root'

     )
)


process.p = cms.Path(
                      process.layer0ElecForETau
                     *process.layer1ElecForETau
                     *process.PFTauForETauPreamble
                     *process.PFTauForETauPat
                     *process.patMET
                     *process.TauRefProducer
                     *process.analyzeETau
                     )


