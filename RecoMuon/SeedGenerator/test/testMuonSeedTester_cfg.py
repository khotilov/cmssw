import FWCore.ParameterSet.Config as cms

process = cms.Process("RecoMuon")
# Messages
#process.load("RecoMuon.Configuration.MessageLogger_cfi")
#process.load("FWCore.MessageService.MessageLogger_cfi")

# Muon Reco
process.load("RecoLocalMuon.Configuration.RecoLocalMuon_cff")

process.load("RecoMuon.Configuration.RecoMuon_cff")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")


process.source = cms.Source("PoolSource",
#                           fileNames = cms.untracked.vstring('/store/relval/CMSSW_2_1_9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/IDEAL_V9_v2/0000/2A00EECC-A185-DD11-93A9-000423D9517C.root')
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/s/stoyan/mc/ingo/cmssw2_1_9-mumin_e_60_300_probev2__1.root')
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('RecoMuons.root')
                               )

process.STAMuonAnalyzer = cms.EDFilter("STAMuonAnalyzer",
    DataType = cms.untracked.string('SimData'),
    StandAloneTrackCollectionLabel = cms.untracked.string('standAloneMuons'),
#    MuonSeedCollectionLabel = cms.untracked.string('MuonSeedTester'),
    MuonSeedCollectionLabel = cms.untracked.string('SETMuonSeed'),
    rootFileName = cms.untracked.string('STAMuonAnalyzer.root'),
)


#process.p = cms.Path(process.MuonSeed*process.standAloneMuons)
#process.p = cms.Path(process.MuonSeedTester*process.standAloneMuons*process.STAMuonAnalyzer)
process.p = cms.Path(process.SETMuonSeed*process.standAloneMuons*process.STAMuonAnalyzer)
#process.this_is_the_end = cms.EndPath(process.out)

process.GlobalTag.globaltag = 'IDEAL_V9::All'


