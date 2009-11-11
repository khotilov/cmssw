import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

# Some generic services and conditions data
process.Timing = cms.Service("Timing")
process.Tracer = cms.Service("Tracer",sourceSeed = cms.untracked.string("$$"))
#process.load("DQM.SiStripCommon.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP3X_V11::All')

# Input files
castor = cms.untracked.vstring()
local  = cms.untracked.vstring()
process.source = cms.Source (
    "PoolSource",
    fileNames = castor,
    )
castor.extend( [
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/FA7139E8-97BD-DE11-A3E2-002618943935.root',
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/BC3224A5-9ABD-DE11-A625-002354EF3BDB.root',
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/8C578DA3-C0BD-DE11-9DEA-0017312A250B.root',
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/7A29EA77-9DBD-DE11-A3BC-0026189438ED.root',
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/3EA8A506-10BE-DE11-BB21-0018F3D09704.root',
    '/store/relval/CMSSW_3_4_0_pre2/RelValQCD_Pt_80_120/GEN-SIM-RECO/STARTUP3XY_V9-v1/0003/04383FF7-9EBD-DE11-8511-0018F3D09616.root',
    ] )
local.extend( [
    'file:/data/bainbrid/jpt/FA7139E8-97BD-DE11-A3E2-002618943935.root',
    'file:/data/bainbrid/jpt/BC3224A5-9ABD-DE11-A625-002354EF3BDB.root',
    'file:/data/bainbrid/jpt/8C578DA3-C0BD-DE11-9DEA-0017312A250B.root',
    'file:/data/bainbrid/jpt/7A29EA77-9DBD-DE11-A3BC-0026189438ED.root',
    'file:/data/bainbrid/jpt/3EA8A506-10BE-DE11-BB21-0018F3D09704.root',
    'file:/data/bainbrid/jpt/04383FF7-9EBD-DE11-8511-0018F3D09616.root',
    ] )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# Identify GenParticles to be used to build GenJets (ie, no neutrinos or BSM)
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 12, 14, 16
    )
process.genParticlesForJets.excludeFromResonancePids = cms.vuint32(12, 14, 16)

# Build reco::GenJets from GenParticles
from RecoJets.JetProducers.ic5GenJets_cfi import iterativeCone5GenJets
process.iterativeCone5GenJetsNoNuBSM = iterativeCone5GenJets.clone()

# ZSP and JPT corrections
process.load("JetMETCorrections.Configuration.ZSPJetCorrections219_cff")
process.load("JetMETCorrections.Configuration.JetPlusTrackCorrections_cff")
#process.JetPlusTrackZSPCorrectorIcone5.VectorialCorrection  = cms.bool(True)
#process.JetPlusTrackZSPCorrectorIcone5.UseResponseInVecCorr = cms.bool(False)

# Analyzer module
process.myanalysis = cms.EDAnalyzer(
    "JPTAnalyzer",
    HistOutFile      = cms.untracked.string('analysis.root'),
    calojets         = cms.string('iterativeCone5CaloJets'),
    zspjets          = cms.string('ZSPJetCorJetIcone5'),
    genjets          = cms.string('iterativeCone5GenJetsNoNuBSM'),
    JetCorrectionJPT = cms.string('JetPlusTrackZSPCorrectorIcone5'),
    )

# Path (GenJets + ZSP + EID + JTA + analysis)
process.p1 = cms.Path(
    process.genParticlesForJets *
    process.iterativeCone5GenJetsNoNuBSM *
    process.ZSPJetCorrections *
    process.ZSPrecoJetAssociationsIcone5 *
    process.myanalysis
    )
