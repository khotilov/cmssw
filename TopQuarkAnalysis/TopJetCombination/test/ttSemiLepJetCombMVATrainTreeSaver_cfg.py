import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TtSemiLepJetCombMVATrainer')
process.MessageLogger.cerr.TtSemiLepJetCombMVATrainer = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_4_0_pre1/RelValTTbar/GEN-SIM-RECO/MC_31X_V9-v1/0008/2C8CD8FE-B6B5-DE11-ACB8-001D09F2905B.root'
    )
)
## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry & conditions
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MC_3XY_V10::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## std sequence for ttGenEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## configure ttGenEventFilters
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEventFilters_cff")
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = False
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = False

## configure jet-parton matching
process.load("TopQuarkAnalysis.TopTools.TtSemiLepJetPartonMatch_cfi")
#process.ttSemiLepJetPartonMatch.partonsToIgnore = ["LepB"]

## configure mva trainer
process.load("TopQuarkAnalysis.TopJetCombination.TtSemiLepJetCombMVATrainTreeSaver_cff")
## change maximum number of jets taken into account per event (default: 4)
#process.ttSemiLepJetPartonMatch .maxNJets = 5
#process.trainTtSemiLepJetCombMVA.maxNJets = process.ttSemiLepJetPartonMatch.maxNJets

## make trainer looper known to the process
from TopQuarkAnalysis.TopJetCombination.TtSemiLepJetCombMVATrainTreeSaver_cff import looper
process.looper = looper

process.p = cms.Path(process.makeGenEvt *
                     process.patDefaultSequence *
                     process.ttSemiLeptonicFilter *
                     process.ttSemiLepJetPartonMatch *
                     process.saveTtSemiLepJetCombMVATrainTree)
