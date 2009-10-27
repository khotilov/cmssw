import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.categories.append('TtSemiLeptonicEvent')
process.MessageLogger.cerr.TtSemiLeptonicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_4_0_pre1/RelValTTbar/GEN-SIM-RECO/MC_31X_V9-v1/0008/2C8CD8FE-B6B5-DE11-ACB8-001D09F2905B.root'
     ),
     skipEvents = cms.untracked.uint32(0)                            
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
process.GlobalTag.globaltag = cms.string('MC_3XY_V12::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## sequences for ttGenEvent and TtSemiLeptonicEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff")
## enable additional per-event printout from the TtSemiLeptonicEvent
#process.ttSemiLepEvent.verbosity = 1
## change maximum number of jets taken into account per event (default: 4)
#process.ttSemiLepEvent.maxNJets = 5
## change jet-parton matching algorithm
#process.ttSemiLepJetPartonMatch.algorithm = "unambiguousOnly"

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttSemiLepEvtBuilder_cff import addTtSemiLepHypotheses
addTtSemiLepHypotheses(process,
                       ["kMaxSumPtWMass", "kMVADisc"]
                       )

## load HypothesisAnalyzer
process.load("TopQuarkAnalysis.Examples.HypothesisAnalyzer_cff")

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTopHypothesis.root')
)

## end path   
process.path = cms.Path(process.patDefaultSequence *
                        process.makeGenEvt *
                        process.makeTtSemiLepEvent *
                        process.analyzeHypotheses)
