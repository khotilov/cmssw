import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
## dump content of TopGenEvent
process.MessageLogger.categories.append('TopGenEvent')

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_3_0_pre4/RelValTTbar/GEN-SIM-RECO/MC_31X_V8-v1/0000/FCC5D1DF-C3A6-DE11-81AA-001D09F28EA3.root'
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
process.GlobalTag.globaltag = cms.string('MC_31X_V8::All')

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

from TopQuarkAnalysis.Examples.TopGenEventAnalyzer_cfi import analyzeTopGenEvent
process.analyzeTopGenEvent = analyzeTopGenEvent

# register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeTopGenEvent.root')
)

## end path   
process.p1 = cms.Path(process.patDefaultSequence *
                      process.makeGenEvt *
                      process.analyzeTopGenEvent)

