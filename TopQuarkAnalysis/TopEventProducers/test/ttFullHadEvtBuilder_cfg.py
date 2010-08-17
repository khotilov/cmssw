import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

## add message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TtFullHadronicEvent')
process.MessageLogger.cerr.TtFullHadronicEvent = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_6_0/RelValTTbar/GEN-SIM-RECO/START36_V4-v1/0014/EEA7EEC1-FC49-DF11-9E91-003048678D9A.root'
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

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")

## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START36_V4::All')

## load magnetic field
process.load("Configuration.StandardSequences.MagneticField_cff")

## std sequence for PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

## std sequence to produce the ttFullHadEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff")
process.ttFullHadEvent.verbosity = 1

## choose which hypotheses to produce
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff import *
addTtFullHadHypotheses(process,
                       ["kKinFit"]
                       )

#removeTtFullHadHypGenMatch(process)

## change maximum number of jets taken into account per event (default: 6)
#from TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff import *
#setForAllTtFullHadHypotheses(process, "maxNJets", 8)

## process path
process.p = cms.Path(process.patDefaultSequence *
                     process.makeGenEvt         *
                     process.makeTtFullHadEvent
                     )

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName     = cms.untracked.string('ttFullHadEvtBuilder.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *'),                      
    dropMetaDataForDroppedData = cms.untracked.bool(True)
)
process.outpath = cms.EndPath(process.out)

## PAT content
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patEventContentNoCleaning
## TQAF content
from TopQuarkAnalysis.TopEventProducers.tqafEventContent_cff import *
process.out.outputCommands += tqafEventContent
