import FWCore.ParameterSet.Config as cms

process = cms.Process("TQAF")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_8_2/RelValTTbar/GEN-SIM-RECO/MC_38Y_V9-v1/0018/E8B5D618-96AF-DF11-835A-003048679070.root'
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

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['mc']

#-------------------------------------------------
# TQAF configuration
# (comment the patDefaultSequence from the path if
# you want to produce TQAF on top of already
# existing PAT objects)
#-------------------------------------------------

## std sequence for PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## std sequence for TQAF
process.load("TopQuarkAnalysis.TopEventProducers.tqafSequences_cff")

## switch jet collection
#from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJetCollection(process, 
#                    cms.InputTag('sisCone5CaloJets'),
#                    doJTA            = True,           
#                    doBTagging       = True,           
#                    jetCorrLabel     = ('SC5','Calo'),
#                    doType1MET       = True,          
#                    genJetCollection = cms.InputTag("sisCone5GenJets"),
#                    doJetID          = True,
#                    jetIdLabel       = "sc5"
#                    )

## process path
process.p = cms.Path(process.patDefaultSequence *
                     process.tqafTtSemiLeptonic
                     )

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    fileName       = cms.untracked.string('tqafOutput.root'),
    SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('drop *'),                      
    dropMetaData   = cms.untracked.string("DROPPED")  ## NONE    for none
                                                      ## DROPPED for drop for dropped data
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
