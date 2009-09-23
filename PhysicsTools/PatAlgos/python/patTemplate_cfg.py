import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_1_2/RelValTTbar/GEN-SIM-RECO/STARTUP31X_V2-v1/0006/0CC00E3B-5A78-DE11-A2AB-000423D94A04.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP31X_V1::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

# Output module configuration
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                                   # save only events passing the full path
                                   SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                                   # save PAT Layer 1 output
                                   outputCommands = cms.untracked.vstring('drop *', *patEventContent ) # you need a '*' to unpack the list of commands 'patEventContent'
                               )
process.outpath = cms.EndPath(process.out)
