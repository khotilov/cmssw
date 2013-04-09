import FWCore.ParameterSet.Config as cms

process = cms.Process("CLEANING")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'file:pickevents.root'
)
)


process.load('MuonAnalysis.PFMuonTools.muonCleaning_cff')

#redo MET
from RecoMET.METProducers.METSigParams_cfi import *


process.pfMetCleaned = cms.EDProducer("METProducer",
                       METSignificance_params,
                       src = cms.InputTag("particleFlowCleaned"),
                       METType = cms.string('PFMET'),
                       alias = cms.string('PFMET'),
                       noHF = cms.bool(False),
                       globalThreshold = cms.double(0.0),
                       InputType = cms.string('PFCandidateCollection'),
                       calculateSignificance = cms.bool(False),
                       jets = cms.InputTag("ak5PFJets") #used for significance calculation
)                       

#make an old PF Muon collection to compare
process.oldPFMuons = cms.EDFilter("PdgIdPFCandidateSelector",
                       src = cms.InputTag("particleFlow"),
                       pdgId = cms.vint32( -13, 13)
)


process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('metOutputPE.root'),
)
             


process.p = cms.Path(process.muonCleaning+process.pfMetCleaned+process.oldPFMuons)
process.e = cms.EndPath(process.output)

