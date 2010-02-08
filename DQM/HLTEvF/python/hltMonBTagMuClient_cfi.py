import FWCore.ParameterSet.Config as cms

hltMonBTagMuClient = cms.EDFilter('HLTMonBTagClient',
    monitorName     = cms.string('HLT/HLTMonBJet'),
    pathName        = cms.string('HLT_BTagMu_Jet10U'),
    storeROOT       = cms.untracked.bool(False),
    outputFile      = cms.untracked.string('HLTMonBTag.root')
)
