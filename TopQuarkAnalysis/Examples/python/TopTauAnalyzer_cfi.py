import FWCore.ParameterSet.Config as cms

#
# module to make simple analyses of tautrons
#
analyzeTau = cms.EDAnalyzer("TopTauAnalyzer",
    input = cms.InputTag("selectedLayer1Taus")
)


