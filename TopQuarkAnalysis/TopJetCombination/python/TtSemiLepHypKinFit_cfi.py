import FWCore.ParameterSet.Config as cms

#
# module to make the kinematic fit hypothesis
#
ttSemiLepHypKinFit = cms.EDProducer("TtSemiLepHypKinFit",
    jets  = cms.InputTag("selectedLayer1Jets"),
    leps  = cms.InputTag("selectedLayer1Muons"),
    mets  = cms.InputTag("selectedLayer1METs"),
    match = cms.InputTag("kinFitTtSemiLepEventHypothesis"),
    status    = cms.InputTag("kinFitTtSemiLepEventHypothesis","Status"),
    partons   = cms.InputTag("kinFitTtSemiLepEventHypothesis","Partons"),
    leptons   = cms.InputTag("kinFitTtSemiLepEventHypothesis","Leptons"),
    neutrinos = cms.InputTag("kinFitTtSemiLepEventHypothesis","Neutrinos")
)


