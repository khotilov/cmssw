import FWCore.ParameterSet.Config as cms

ecalEndcapRecHitsValidation = cms.EDFilter("EcalEndcapRecHitsValidation",
    EEdigiCollection = cms.InputTag("simEcalDigis","eeDigis"),
    EEuncalibrechitCollection = cms.InputTag("ecalWeightUncalibRecHit","EcalUncalibRecHitsEE"),
    verbose = cms.untracked.bool(True)
)



