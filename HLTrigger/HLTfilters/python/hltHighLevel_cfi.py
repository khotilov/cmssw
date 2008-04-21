import FWCore.ParameterSet.Config as cms

hltHighLevel = cms.EDFilter("HLTHighLevel",
    HLTPaths = cms.vstring(), ## provide list of HLT paths you want

    andOr = cms.bool(True), ## true: OR of those on your list 

    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
)


