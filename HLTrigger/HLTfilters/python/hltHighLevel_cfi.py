import FWCore.ParameterSet.Config as cms

hltHighLevel = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
    HLTPaths = cms.vstring(),           # provide list of HLT paths (or patterns) you want
    andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.untracked.bool(True)    # throw exception on unknown path names
)
