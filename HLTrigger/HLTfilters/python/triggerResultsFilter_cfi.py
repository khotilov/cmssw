import FWCore.ParameterSet.Config as cms

triggerResultsFilter = cms.EDFilter('TriggerResultsFilter',
    hltResults    = cms.InputTag('TriggerResults'),     # HLT results   - set to empty to ignore HLT
    l1tResults    = cms.InputTag('hltGtDigis'),         # L1 GT results - set to empty to ignore L1
    l1tIgnoreMask = cms.bool(False),                    # use L1 mask
    daqPartitions = cms.uint32(0x01),                   # used by the definition of the L1 mask
    throw         = cms.bool(True),                     # throw exception on unknown trigger names
    triggerConditions = cms.string('HLT_*')
)
