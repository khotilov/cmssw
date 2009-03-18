# cfg file to test the online producer of L1GtTriggerMenuRcd

import FWCore.ParameterSet.Config as cms

process = cms.Process("L1ConfigWritePayloadDummy")

# number of events and source
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    interval = cms.uint64(1)
)

# TSC key
process.load("CondTools.L1Trigger.L1SubsystemKeysOnline_cfi")
process.L1SubsystemKeysOnline.tscKey = cms.string( 'TSC_000618_090304_MIDWEEK2008_GTgt20090bst30_GMTstartup3_DTTF_DT_MI' )

# Subclass of L1ObjectKeysOnlineProdBase.
process.load("L1TriggerConfig.L1GtConfigProducers.l1GtTscObjectKeysOnline_cfi")
process.l1GtTscObjectKeysOnline.subsystemLabel = cms.string('')

# Generate dummy L1TriggerKeyList
process.load("CondTools.L1Trigger.L1TriggerKeyListDummy_cff")

# Get configuration data from OMDS.  This is the subclass of L1ConfigOnlineProdBase.
process.load("L1TriggerConfig.L1GtConfigProducers.l1GtTriggerMenuOnline_cfi")


process.getter = cms.EDAnalyzer("EventSetupRecordDataGetter",
   toGet = cms.VPSet(cms.PSet(
   record = cms.string('L1GtTriggerMenuRcd'),
   data = cms.vstring('L1GtTriggerMenu')
   )),
   verbose = cms.untracked.bool(True)
)

process.p = cms.Path(process.getter)

# Message Logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
process.MessageLogger.cout.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.debugModules = cms.untracked.vstring('*')

