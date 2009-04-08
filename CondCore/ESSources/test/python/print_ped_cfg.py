import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0),
        authenticationPath = cms.untracked.string('.')
    ),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('PedestalsRcd'),
        tag = cms.string('mytest')
    )),
    connect = cms.string('sqlite_file:test.db')
)
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(16),
    numberEventsInRun = cms.untracked.uint32(1)
)
process.maxEvents=cms.untracked.PSet(input=cms.untracked.int32(5))
process.prod = cms.EDAnalyzer("PedestalsAnalyzer")

process.p = cms.Path(process.prod)


