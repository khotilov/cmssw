import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = 'oracle://devdb10/cms_xiezhen_dev'
process.CondDBCommon.DBParameters.authenticationPath = '/afs/cern.ch/user/x/xiezhen'

process.source = cms.Source("EmptyIOVSource",
    lastValue = cms.uint64(1),
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    interval = cms.uint64(1)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('mySiStripNoisesRcd'),
        tag = cms.string('noise_tag')
    ))
)

process.mytest = cms.EDAnalyzer("writeBlob")

process.p = cms.Path(process.mytest)


