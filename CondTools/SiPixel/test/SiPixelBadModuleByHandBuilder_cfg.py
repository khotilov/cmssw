import FWCore.ParameterSet.Config as cms

process = cms.Process("ICALIB")
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.source = cms.Source("EmptyIOVSource",
    lastValue = cms.uint64(84000),
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(84000),
    interval = cms.uint64(1)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('')
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:Quality_v03.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityRcd'),
        tag = cms.string('SiPixelQuality_v03')
    ))
)

process.prod = cms.EDFilter("SiPixelBadModuleByHandBuilder",
    BadModuleList = cms.untracked.VPSet(cms.PSet(
        errortype = cms.string('whole'),
        detid = cms.uint32(302197784)
    ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(302195232)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(302123296)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(302127136)
        ), 
        cms.PSet(
            errortype = cms.string('tbmA'),
            detid = cms.uint32(302125076)
        ), 
        cms.PSet(
            errortype = cms.string('tbmB'),
            detid = cms.uint32(302126364)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(302188552)
        ), 
        cms.PSet(
            errortype = cms.string('tbmA'),
            detid = cms.uint32(302121992)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(302126596)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344014340)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344014344)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344014348)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344019460)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344019464)
        ), 
        cms.PSet(
            errortype = cms.string('whole'),
            detid = cms.uint32(344019468)
        ), 
        cms.PSet(
            errortype = cms.string('none'),
            detid = cms.uint32(302187268),
            badroclist = cms.vuint32(6)
        ),
        cms.PSet(
            errortype = cms.string('none'),
            detid = cms.uint32(302195472),
            badroclist = cms.vuint32(0)
        ),
        cms.PSet(
            errortype = cms.string('none'),
            detid = cms.uint32(302128136),
            badroclist = cms.vuint32(3)
         )),
    Record = cms.string('SiPixelQualityRcd'),
    SinceAppendMode = cms.bool(True),
    IOVMode = cms.string('Run'),
    printDebug = cms.untracked.bool(True),
    doStoreOnDB = cms.bool(True)
)

#process.print = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.prod)
#process.ep = cms.EndPath(process.print)


