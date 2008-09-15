import FWCore.ParameterSet.Config as cms

process = cms.Process("PixelDBReader")
process.load("Geometry.TrackerSimData.trackerSimGeometryXML_cfi")

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("CondTools.SiPixel.SiPixelGainCalibrationService_cfi")

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource",
    numberEventsInRun = cms.untracked.uint32(10),
    firstRun = cms.untracked.uint32(1)
)

process.Timing = cms.Service("Timing")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
    ignoreTotal = cms.untracked.int32(0)
)

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
    process.CondDBCommon,
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelGainCalibrationForHLTRcd'),
        tag = cms.string('V2_trivial_TBuffer_hlt')
    ))
)

process.prefer("PoolDBESSource")
process.SiPixelCondObjForHLTReader = cms.EDFilter("SiPixelCondObjForHLTReader",
    process.SiPixelGainCalibrationServiceParameters,
    fileName = cms.string('histos_HLT.root')
)

#process.print = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.SiPixelCondObjForHLTReader)
#process.ep = cms.EndPath(process.print)
process.CondDBCommon.connect = 'sqlite_file:prova.db'
process.CondDBCommon.DBParameters.messageLevel = 2
process.CondDBCommon.DBParameters.authenticationPath = ''


