import FWCore.ParameterSet.Config as cms

process = cms.Process("GeometryTest")
process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.load("DetectorDescription.OfflineDBLoader.test.cmsIdealGeometryForWrite_cfi")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptyIOVSource",
                                lastRun = cms.untracked.uint32(1),
                                timetype = cms.string('runnumber'),
                                firstRun = cms.untracked.uint32(1),
                                interval = cms.uint32(1)
                            )

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                              process.CondDBCommon,
                                              toPut = cms.VPSet(cms.PSet(
            record = cms.string('IdealGeometryRecord'),
                    tag = cms.string('IdealGeometry01')
                ))
                                          )

process.load = cms.EDFilter("WriteOneGeometryFromXML",
                                rotNumSeed = cms.int32(0),
                                dumpSpecs = cms.untracked.bool(False),
                                dumpGeoHistory = cms.untracked.bool(False),
                                dumpPosInfo = cms.untracked.bool(False)
                            )

process.myprint = cms.OutputModule("AsciiOutputModule")

process.p1 = cms.Path(process.load)
process.ep = cms.EndPath(process.myprint)
process.CondDBCommon.connect = 'sqlite_file:testIdeal.db'
process.CondDBCommon.DBParameters.messageLevel = 0
process.CondDBCommon.DBParameters.authenticationPath = './'
