import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.MessageLogger=cms.Service("MessageLogger",
                                  cout=cms.untracked.PSet(threshold=cms.untracked.string('INFO')),
                                  destinations=cms.untracked.vstring("cout")
                                  )

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.CondDBCommon.connect = cms.string('oracle://cms_orcoff_prep/CMS_COND_31X_DQM_SUMMARY')
process.CondDBCommon.DBParameters.authenticationPath = cms.untracked.string('/build/diguida/conddb')
process.CondDBCommon.BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
process.CondDBCommon.DBParameters.messageLevel = cms.untracked.int32(1) #3 for high verbosity

process.source = cms.Source("EmptyIOVSource", #needed to EvSetup in order to load data
                            timetype = cms.string('runnumber'),
                            firstValue = cms.uint64(1),
                            lastValue = cms.uint64(1),
                            interval = cms.uint64(1)
                            )

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
                                          process.CondDBCommon,
                                          timetype = cms.untracked.string('runnumber'),
                                          toPut = cms.VPSet(cms.PSet(record = cms.string('GeometryFile'),
                                                                     tag = cms.string('ROOTFILE_Test') 
                                                                     )
                                                            ),
                                          logconnect = cms.untracked.string('sqlite_file:ROOTFILE_TestLog.db')                                     
                                          )

process.dqmReferenceHistogramRootFileTest = cms.EDAnalyzer("DQMReferenceHistogramRootFilePopConAnalyzer",
                                           record = cms.string('GeometryFile'),
                                           loggingOn = cms.untracked.bool(True), #always True, needs to create the log db
                                           SinceAppendMode = cms.bool(True),
                                           Source = cms.PSet(ROOTFile = cms.untracked.string("salvo.root"),
                                                             firstSince = cms.untracked.uint64(1), #1, 43434, 46335, 51493, 51500
                                                             debug = cms.untracked.bool(True)
                                                             )
                                           )

process.p = cms.Path(process.dqmReferenceHistogramRootFileTest)
