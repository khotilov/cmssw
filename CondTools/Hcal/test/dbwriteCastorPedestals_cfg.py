import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.MessageLogger=cms.Service("MessageLogger",
                              destinations=cms.untracked.vstring("cout"),
                              cout=cms.untracked.PSet(
                              treshold=cms.untracked.string("INFO")
                              )
)

process.load("CondCore.DBCommon.CondDBCommon_cfi")

#process.CondDBCommon.connect = cms.string('oracle://cms_orcoff_prep/CMS_COND_30X_HCAL')
#process.CondDBCommon.DBParameters.authenticationPath = cms.untracked.string('./authentication.xml')
process.CondDBCommon.connect = cms.string('sqlite_file:testExample.db')
process.CondDBCommon.DBParameters.authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')

process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    interval = cms.uint64(1)
)

process.es_ascii = cms.ESSource("CastorTextCalibrations",
    input = cms.VPSet(cms.PSet(
        object = cms.string('Pedestals'),
        file = cms.FileInPath('CondFormats/CastorObjects/data/castor_pedestals_test.txt')
    ))
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    timetype = cms.untracked.string('runnumber'),
    logconnect= cms.untracked.string('sqlite_file:log.db'),
#    logconnect= cms.untracked.string('oracle://cms_orcoff_prep/CMS_COND_30X_POPCONLOG'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('CastorPedestalsRcd'),
        tag = cms.string('castor_pedestals_v1.0_test')
         ))
)

process.mytest = cms.EDAnalyzer("CastorPedestalsPopConAnalyzer",
    record = cms.string('CastorPedestalsRcd'),
    loggingOn= cms.untracked.bool(True),
    SinceAppendMode=cms.bool(True),
    Source=cms.PSet(
#    firstSince=cms.untracked.double(300) 
    IOVRun=cms.untracked.uint32(1)
    )
)

process.p = cms.Path(process.mytest)
