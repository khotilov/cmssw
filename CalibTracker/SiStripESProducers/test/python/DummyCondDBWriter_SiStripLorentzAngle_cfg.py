# The following comments couldn't be translated into the new config version:

# upload to database 

#string timetype = "timestamp"    

import FWCore.ParameterSet.Config as cms

process = cms.Process("Builder")

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring(''),
    LorentzAngleBuilder = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('INFO')
    ),
    destinations = cms.untracked.vstring('LorentzAngleBuilder.log')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.source = cms.Source("EmptySource",
    numberEventsInRun = cms.untracked.uint32(1),
    firstRun = cms.untracked.uint32(1)
)

process.load("CalibTracker.SiStripESProducers.fake.SiStripLorentzAngleFakeESSource_cfi")
process.load("CalibTracker.SiStripESProducers.DBWriter.SiStripLorentzAngleDummyDBWriter_cfi")

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(2),
        authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')
    ),
    timetype = cms.untracked.string('runnumber'),
<<<<<<< DummyCondDBWriter_SiStripLorentzAngle_cfg.py
    connect = cms.string('sqlite_file:LA_test_dummy.db'),
=======
    connect = cms.string('sqlite_file:dbfile.db'),
>>>>>>> 1.2
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('SiStripLorentzAngleRcd'),
        tag = cms.string('SiStripLorentzAngle_Fake_30X')
    ))
)

process.SiStripLorentzAngleGenerator.TIB_EstimatedValue = cms.double(0.0174)
process.SiStripLorentzAngleGenerator.TOB_EstimatedValue = cms.double(0.0222)
process.SiStripLorentzAngleGenerator.PerCent_Err = cms.double(20)

process.p1 = cms.Path(process.siStripLorentzAngleDummyDBWriter)


