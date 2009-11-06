
import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIMRPRESCALE")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("FWCore.Modules.preScaler_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
#process.source = cms.Source("NewEventStreamFileReader",
    fileNames = cms.untracked.vstring(
        '/store/data/BeamCommissioning09/RandomTriggers/RAW/v1/000/119/805/443B111A-A0CA-DE11-848B-001617C3B77C.root'
#        '/store/data/Commissioning09/Test/RECO/v1/000/084/127/6CE9285A-A33A-DE11-8219-001D09F24600.root'
#       '/store/data/MWGR_21/Express/000/096/890/MWGR_21.00096890.0010.Express.storageManager.07.0000.dat'
#        '/store/data/MWGR_19/Express/000/084/127/MWGR_19.00084127.0082.Express.storageManager.03.0000.dat'
    )
)

# output module
#
process.load("Configuration.EventContent.EventContentCosmics_cff")



process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/DataOps/python/rawprescaleskimmer.py,v $'),
    annotation = cms.untracked.string('Test Skim Thingy')
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) ) ## default is false



process.prescale5 = process.preScaler.clone(
    prescaleFactor = 50,
    prescaleOffset = 0,
)
process.prescale1 = process.preScaler.clone(
    prescaleFactor = 10,
    prescaleOffset = 0,
)
process.prescale10 = process.preScaler.clone(
    prescaleFactor = 100,
    prescaleOffset = 0,
)
process.prescale50 = process.preScaler.clone(
    prescaleFactor = 500,
    prescaleOffset = 0,
)
process.prescale200 = process.preScaler.clone(
    prescaleFactor = 2000,
    prescaleOffset = 0,
)


process.fakeSkimOut1 = cms.OutputModule("PoolOutputModule",
    process.AODEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('prescalepath5')
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('PreScaleThingy50'),
        dataTier = cms.untracked.string('RAW')
    ),
    fileName = cms.untracked.string('PreScaleThingy5.root')
) 


process.fakeSkimOut2 = cms.OutputModule("PoolOutputModule",
    process.AODEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('prescalepath1')
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('PreScaleThingy10'),
        dataTier = cms.untracked.string('RAW')
    ),
    fileName = cms.untracked.string('PreScaleThingy1.root')
)

 
process.fakeSkimOut3 = cms.OutputModule("PoolOutputModule",
    process.AODEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('prescalepath10')
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('PreScaleThingy100'),
        dataTier = cms.untracked.string('RAW')
    ),
    fileName = cms.untracked.string('PreScaleThingy10.root')
)

 
process.fakeSkimOut4 = cms.OutputModule("PoolOutputModule",
    process.AODEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('prescalepath50')
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('PreScaleThingy50'),
        dataTier = cms.untracked.string('RAW')
    ),
    fileName = cms.untracked.string('PreScaleThingy50.root')
)

 
process.fakeSkimOut5 = cms.OutputModule("PoolOutputModule",
    process.FEVTEventContent,
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('prescalepath200')
    ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('PreScaleThingy2000'),
        dataTier = cms.untracked.string('RAW')
    ),
    fileName = cms.untracked.string('PreScaleThingy200.root')
)


process.endjob_step = cms.Path(process.endOfProcess)

process.prescalepath1 = cms.Path(process.prescale1)
process.prescalepath5 = cms.Path(process.prescale5)
process.prescalepath10 = cms.Path(process.prescale10)
process.prescalepath50 = cms.Path(process.prescale50)
process.prescalepath200 = cms.Path(process.prescale200)


process.e1=cms.EndPath(process.fakeSkimOut1)
process.e2=cms.EndPath(process.fakeSkimOut2)
process.e3=cms.EndPath(process.fakeSkimOut3)
process.e4=cms.EndPath(process.fakeSkimOut4)
process.e5=cms.EndPath(process.fakeSkimOut5)



#process.allPath = cms.Path( process.RawToDigi * process.reconstructionCosmics * process.DQMOfflineCosmics * process.MEtoEDMConverter * process.eventAuxiliaryHistoryProducer )

#avoid frequent exceptions in 226:
#process.hcalMonitor.TrigPrimMonitor = False

process.schedule = cms.Schedule( process.endjob_step,process.prescalepath1,process.prescalepath10,process.prescalepath5,process.prescalepath50,process.prescalepath200,process.e1,process.e2,process.e3,process.e4,process.e5 )

