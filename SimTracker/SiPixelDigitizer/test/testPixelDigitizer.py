# python script Pixel Digitizer testing:
# author: V. Cuplov
#
# Process available:
# 1) simulation of hits in the Tracker
# 2) digitization + validation
# 3) reconstruction or clusterization + validation 
# 4) tracks + validation
#
##############################################################################

import FWCore.ParameterSet.Config as cms

process = cms.Process("DigiTest")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Validation.TrackerHits.trackerHitsValidation_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimTracker.Configuration.SimTracker_cff")
process.load("Validation.TrackerDigis.trackerDigisValidation_cff")
process.load("SimG4Core.Configuration.SimG4Core_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

## For 31X use:
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
## for older releases
#process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")

process.load("Validation.TrackerRecHits.trackerRecHitsValidation_cff")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("Validation.TrackingMCTruth.trackingTruthValidation_cfi")
process.load("Validation.RecoTrack.TrackValidation_cff")
process.load("Validation.RecoTrack.SiTrackingRecHitsValid_cff")
process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")
process.load("EventFilter.SiPixelRawToDigi.SiPixelDigiToRaw_cfi")
process.load("EventFilter.SiPixelRawToDigi.SiPixelRawToDigi_cfi")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50)
)

#process.MessageLogger = cms.Service("MessageLogger",
#    debugModules = cms.untracked.vstring('PixelDigisTest'),
#    destinations = cms.untracked.vstring('cout'),
#    destinations = cms.untracked.vstring("log","cout"),
#    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('ERROR')
#    )
#    log = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#    )
#)


# get the files from DBS:
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/B41EF45B-D14E-DE11-BC68-001D09F24600.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/682FF80C-524F-DE11-9159-001D09F2543D.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/0E9B84FD-D24E-DE11-94D9-001617C3B70E.root')
#)

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource", fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( ['/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/B41EF45B-D14E-DE11-BC68-001D09F24600.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/682FF80C-524F-DE11-9159-001D09F2543D.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-RECO/IDEAL_31X_v1/0007/0E9B84FD-D24E-DE11-94D9-001617C3B70E.root'] );

secFiles.extend( ['/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/B479D4A7-D14E-DE11-A78F-001D09F250AF.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/A6F81696-D04E-DE11-8B2E-001617C3B6DC.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/A625B90A-D34E-DE11-A4D3-000423D6A6F4.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/5441A448-D14E-DE11-B5BC-001D09F295A1.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/4C9657A5-D34E-DE11-9C9D-001D09F248F8.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/3EA6B3D7-514F-DE11-8DFD-000423D6A6F4.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/2C49163C-D24E-DE11-AA16-001D09F2438A.root','/store/relval/CMSSW_3_1_0_pre9/RelValSingleMuPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_31X_v1/0007/1E69E467-D04E-DE11-8A75-001617E30D4A.root'] )


#process.PoolDBESSource = cms.ESSource("PoolDBESSource",
#    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
#    DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0),
#        authenticationPath = cms.untracked.string('')
#    ),
#    timetype = cms.string('runnumber'),
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('SiPixelQualityRcd'),
#        tag = cms.string('SiPixelBadModule_test')
#    )),
#    connect = cms.string('sqlite_file:test.db')
#)


# Choose the global tag here:
process.GlobalTag.globaltag = 'MC_31X_V1::All'

process.o1 = cms.OutputModule("PoolOutputModule",
                              outputCommands = cms.untracked.vstring('drop *','keep *_*_*_DigiTest'),
            fileName = cms.untracked.string('file:/afs/cern.ch/user/v/vesna/DigitizerWork/CMSSW_3_2_1/src/SimTracker/SiPixelDigitizer/test/Digis.root')  
)

process.Timing = cms.Service("Timing")

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

# Simulation of hits in the Tracker:
#process.simhits = cms.Sequence(process.g4SimHits*process.trackerHitsValidation)

# Digitization:
process.digis = cms.Sequence(process.simSiPixelDigis*process.pixelDigisValid)

# Reconstruction or Clusterization:
#process.rechits = cms.Sequence(process.pixeltrackerlocalreco*process.pixRecHitsValid)
process.rechits = cms.Sequence(process.pixeltrackerlocalreco*process.siStripMatchedRecHits*process.pixRecHitsValid)

#process.simSiPixelDigis.ChargeVCALSmearing=True
process.tracks = cms.Sequence(process.offlineBeamSpot*process.recopixelvertexing*process.trackingParticles*process.trackingTruthValid*process.ckftracks*process.trackerRecHitsValidation)
process.trackinghits = cms.Sequence(process.TrackRefitter*process.trackingRecHitsValid)

# To use when you want to look at Residuals and Pulls:
#process.p1 = cms.Path(process.mix*process.digis*process.siPixelRawData*process.siPixelDigis*process.rechits*process.tracks*process.trackinghits)

#process.p1 = cms.Path(process.mix*process.digis*process.siPixelRawData*process.siPixelDigis*process.pixeltrackerlocalreco)

#This process to get events to be used to make DQM occupancy maps (cmsRun DQM_Pixel_digi.py after you ran testPixelDigitizer.py):
process.p1 = cms.Path(process.mix*process.simSiPixelDigis*process.siPixelRawData*process.siPixelDigis)

# Look at cluster charge:
#process.p1 = cms.Path(process.mix*process.digis*process.siPixelRawData*process.siPixelDigis*process.pixeltrackerlocalreco)

#This process is to run the digis validation:
#process.p1 = cms.Path(process.mix*process.simSiPixelDigis*process.pixelDigisValid)

#This process is to run the digitizer:
#process.p1 = cms.Path(process.mix*process.simSiPixelDigis)

#process.p1 = cms.Path(process.mix*process.digis*process.siPixelRawData*process.siPixelDigis)

process.outpath = cms.EndPath(process.o1)
process.g4SimHits.Generator.HepMCProductLabel = 'source'

