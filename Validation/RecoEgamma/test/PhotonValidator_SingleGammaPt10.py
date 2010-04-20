import FWCore.ParameterSet.Config as cms

process = cms.Process("TestPhotonValidator")
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.load("Validation.RecoEgamma.photonValidationSequence_cff")
process.load("Validation.RecoEgamma.photonPostprocessing_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_3XY_V26::All'

process.DQMStore = cms.Service("DQMStore");
process.load("DQMServices.Components.DQMStoreStats_cfi")
from DQMServices.Components.DQMStoreStats_cfi import *
dqmStoreStats.runOnEndJob = cms.untracked.bool(True)



process.maxEvents = cms.untracked.PSet(
#input = cms.untracked.int32(10)
)



from Validation.RecoEgamma.photonValidationSequence_cff import *
from Validation.RecoEgamma.photonPostprocessing_cfi import *

photonValidation.OutputMEsInRootFile = True
photonValidation.OutputFileName = 'PhotonValidationRelVal357_SingleGammaPt10.root'

photonPostprocessing.batch = cms.bool(True)
photonPostprocessing.InputFileName = photonValidation.OutputFileName

process.source = cms.Source("PoolSource",
noEventSort = cms.untracked.bool(True),
duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            
    fileNames = cms.untracked.vstring(
# official RelVal 357 single Photons pt=10GeV
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0013/EA921D90-6949-DF11-A25A-002618943868.root',
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0012/E610BBA2-3F49-DF11-BD52-003048678B44.root'
    ),
                            
    secondaryFileNames = cms.untracked.vstring(
# official RelVal 357 single Photons pt=10GeV    
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_3XY_V26-v1/0013/BCFFC10C-6949-DF11-BBAA-0030486790B8.root',
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_3XY_V26-v1/0012/D646041A-4049-DF11-A0CC-001A92971BB4.root',
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_3XY_V26-v1/0012/BC26C6A1-3F49-DF11-AD5B-003048D15CC0.root',
        '/store/relval/CMSSW_3_5_7/RelValSingleGammaPt10/GEN-SIM-DIGI-RAW-HLTDEBUG/MC_3XY_V26-v1/0012/9AF01406-3F49-DF11-A534-0030486790FE.root' 
    )
 )


photonPostprocessing.rBin = 48

## For single gamma pt =10
photonValidation.eMax  = 100
photonValidation.etMax = 50
photonValidation.etScale = 0.20
photonPostprocessing.eMax  = 100
photonPostprocessing.etMax = 50



process.FEVT = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring("keep *_MEtoEDMConverter_*_*"),
    fileName = cms.untracked.string('pippo.root')
)



process.p1 = cms.Path(process.tpSelection*process.photonValidationSequence*process.photonPostprocessing*process.dqmStoreStats)
process.schedule = cms.Schedule(process.p1)



