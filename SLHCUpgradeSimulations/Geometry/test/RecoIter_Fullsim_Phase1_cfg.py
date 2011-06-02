# Auto generated configuration file
# using: 
# Revision: 1.172.2.5 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 -s RECO -n 100 --conditions DESIGN_36_V10::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --beamspot Gauss --fileout file:reco.root --filein file:raw.root --python_filename RecoMuon_Fullsim_cfg.py --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("SLHCUpgradeSimulations.Geometry.mixLowLumPU_Phase1_R39F16_cff")
process.load("SLHCUpgradeSimulations.Geometry.Phase1_R39F16_cmsSimIdealGeometryXML_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('SLHCUpgradeSimulations.Geometry.Digi_Phase1_R39F16_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.13 $'),
    annotation = cms.untracked.string('step2 nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)

)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValFourMuons/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/EC71F6CD-8980-E011-AEF0-001A92971B64.root'
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/F4C5F6C6-8980-E011-ACBA-003048678FD6.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/E053363C-F27F-E011-9A0B-003048678B92.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/B41CEE30-0880-E011-B141-003048D3C010.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/A82A6871-0F80-E011-A989-0018F3D09630.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/A44BD2FF-0380-E011-A9C3-003048678FF8.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/6E6218B3-F27F-E011-BDDC-002354EF3BE0.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/58F46179-0280-E011-A314-002618943925.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/50390620-F37F-E011-AF35-00304867BFBC.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/36D18CAA-F57F-E011-9F80-0018F3D09628.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/2EDE4C37-F47F-E011-8896-003048679150.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/1A281128-F57F-E011-B2C8-002354EF3BDA.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/12885030-F67F-E011-A6D1-0018F3D0970E.root',
	#'/store/relval/CMSSW_4_2_3_SLHC_pre1/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_20110605_special-v1/0026/0608F219-F77F-E011-A026-003048D3FC94.root' 
    )
)
# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:valid_reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)
#I'm only interested in the validation stuff
process.output.outputCommands = cms.untracked.vstring('drop *','keep *_MEtoEDMConverter_*_*')

#process.output = cms.OutputModule("PoolOutputModule",
#         outputCommands = process.AODSIMEventContent.outputCommands,
#         fileName = cms.untracked.string(
#		'file:/uscms_data/d2/brownson/slhc/quadMuon_RECO.root')
#)


# Additional output definition

# Other statements
#process.GlobalTag.globaltag = 'MC_42_V10::All'
process.GlobalTag.globaltag = 'DESIGN42_V11::All'

### PhaseI Geometry and modifications ###############################################
process.Timing =  cms.Service("Timing")
## no playback when doing digis
#process.mix.playback = True
#process.MessageLogger.destinations = cms.untracked.vstring("detailedInfo_fullph1geom")

### if pileup we need to set the number
process.mix.input.nbPileupEvents = cms.PSet(
  averageNumber = cms.double(50.0)
)
### if doing inefficiency at <PU>=50
process.simSiPixelDigis.AddPixelInefficiency = 20
## also for strips TIB inefficiency if we want
## TIB1,2 inefficiency at 20%
#process.simSiStripDigis.Inefficiency = 20
## TIB1,2 inefficiency at 50%
#process.simSiStripDigis.Inefficiency = 30
## TIB1,2 inefficiency at 99% (i.e. dead)
#process.simSiStripDigis.Inefficiency = 40

process.load("SLHCUpgradeSimulations.Geometry.fakeConditions_Phase1_R39F16_cff")
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_phase1_cff")

process.ctfWithMaterialTracks.TTRHBuilder = 'WithTrackAngle'
process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(False)
process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(False)
process.PixelCPEGenericESProducer.Upgrade = cms.bool(True)
process.PixelCPEGenericESProducer.SmallPitch = False
process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
process.PixelCPEGenericESProducer.DoCosmics = False

## CPE for other steps
process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.newPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.secPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.thPixelRecHits.CPE = cms.string('PixelCPEGeneric')
process.preFilterZeroStepTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.preFilterStepOneTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.secWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.thWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fourthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')
process.fifthWithMaterialTracks.TTRHBuilder = cms.string('WithTrackAngle')

# Need these lines to stop some errors about missing siStripDigis collections.
# should add them to fakeConditions_Phase1_cff
process.MeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.MeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.MeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.newMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.newMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.newMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.newMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.newMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.secMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.secMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.secMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.secMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.secMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.thMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.thMeasurementTracker.UseStripModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UseStripAPVFiberQualityDB   = cms.bool(False)
process.thMeasurementTracker.UseStripStripQualityDB      = cms.bool(False)
process.thMeasurementTracker.UsePixelModuleQualityDB     = cms.bool(False)
process.thMeasurementTracker.UsePixelROCQualityDB        = cms.bool(False)
process.fourthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()
process.fifthMeasurementTracker.inactiveStripDetectorLabels = cms.VInputTag()

### Now Validation and other user functions #########################################
process.load("Validation.RecoTrack.cutsTPEffic_cfi")
process.load("Validation.RecoTrack.cutsTPFake_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load('SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi')
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.load('Configuration.StandardSequences.Validation_cff')
### look look at OOTB generalTracks and high purity collections
### for high purity also look at 6 and 8 hit requirements
### some definitions in Validation/RecoTrack/python/TrackValidation_cff.py

import PhysicsTools.RecoAlgos.recoTrackSelector_cfi

process.cutsRecoTracksHpw6hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw6hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw6hits.minHit=cms.int32(6)

process.cutsRecoTracksHpw8hits = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpw8hits.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpw8hits.minHit=cms.int32(8)

process.cutsRecoTracksHpwbtagc = PhysicsTools.RecoAlgos.recoTrackSelector_cfi.recoTrackSelector.clone()
process.cutsRecoTracksHpwbtagc.quality=cms.vstring("highPurity")
process.cutsRecoTracksHpwbtagc.minHit=cms.int32(8)
process.cutsRecoTracksHpwbtagc.ptMin = cms.double(1.0)

process.trackValidator.label=cms.VInputTag(cms.InputTag("generalTracks"),
                                           cms.InputTag("cutsRecoTracksHp"),
                                           cms.InputTag("cutsRecoTracksHpwbtagc"),
                                           cms.InputTag("cutsRecoTracksZeroHp"),
                                           cms.InputTag("cutsRecoTracksFirstHp")
#                                           cms.InputTag("cutsRecoTracksSecondHp"),
#                                           cms.InputTag("cutsRecoTracksThirdHp")
                                           )
#process.trackValidator.associators = ['TrackAssociatorByHits']
process.trackValidator.associators = cms.vstring('quickTrackAssociatorByHits')
process.trackValidator.UseAssociators = True
## options to match with 363 histos for comparison
process.trackValidator.histoProducerAlgoBlock.nintEta = cms.int32(20)
process.trackValidator.histoProducerAlgoBlock.nintPt = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.maxPt = cms.double(200.0)
process.trackValidator.histoProducerAlgoBlock.useLogPt = cms.untracked.bool(True)
process.trackValidator.histoProducerAlgoBlock.minDxy = cms.double(-3.0)
process.trackValidator.histoProducerAlgoBlock.maxDxy = cms.double(3.0)
process.trackValidator.histoProducerAlgoBlock.nintDxy = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.minDz = cms.double(-10.0)
process.trackValidator.histoProducerAlgoBlock.maxDz = cms.double(10.0)
process.trackValidator.histoProducerAlgoBlock.nintDz = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.maxVertpos = cms.double(5.0)
process.trackValidator.histoProducerAlgoBlock.nintVertpos = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.minZpos = cms.double(-10.0)
process.trackValidator.histoProducerAlgoBlock.maxZpos = cms.double(10.0)
process.trackValidator.histoProducerAlgoBlock.nintZpos = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.phiRes_rangeMin = cms.double(-0.003)
process.trackValidator.histoProducerAlgoBlock.phiRes_rangeMax = cms.double(0.003)
process.trackValidator.histoProducerAlgoBlock.phiRes_nbin = cms.int32(100)
process.trackValidator.histoProducerAlgoBlock.cotThetaRes_rangeMin = cms.double(-0.01)
process.trackValidator.histoProducerAlgoBlock.cotThetaRes_rangeMax = cms.double(+0.01)
process.trackValidator.histoProducerAlgoBlock.cotThetaRes_nbin = cms.int32(120)
process.trackValidator.histoProducerAlgoBlock.dxyRes_rangeMin = cms.double(-0.01)
process.trackValidator.histoProducerAlgoBlock.dxyRes_rangeMax = cms.double(0.01)
process.trackValidator.histoProducerAlgoBlock.dxyRes_nbin = cms.int32(100)

process.slhcTracksValidation = cms.Sequence(process.cutsRecoTracksHp*
                                 process.cutsRecoTracksHpwbtagc*
                                 process.cutsRecoTracksZeroHp*
                                 process.cutsRecoTracksFirstHp*
#                                 process.cutsRecoTracksSecondHp*
#                                 process.cutsRecoTracksThirdHp*
                                 process.trackValidator)

############ end John's changes ###########################
process.ReadLocalMeasurement = cms.EDAnalyzer("StdHitNtuplizer",
   src = cms.InputTag("siPixelRecHits"),
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   ### if using simple (non-iterative) or old (as in 1_8_4) tracking
   trackProducer = cms.InputTag("generalTracks"),
   OutputFile = cms.string("stdgrechitfullph1g_ntuple.root"),
   ### for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('g4SimHitsTrackerHitsPixelBarrelLowTof',
                         'g4SimHitsTrackerHitsPixelBarrelHighTof',
                         'g4SimHitsTrackerHitsPixelEndcapLowTof',
                         'g4SimHitsTrackerHitsPixelEndcapHighTof')
)
process.anal = cms.EDAnalyzer("EventContentAnalyzer")

## need this at the end as the validation config redefines random seed with just mix
#process.load("IOMC.RandomEngine.IOMC_cff")

### back to standard job commands ##################################################
process.pdigi.remove(process.simCastorDigis)

# Path and EndPath definitions
process.mix_step 		= cms.Path(process.mix)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)

process.reconstruction_step 	= cms.Path(process.trackerlocalreco*
						process.offlineBeamSpot+
                                                process.recopixelvertexing*process.ckftracks_wodEdXandSteps2345)
#                                                process.recopixelvertexing*process.ckftracks_wodEdXandSteps4and5)
process.debug_step 		= cms.Path(process.anal)
process.validation_step 	= cms.Path(process.cutsTPEffic*
						process.cutsTPFake*
						process.slhcTracksValidation)
process.user_step 		= cms.Path(process.ReadLocalMeasurement)
process.endjob_step 		= cms.Path(process.endOfProcess)
process.out_step 		= cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.validation_step,process.user_step,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.reconstruction_step,process.validation_step,process.user_step,process.endjob_step,process.out_step)
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.out_step)

