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
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.8 $'),
    annotation = cms.untracked.string('step2 nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)

)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/FAE740C0-698F-E011-B68F-0030487A1884.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/F83F7254-D68E-E011-BDC1-0030487A1FEC.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/D8DF747F-D38E-E011-B1FB-0030487CAEAC.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/A4E7FE3A-D88E-E011-B82A-0030487CD716.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/A2EBE7A8-D48E-E011-B2BF-0030487CAEAC.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/76F3D653-D78E-E011-8E7C-0030487CD718.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/66486529-D98E-E011-8143-0030487CD6DA.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/54D8B486-D68E-E011-A2DE-0030487CD6DA.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/4E44DA77-D78E-E011-997B-003048F118D2.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/469D9959-D68E-E011-B661-00304879EDEA.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/30478A0F-D88E-E011-A588-0030487CD6D8.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/2854E6C9-D58E-E011-A5E5-0030487CD17C.root',
       '/store/relval/CMSSW_4_2_3_SLHC2/RelValTTbar_Tauola/GEN-SIM/DESIGN42_V11_110603_special-v1/0000/12705C98-D18E-E011-9D2D-0030487C6A66.root' 
    )
)
# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    #outputCommands = cms.untracked.vstring('keep *','drop *_mix_*_*'),
    fileName = cms.untracked.string('file:reco.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    )
)
#I'm only interested in the validation stuff
#process.output.outputCommands = cms.untracked.vstring('drop *','keep *_MEtoEDMConverter_*_*')

#process.output = cms.OutputModule("PoolOutputModule",
#         outputCommands = process.AODSIMEventContent.outputCommands,
#         fileName = cms.untracked.string(
#		'file:/uscms_data/d2/brownson/slhc/quadMuon_RECO.root')
#)


# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'DESIGN42_V11::All'

### PhaseI Geometry and modifications ###############################################
#process.load("SLHCUpgradeSimulations.Geometry.Phase1_R39F16_cmsSimIdealGeometryXML_cff")
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

process.muons.TrackerKinkFinderParameters.TrackerRecHitBuilder = cms.string('WithTrackAngle')
# The SeedMergerPSet should be added to the following file for Phase 1
# RecoTracker/SpecialSeedGenerators/python/CombinatorialSeedGeneratorForCosmicsRegionalReconstruction_cfi.py
# but pixel layers are not used here for cosmic TODO: need Maria and Jan to do appropriate thing here
process.regionalCosmicTrackerSeeds.SeedMergerPSet = cms.PSet(
	mergeTriplets = cms.bool(False),
	ttrhBuilderLabel = cms.string( "PixelTTRHBuilderWithoutAngle" ),
	addRemainingTriplets = cms.bool(False),
	layerListName = cms.string( "PixelSeedMergerQuadruplets" )
	)
process.regionalCosmicTracks.TTRHBuilder = cms.string('WithTrackAngle')


## when using the SV producer fix from later CMSSW_4_2_1 tag
process.secondaryVertexTagInfos.beamSpotTag = cms.InputTag('offlineBeamSpot')
process.ghostTrackVertexTagInfos.beamSpotTag = cms.InputTag('offlineBeamSpot')

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
## already in fake conditions don't add here (might overide?)
#process.load("RecoVertex.BeamSpotProducer.BeamSpotFakeParameters_cfi")

## need this at the end as the validation config redefines random seed with just mix
process.load("IOMC.RandomEngine.IOMC_cff")

### back to standard job commands ##################################################
process.pdigi.remove(process.simCastorDigis)
process.DigiToRaw.remove(process.castorRawData)

process.DigiToRaw.remove(process.siPixelRawData)
process.RawToDigi.remove(process.siPixelDigis)

process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)

# Path and EndPath definitions
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)

process.reconstruction_step 	= cms.Path(process.reconstruction)
process.mix_step 		= cms.Path(process.mix)
#process.reconstruction_step 	= cms.Path(process.trackerlocalreco*
#						process.offlineBeamSpot+
#                                                process.recopixelvertexing*process.ckftracks_wodEdXandSteps4and5)
process.debug_step 		= cms.Path(process.anal)
#process.validation_step 	= cms.Path(process.cutsTPEffic*
#						process.cutsTPFake*
#						process.slhcTracksValidation)
process.user_step 		= cms.Path(process.ReadLocalMeasurement)
process.endjob_step 		= cms.Path(process.endOfProcess)
process.out_step 		= cms.EndPath(process.output)

# Schedule definition
#process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.mix_step,process.reconstruction_step,process.validation_step,process.user_step,process.endjob_step,process.out_step)
#process.schedule = cms.Schedule(process.mix_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.out_step)
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.out_step)

