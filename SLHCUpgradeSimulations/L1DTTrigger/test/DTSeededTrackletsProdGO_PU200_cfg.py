import FWCore.ParameterSet.Config as cms

################################################################################
maxEvNr = 100

PT = '10'
PU = 'PU200'

#Stubs = "fromSimHits"
Stubs = "fromPixelDigis"
#Stubs = "fromTrackTriggerHits"

USE_TSTHETA = True    ### True/False!
USE_roughTHETA = False


################################################################################
 
use_tstheta_label = ""
if USE_TSTHETA == True:
        use_tstheta_label = "_tsTheta"
if USE_roughTHETA == True:
        use_theta_label = "_roughTheta"

inputDir = "/lustre/cmswork/zotto/SLHCdata/"
outputDir = "/lustre/cmswork/zotto/SLHCdata/"
#inputDir = "/lustre/cmsdata/zotto/SLHCdata/"
#outputDir = "/lustre/cmsdata/zotto/SLHCdata/"
#outpuDir ="."

#singleMuonPt = PT


rootFileNamesIn = cms.untracked.vstring(
      'file:' 
      +   inputDir + 
     'FullSim_Pt'+ PT + '_' + PU + '.root')
     
rootFileNameOut = cms.untracked.string(
    'file:' 
 + outputDir  +
     'out_mu_'+ PT + 'GeV_' + PU + '.root')
#     'out_mu_PT_GeV_' + PU + '.root') 
	 
rootFileNameHis = cms.untracked.string(
     'h_out_mu_'+ PT + 'GeV_' + PU + '.root')
#     'h_out_mu_PT_GeV_' + PU + '.root') 

asciiFileName = cms.untracked.string(
#     'file:' + outputDir  +
     'out_mu_'+ PT + 'GeV_' + PU + '.txt')
#			'out_mu_PT_GeV_' + PU + '.txt')
#
msgFileNames = cms.untracked.vstring(
#     'file:' + outputDir  +
     'out_mu_'+ PT + 'GeV_' + PU + '.msg')
#			'out_mu_PT_GeV_' + PU + '.msg')

    
# rootFileNameOut = cms.untracked.string(
#     outputDir + 'PT'+ 'GeVPtMuons_DTSeededTracklets_fromPixelDigis' 
#     + use_tstheta_label + '.root'
# ) 
	 
rootFileNamePlots = cms.untracked.string(
 #    outputDir + 
     'Plots_' + PT + 'GeVMuons' + use_tstheta_label + '.root'
) 

patternFileName = cms.untracked.string(
#     outputDir + 
     PT + 'GeV_Pt_Muons_ascii' + '.pat')

# asciiFileName = cms.untracked.string(
#      outputDir + 'PT' + 'GeVPtMuons_DTSeededTracklets_fromPixelDigis'
#      + use_tstheta_label + '.outAscii'
# )

# msgFileNames = cms.untracked.vstring(
#      outputDir + 'PT' + 'GeVPtMuons_DTSeededTracklets_fromPixelDigis'
#      + use_tstheta_label + '.msg'
# )

################################################################################

process = cms.Process("MuonsUpgrade")

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# Other statements --------------------------------------------------------------
#process.GlobalTag.globaltag = 'MC_42_V12::All'
#process.GlobalTag.globaltag = "DESIGN42_V11::All"

process.load("SLHCUpgradeSimulations.L1DTTrigger.DTstackedTracker_specific_cff")
process.load("SLHCUpgradeSimulations.L1DTTrigger.DTSeededTrackletsProd_cff")
# process.load("SLHCUpgradeSimulations.L1DTTrigger.DTSeededTrackletsPlots_cfi")
process.load("SLHCUpgradeSimulations.L1DTTrigger.DTSeededTrackletsOutMod_cff")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvNr))
process.MessageLogger.destinations = msgFileNames
process.source.fileNames = rootFileNamesIn

# to produce, in case, collection of L1MuDTTrack objects:    
#process.dttfDigis = cms.Path(process.simDttfDigis)

# muon specific ----------------------------------------------------------------
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("SimMuon/DTDigitizer/muonDTDigis_cfi")
#process.load("L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfig_cff")
## # tracker specific -------------------------------------------------------------
process.load('SLHCUpgradeSimulations.Geometry.Longbarrel_cmsSimIdealGeometryXML_cff')
process.load("SLHCUpgradeSimulations.Utilities.StackedTrackerGeometry_cfi")
process.load('SLHCUpgradeSimulations.Geometry.fakeConditions_Longbarrel_cff')
process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_Longbarrel_cff")
process.load("SLHCUpgradeSimulations.Geometry.upgradeTracking_longbarrel_cff")
process.load("SimTracker.Configuration.SimTracker_cff") # era commentato!!!!!!!!!!!!!!!!

if USE_TSTHETA == True:
	process.DTL1slhcProd.use_TSTheta = cms.untracked.bool(True)
if USE_roughTHETA == True:
        process.DTL1slhcProd.use_roughTheta = cms.untracked.bool(True)
process.DTL1slhcProd.rootFileNameOut = rootFileNameOut
process.DTL1slhcProd.asciiFileName   = asciiFileName
process.DTL1slhcProd.rootFileNameHis = rootFileNameHis
process.DTL1slhcProd.patternAsciiFileName = patternFileName 
#-----------------------------------------------
process.DTL1slhcProd.singleMuonPt = cms.untracked.double(140)
#-----------------------------------------------
process.DTL1slhcProd.min_invRb = cms.untracked.double(0.00000045)
process.DTL1slhcProd.max_invRb = cms.untracked.double(1.)
# Assumed radius of boundary surface of the magnetic field:
process.DTL1slhcProd.Erre      = cms.untracked.double(360.0)
# Correction factor computing Pt using third "DT muon" method:
process.DTL1slhcProd.station2_correction = cms.untracked.double(1.23)
process.DTL1slhcProd.third_method_accurate = cms.untracked.bool(False) 
# To apply PT cut to stubs
process.DTL1slhcProd.magneticFieldStrength = cms.untracked.double(4.0112)
process.DTL1slhcProd.ptThreshold           = cms.untracked.double(10.0)

process.scms = cms.Path(process.DTL1slhcProd)


#process.DTL1slhcPlots.rootFileNamePlots = rootFileNamePlots
#process.DTL1slhcPlots.desert = cms.untracked.double(0.005)
#process.DTL1slhcPlots.singleMuonPt = singleMuonPt
#process.DTL1slhcPlots.singleMuonPt = cms.untracked.double(140)
#process.DTL1slhcPlots.optimize_plot = True
#process.plots = cms.Path(process.DTL1slhcPlots)

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( 
    Rethrow = cms.untracked.vstring('ProductNotFound') 
)

process.Out.fileName = rootFileNameOut
process.end = cms.EndPath(process.Out)

process.schedule = cms.Schedule(process.scms, process.end)
#process.schedule = cms.Schedule(process.scms, process.plots, process.end)
#process.schedule = cms.Schedule(process.dttfDigis, process.scms, process.end)


