####################################################
# SLHCUpgradeSimulations                           #
# Configuration file for Full Workflow             #
# Step 2 (again)                                   #
# Understand if everything is fine with            #
# L1TkCluster e L1TkStub                           #
####################################################
# Nicola Pozzobon                                  #
# CERN, August 2012                                #
####################################################

#################################################################################################
# import of general framework
#################################################################################################
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process('ValidatePixelDigiMaps')

#################################################################################################
# global tag
#################################################################################################
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')

#################################################################################################
# load the specific tracker geometry
#################################################################################################
process.load('Configuration.Geometry.GeometryExtendedPhase2TkLB4LPS_2L2SReco_cff')
process.load('Configuration.Geometry.GeometryExtendedPhase2TkLB4LPS_2L2S_cff')
#process.load("SLHCUpgradeSimulations.Utilities.StackedTrackerGeometry_cfi")

#process.TrackerDigiGeometryESModule.applyAlignment = False
#process.TrackerDigiGeometryESModule.fromDDD = cms.bool(True)
#process.TrackerGeometricDetESModule.fromDDD = cms.bool(True)

#################################################################################################
# load the magnetic field
#################################################################################################
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('Configuration.StandardSequences.MagneticField_40T_cff')
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True

#################################################################################################
# define the source and maximum number of events to generate and simulate
#################################################################################################
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
#    fileNames = cms.untracked.vstring('file:TenMuPt_0_50_ExtendedPhase2BE_500_L1TT.root')
#    fileNames = cms.untracked.vstring('file:TenMuPt_0_50_ExtendedPhase2BE_500_RAW2DIGI_L1Reco_RECO.root')
     fileNames = cms.untracked.vstring('file:TenMuPt_0_50_ExtendedPhase2TkLB4PS22S_5000_DIGI_L1_DIGI2RAW_L1TT_RECO.root')

)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

#################################################################################################
# load the analyzer
#################################################################################################
process.ValidatePixelDigiMaps = cms.EDAnalyzer("ValidatePixelDigiMaps",
#    DebugMode = cms.bool(True)
)

#################################################################################################
# define output file and message logger
#################################################################################################
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('file:ValidatePixelDigiMaps_ExtendedPhase2TkLB4PS22S.root')
)

#################################################################################################
# define the final path to be fed to cmsRun
#################################################################################################
process.p = cms.Path( process.ValidatePixelDigiMaps )

