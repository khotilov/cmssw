import FWCore.ParameterSet.Config as cms

from FWCore.MessageLogger.MessageLogger_cfi import *
process = cms.Process("IGUANA")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("VisFramework.VisFrameworkBase.VisConfigurationService_cff")
process.source = cms.Source("PoolSource",
                            # For NAF tutorial
                            fileNames = cms.untracked.vstring('/store/data/Commissioning08/Cosmics/RECO/CRAFT_V3P_TrackerPointing_v1/0010/E8D35F0C-03A5-DD11-85BD-003048D15E52.root')                  
                            # if you want to try it on lxplus use this instead
                            #fileNames = cms.untracked.vstring('/store/data/Commissioning08/Cosmics/RAW-RECO/CRAFT_ALL_V12_229_Tosca090322_ReReco_FromTrackerPointing_v1/0005/EEEB64C0-3E37-DE11-84A1-001A92971AA8.root')
)
process.VisConfigurationService.EnabledTwigs = ('* Muon *',
                                                '* Track *',
                                                '*Segment*',
                                                '* Muon *',
                                                '* Tracking RecHit *',
                                                '/Objects/CMS Event and Detector/Magnet',
                                                '/Objects/Event Collections/Run and Event Number')

