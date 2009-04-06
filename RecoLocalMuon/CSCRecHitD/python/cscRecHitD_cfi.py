import FWCore.ParameterSet.Config as cms

# parameters for CSC rechit building
from RecoLocalMuon.CSCRecHitD.cscRecHitD_cff import *
csc2DRecHits = cms.EDProducer("CSCRecHitDProducer",
    #
    #    Parameters for coordinate and uncertainty calculations
    #    Data and MC parameters are (still) different
    #    Needs tuning
    #
    cscRecHitDParameters,
    #
    #    Parameters for strip hits
    #
    CSCStripPeakThreshold = cms.double(10.0),
    CSCStripClusterChargeCut = cms.double(25.0),
    CSCStripxtalksOffset = cms.double(0.03),
    #
    #    How to find SCA peak time?
    #                              
    UseAverageTime = cms.bool(False),
    UseParabolaFit = cms.bool(False),
    UseFourPoleFit = cms.bool(True),                       
    #
    #    Parameters for wire hits
    CSCWireClusterDeltaT = cms.int32(1),
    #
    #    Calibration info:
    CSCUseCalibrations = cms.bool(True),
    #    Pedestal treatment
    CSCUseStaticPedestals = cms.bool(False),
    CSCNoOfTimeBinsForDynamicPedestal = cms.int32(2),
    #
    #    Which digis:
    #
    #  When using data from unpacker
    wireDigiTag = cms.InputTag("muonCSCDigis","MuonCSCWireDigi"),
    stripDigiTag = cms.InputTag("muonCSCDigis","MuonCSCStripDigi"),
    #  When using data from simulation
    #    wireDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi"),
    #    stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi"),
                              
    #
    #    Parameters which are not used currently
    #
    CSCDebug = cms.untracked.bool(False),
    readBadChannels = cms.bool(False),
    readBadChambers = cms.bool(False),
    #  To be set once wire digis have proper timing info:
    CSCstripWireDeltaTime = cms.int32(8),
    # to be deleted
    CSCStripClusterSize = cms.untracked.int32(3)
)


