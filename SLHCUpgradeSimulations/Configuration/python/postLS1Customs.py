
import FWCore.ParameterSet.Config as cms

from muonCustoms import customise_csc_geom_cond_digi


def digiEventContent(process):
    #extend the event content

    alist=['RAWSIM','FEVTDEBUG','FEVTDEBUGHLT','GENRAW','RAWSIMHLT','FEVT']
    for a in alist:
        b=a+'output'
        if hasattr(process,b):
            getattr(process,b).outputCommands.append('keep *_simMuonCSCDigis_*_*')
            getattr(process,b).outputCommands.append('keep *_simMuonRPCDigis_*_*')
            getattr(process,b).outputCommands.append('keep *_simHcalUnsuppressedDigis_*_*')

    return process    

def digiCustomsRelVal(process):
    #deal with csc
    process=customise_csc_geom_cond_digi(process)
    process=digiEventContent(process)
    process.digi2raw_step.remove(process.cscpacker)
    return process

def digiCustoms(process):
    process=turnOffXFrame(process)
    process=digiCustomsRelVal(process)
    return process


def hltCustoms(process):
    process.CSCGeometryESModule.useGangedStripsInME1a = False

    process.hltCsc2DRecHits.readBadChannels = cms.bool(False)
    process.hltCsc2DRecHits.CSCUseGasGainCorrection = cms.bool(False)

    # Switch input for CSCRecHitD to  s i m u l a t e d  digis

    process.hltCsc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    process.hltCsc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")

    return process

def recoCustoms(process):

    # ME1/1A is  u n g a n g e d  Post-LS1

    process.CSCGeometryESModule.useGangedStripsInME1a = False

    # Turn off some flags for CSCRecHitD that are turned ON in default config

    process.csc2DRecHits.readBadChannels = cms.bool(False)
    process.csc2DRecHits.CSCUseGasGainCorrection = cms.bool(False)

    # Switch input for CSCRecHitD to  s i m u l a t e d  digis

    process.csc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    process.csc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")

    return process

def recoOutputCustoms(process):

    alist=['AODSIM','RECOSIM','FEVTSIM','FEVTDEBUG','FEVTDEBUGHLT','RECODEBUG','RAWRECOSIMHLT','RAWRECODEBUGHLT']
    for a in alist:
        b=a+'output'
        if hasattr(process,b):
            getattr(process,b).outputCommands.append('keep *_simMuonCSCDigis_*_*')
            getattr(process,b).outputCommands.append('keep *_simMuonRPCDigis_*_*')
            getattr(process,b).outputCommands.append('keep *_simHcalUnsuppressedDigis_*_*')
            getattr(process,b).outputCommands.append('keep *_rawDataCollector_*_*')
    return process

            
