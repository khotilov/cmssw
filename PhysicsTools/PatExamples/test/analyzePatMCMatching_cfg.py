import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

#------------------------------------------------------------------------------------------------------
# To configure the Matching, we have to configure the PAT-Workflow starting from the patDefaultSequence:
#------------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.patTemplate_cfg import *

## prepare several clones of match associations for status 1, 3 and in flight muons (status -1)
## add inFlightMuons
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.inFlightMuons = cms.EDProducer("PATGenCandsFromSimTracksProducer",
        src           = cms.InputTag("g4SimHits"),   ## use "famosSimHits" for FAMOS
        setStatus     = cms.int32(-1),
        particleTypes = cms.vstring("mu+"),          ## picks also mu-, of course
        filter        = cms.vstring("pt > 0.5"),     ## just for testing
        makeMotherLink = cms.bool(True),
        writeAncestors = cms.bool(True),             ## save also the intermediate GEANT ancestors of the muons
        genParticles   = cms.InputTag("genParticles"),
)

process.muMatch3 = process.muonMatch.clone(mcStatus = [3]) # hard scattering
process.muMatch1 = process.muonMatch.clone(mcStatus = [1]) # stable
process.muMatchF = process.muonMatch.clone(mcStatus = [-1], matched = cms.InputTag("inFlightMuons"))

process.muMatch1.checkCharge = False
process.muMatch3.checkCharge = False

#process.muMatch3.resolveByMatchQuality = True
#process.muMatch1.resolveByMatchQuality = True

process.muMatch3.maxDeltaR = 0.05
process.muMatch3.maxDPtRel = 0.1

process.muMatch1.maxDeltaR = 0.05
process.muMatch1.maxDPtRel = 0.1

process.muMatchF.maxDeltaR = 0.3
process.muMatchF.maxDPtRel = 0.2


process.muonMatchByPt = cms.EDProducer("MCMatcherByPt", # cut on deltaR, deltaPt/Pt; pick best by deltaPt
    src     = cms.InputTag("muons"), # RECO objects to match  
    matched = cms.InputTag("genParticles"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
    maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
)


## add the new matches to the default sequence
process.patDefaultSequence.replace(process.muonMatch,
                                   process.muMatch1 +
                                   process.muMatch3 +
                                   process.muMatchF
                                  #+ process.muonMatchByPt
)

process.patMuons.genParticleMatch = cms.VInputTag(
    cms.InputTag("muMatch3"),
    cms.InputTag("muMatch1"),
    cms.InputTag("muMatchF")
    #, cms.InputTag("muonMatchByPt")
)


#-----------------------------------------
# As usual add those two usefull things:
#----------------------------------------

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('analyzePatMCMatching.root')
)

process.MessageLogger = cms.Service("MessageLogger")


#----------------------------------------------------------------------
# Finally let's analyze the matching and run all that in correct order:
#----------------------------------------------------------------------

process.analyzePatMCMatching = cms.EDAnalyzer("PatMCMatching",
  muonSrc     = cms.untracked.InputTag("cleanPatMuons")                                             
)

process.out.outputCommands = cms.untracked.vstring('keep *') 
process.outpath.remove(process.out)


process.p = cms.Path(process.inFlightMuons + process.patDefaultSequence + process.analyzePatMCMatching)


#----------------------------------------------------------------------
# Change the input file to compare
#----------------------------------------------------------------------
process.maxEvents.input = -1

process.source.fileNames =[
# 'rfio:///castor/cern.ch/cms/store/relval/CMSSW_3_6_2/RelValJpsiMM/GEN-SIM-RECO/START36_V10-v1/0002/04DF42BB-F970-DF11-9B1D-00304867C1B0.root'

'rfio:///castor/cern.ch/cms/store/relval/CMSSW_3_6_2/RelValZMM/GEN-SIM-RECO/START36_V10-v1/0002/16F4C9D1-1B71-DF11-B488-0018F3D095FE.root'
 ]
del(process.out)
del(process.outpath)
