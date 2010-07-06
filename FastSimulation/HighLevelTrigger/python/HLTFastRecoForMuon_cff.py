import FWCore.ParameterSet.Config as cms

# RecoMuon flux ##########################################################
# L2 seeds from L1 input
# module hltL2MuonSeeds = L2MuonSeeds from "RecoMuon/L2MuonSeedGenerator/data/L2MuonSeeds.cfi"
# replace hltL2MuonSeeds.GMTReadoutCollection = l1extraParticles
# replace hltL2MuonSeeds.InputObjects = l1extraParticles
# L3 regional reconstruction
from FastSimulation.Muons.L3Muons_cff import *
import FastSimulation.Muons.L3Muons_cfi
hltL3Muons = FastSimulation.Muons.L3Muons_cfi.L3Muons.clone()
hltL3Muons.L3TrajBuilderParameters.TrackTransformer.TrackerRecHitBuilder = 'WithoutRefit'
hltL3Muons.L3TrajBuilderParameters.TrackerRecHitBuilder = 'WithoutRefit'
hltL3Muons.TrackLoaderParameters.beamSpot = cms.InputTag("offlineBeamSpot")

# L3 regional seeding, candidating, tracking
#--the two below have to be picked up from confDB: 
# from FastSimulation.Muons.TSGFromL2_cfi import *
# from FastSimulation.Muons.HLTL3TkMuons_cfi import *
from FastSimulation.Muons.TrackCandidateFromL2_cfi import *
from FastSimulation.Muons.TSGFromL2_cfi import *
hltL3TrajectorySeed = FastSimulation.Muons.TSGFromL2_cfi.hltL3TrajectorySeed.clone()

