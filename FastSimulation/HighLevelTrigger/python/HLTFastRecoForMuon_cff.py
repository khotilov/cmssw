import FWCore.ParameterSet.Config as cms

# RecoMuon flux ##########################################################
# L2 seeds from L1 input
# module hltL2MuonSeeds = L2MuonSeeds from "RecoMuon/L2MuonSeedGenerator/data/L2MuonSeeds.cfi"
# replace hltL2MuonSeeds.GMTReadoutCollection = l1extraParticles
# replace hltL2MuonSeeds.InputObjects = l1extraParticles
# L3 regional reconstruction
from FastSimulation.Muons.L3Muons_cff import *
#import FastSimulation.Muons.L3Muons_cfi
#hltL3Muonstemp = FastSimulation.Muons.L3Muons_cfi.L3Muons.clone()
#hltL3Muonstemp.L3TrajBuilderParameters.TrackTransformer.TrackerRecHitBuilder = 'WithoutRefit'
#hltL3Muonstemp.L3TrajBuilderParameters.TrackerRecHitBuilder = 'WithoutRefit'
#hltL3Muonstemp.TrackLoaderParameters.beamSpot = cms.InputTag("offlineBeamSpot")

# L3 regional seeding, candidating, tracking
#--the two below have to be picked up from confDB: 
# from FastSimulation.Muons.TSGFromL2_cfi import *
# from FastSimulation.Muons.HLTL3TkMuons_cfi import *
#from FastSimulation.Muons.TrackCandidateFromL2_cfi import *
#from FastSimulation.Muons.TSGFromL2_cfi import *
#hltL3TrajectorySeed = FastSimulation.Muons.TSGFromL2_cfi.hltL3TrajectorySeed.clone()

import FastSimulation.Muons.TSGFromL2_cfi as TSG
#from FastSimulation.Muons.TSGFromL2_cfi import OIStatePropagators as OIProp
from FastSimulation.Muons.TSGFromL2_cfi import OIHitPropagators as OIHProp
## Make three individual seeders
## OIState can be taken directly from configuration
#hltL3TrajSeedOIState = TSG.l3seeds("OIState")
#hltL3TrajSeedOIState.ServiceParameters.Propagators = cms.untracked.vstring()
#OIProp(hltL3TrajSeedOIState,hltL3TrajSeedOIState.TkSeedGenerator)
hltL3TrajSeedOIHit = TSG.l3seeds("OIHitCascade")
hltL3TrajSeedOIHit.ServiceParameters.Propagators = cms.untracked.vstring()
OIHProp(hltL3TrajSeedOIHit,hltL3TrajSeedOIHit.TkSeedGenerator.iterativeTSG)
hltL3TrajSeedIOHit = TSG.l3seeds("IOHitCascade")

## Make one TrackCand for each seeder
from FastSimulation.Muons.TrackCandidateFromL2_cfi import *
hltL3TrackCandidateFromL2OIState = FastSimulation.Muons.TrackCandidateFromL2_cfi.hltL3TrackCandidateFromL2.clone()
hltL3TrackCandidateFromL2OIState.SeedProducer = "hltL3TrajSeedOIState"
hltL3TrackCandidateFromL2OIHit = FastSimulation.Muons.TrackCandidateFromL2_cfi.hltL3TrackCandidateFromL2.clone()
hltL3TrackCandidateFromL2OIHit.SeedProducer = "hltL3TrajSeedOIHit"    
hltL3TrackCandidateFromL2IOHit = FastSimulation.Muons.TrackCandidateFromL2_cfi.hltL3TrackCandidateFromL2.clone()
hltL3TrackCandidateFromL2IOHit.SeedProducer = "hltL3TrajSeedIOHit"
hltL3TrackCandidateFromL2NoVtx = FastSimulation.Muons.TrackCandidateFromL2_cfi.hltL3TrackCandidateFromL2.clone()
hltL3TrackCandidateFromL2NoVtx.SeedProducer = "hltL3TrajectorySeedNoVtx"
