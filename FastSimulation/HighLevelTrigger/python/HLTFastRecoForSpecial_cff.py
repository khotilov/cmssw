import FWCore.ParameterSet.Config as cms

import FastSimulation.Tracking.TrajectorySeedProducer_cfi
# Pixel track find for minbias
#   very low pt threshold. Even if only high-pt tracks
#   are selected, the low-Pt might be wanted to check 
#   isolation of the high-Pt track.
pixelTripletSeedsForMinBias = FastSimulation.Tracking.TrajectorySeedProducer_cfi.trajectorySeedProducer.clone()
import FastSimulation.Tracking.HLTPixelTracksProducer_cfi
hltPixelTracksForMinBias = FastSimulation.Tracking.HLTPixelTracksProducer_cfi.hltPixelTracks.clone()
from FastSimulation.Tracking.PixelTracksProducer_cfi import *
import HLTrigger.HLTfilters.hltBool_cfi
hltFilterTriggerType = HLTrigger.HLTfilters.hltBool_cfi.hltBool.clone()
#--- Changes needed for HLT.cff with 1E32 menu ---#
# replace hltL1sHcalPhiSym.L1SeedsLogicalExpression = "L1_ZeroBias"
# replace hltL1sEcalPhiSym.L1SeedsLogicalExpression = "L1_ZeroBias"
# Raw data don;t exist in fast simulation -> dummy sequence for now
from FastSimulation.HighLevelTrigger.DummyModule_cfi import *
import FastSimulation.HighLevelTrigger.DummyModule_cfi
#sequence HLTIsoTrRegFEDSelection = { dummyModule }
hltSiStripRegFED = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
import FastSimulation.HighLevelTrigger.DummyModule_cfi
hltEcalRegFED = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
import FastSimulation.HighLevelTrigger.DummyModule_cfi
hltSubdetFED = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
import FastSimulation.HighLevelTrigger.DummyModule_cfi
hcalFED = FastSimulation.HighLevelTrigger.DummyModule_cfi.dummyModule.clone()
HLTPixelTrackingForMinBiasSequence = cms.Sequence(pixelTripletSeedsForMinBias*hltPixelTracksForMinBias)
pixelTripletSeedsForMinBias.numberOfHits = [3]
pixelTripletSeedsForMinBias.firstHitSubDetectorNumber = [2]
pixelTripletSeedsForMinBias.firstHitSubDetectors = [1, 2]
pixelTripletSeedsForMinBias.secondHitSubDetectorNumber = [2]
pixelTripletSeedsForMinBias.secondHitSubDetectors = [1, 2]
pixelTripletSeedsForMinBias.thirdHitSubDetectorNumber = [2]
pixelTripletSeedsForMinBias.thirdHitSubDetectors = [1, 2]
pixelTripletSeedsForMinBias.seedingAlgo = ['PixelTriplet']
pixelTripletSeedsForMinBias.originpTMin = [0.2]
pixelTripletSeedsForMinBias.pTMin = [0.2]
hltPixelTracksForMinBias.SeedProducer = cms.InputTag("pixelTripletSeedsForMinBias","PixelTriplet")


