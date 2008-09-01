# AlCaReco for muon based alignment using beam-halo muons

import FWCore.ParameterSet.Config as cms

# import HLTrigger.HLTfilters.hltHighLevel_cfi
# ALCARECOMuAlBeamHaloHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# ALCARECOMuAlBeamHaloHLT.HLTPaths = ['HLT_CSCBeamHalo', 'HLT_CSCBeamHaloRing2or3']

from Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi import *
from RecoMuon.DetLayers.muonDetLayerGeometry_cfi import *

import RecoLocalMuon.CSCSegment.cscSegments_cfi
cscSegmentsALCARECOBH = RecoLocalMuon.CSCSegment.cscSegments_cfi.cscSegments.clone()
cscSegmentsALCARECOBH.algo_type = 4 # Choice of the building algo: 1 SK, 2 TC, 3 DF, 4 ST, ...
cscSegmentsALCARECOBH.inputObjects = 'csc2DRecHits'

import RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi
CosmicMuonSeedALCARECOBH = RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi.CosmicMuonSeed.clone()
CosmicMuonSeedALCARECOBH.EnableDTMeasurement = False
CosmicMuonSeedALCARECOBH.CSCRecSegmentLabel = 'cscSegmentsALCARECOBH'

import RecoMuon.CosmicMuonProducer.cosmicMuons_cfi
cosmicMuonsALCARECOBH = RecoMuon.CosmicMuonProducer.cosmicMuons_cfi.cosmicMuons.clone()
cosmicMuonsALCARECOBH.MuonSeedCollectionLabel = 'CosmicMuonSeedALCARECOBH'
cosmicMuonsALCARECOBH.TrajectoryBuilderParameters.EnableDTMeasurement = False
cosmicMuonsALCARECOBH.TrajectoryBuilderParameters.CSCRecSegmentLabel = 'cscSegmentsALCARECOBH'

reconstructAsCosmicMuonsALCARECOBH = cms.Sequence(cscSegmentsALCARECOBH * CosmicMuonSeedALCARECOBH * cosmicMuonsALCARECOBH)

ALCARECOMuAlBeamHalo = cms.EDFilter("AlignmentCSCBeamHaloSelectorModule",
    filter = cms.bool(True),
    src = cms.InputTag("cosmicMuonsALCARECOBH"),
    minStations = cms.uint32(0), # no "energy cut" yet
    minHitsPerStation = cms.uint32(4)
)

# seqALCARECOMuAlBeamHalo = cms.Sequence(ALCARECOMuAlBeamHaloHLT + reconstructAsCosmicMuonsALCARECOBH * ALCARECOMuAlBeamHalo)
seqALCARECOMuAlBeamHalo = cms.Sequence(reconstructAsCosmicMuonsALCARECOBH * ALCARECOMuAlBeamHalo)

