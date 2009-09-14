import FWCore.ParameterSet.Config as cms
import RecoBTag.SoftLepton.muonSelection

# SoftLeptonTagInfo producer for tagging caloJets with global muons 
softMuonTagInfos = cms.EDFilter("SoftLepton",
    jets = cms.InputTag("ak5CaloJets"),
    leptons = cms.InputTag("muons"),
    primaryVertex = cms.InputTag("offlinePrimaryVertices"),
    
    refineJetAxis = cms.uint32(0),          # use calorimetric jet direction by default

    leptonQualityCut = cms.double(0.5),
    leptonDeltaRCut = cms.double(0.4),      # lepton distance from jet axis
    leptonChi2Cut = cms.double(9999.0),     # no cut on lepton's track's chi2/ndof
    
    muonSelection = RecoBTag.SoftLepton.muonSelection.AllGlobalMuons
)
