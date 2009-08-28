import FWCore.ParameterSet.Config as cms

muonHistManager = cms.PSet(
    pluginName = cms.string('muonHistManager'),
    pluginType = cms.string('MuonHistManager'),
      
    muonSource = cms.InputTag('cleanLayer1MuonsSel'),
    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    genParticleSource = cms.InputTag('genParticles'),

    dqmDirectory_store = cms.string('MuonQuantities'),

    #requireGenMuonMatch = cms.bool(True),
    requireGenMuonMatch = cms.bool(False),

    #normalization = cms.string("muons"),
    normalization = cms.string("events"),

    #makeIsoPtCtrlHistograms = cms.bool(False),
    makeIsoPtCtrlHistograms = cms.bool(True),
    makeIsoPtConeSizeDepHistograms = cms.bool(False)
)
