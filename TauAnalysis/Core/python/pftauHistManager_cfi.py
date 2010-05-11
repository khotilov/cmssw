import FWCore.ParameterSet.Config as cms

tauHistManager = cms.PSet(    
    pluginName = cms.string('tauHistManager'),
    pluginType = cms.string('TauHistManager'),
      
    tauSource = cms.InputTag('cleanPatTaus'),
    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    jetSource = cms.InputTag('selectedLayer1JetsEt20Cumulative'),
    genParticleSource = cms.InputTag('genParticles'),
  
    tauIndicesToPlot = cms.string('all'),

    dqmDirectory_store = cms.string('TauQuantities'),

    #requireGenTauMatch = cms.bool(True),
    requireGenTauMatch = cms.bool(False),

    skipPdgIdsGenParticleMatch = cms.vint32( 12, 14, 16 ),
    
    #normalization = cms.string("taus"),
    normalization = cms.string("events"),

    makeIsoPtCtrlHistograms = cms.bool(False),
    makeIsoPtConeSizeDepHistograms = cms.bool(False),

    makeTauIdEfficiencyHistograms = cms.vstring(
        #"effByStandardChainZTTsim"
    ),

    makeTauFakeRateHistograms = cms.vstring(
        #"frByStandardChainMuEnrichedQCDsim",
        #"frByStandardChainDiJetHighPtsim",
        #"frByStandardChainDiJetSecondPtsim",    
        #"frByStandardChainWJetssim"
    )
)
