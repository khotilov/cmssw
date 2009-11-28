import FWCore.ParameterSet.Config as cms

caloMEtHistManager = cms.PSet(
    pluginName = cms.string('caloMEtHistManager'),
    pluginType = cms.string('CaloMEtHistManager'),
      
    metSource = cms.InputTag('layer1METs'),
    #metSignificanceSource = cms.InputTag('met'),
    metSignificanceSource = cms.InputTag('metsignificance'),
    
    dqmDirectory_store = cms.string('CaloMEtQuantities')
)
