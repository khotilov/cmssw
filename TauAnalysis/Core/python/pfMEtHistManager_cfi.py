import FWCore.ParameterSet.Config as cms

pfMEtHistManager = cms.PSet(
    pluginName = cms.string('pfMEtHistManager'),
    pluginType = cms.string('PFMEtHistManager'),
      
    metSource = cms.InputTag('patPFMETs'),
    metSignificanceSource = cms.InputTag(''),
    
    dqmDirectory_store = cms.string('PFMEtQuantities')
)
