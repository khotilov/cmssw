import FWCore.ParameterSet.Config as cms

electronHistManager = cms.PSet(
  name = cms.string('electronHistManager'),
  type = cms.string('ElectronHistManager'),
      
  electronSource = cms.InputTag('allLayer1ElectronsSel'),
  vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),

  dqmDirectory_store = cms.string('ElectronQuantities'),

  #requireGenElectronMatch = cms.bool(True)
  requireGenElectronMatch = cms.bool(False)
)
