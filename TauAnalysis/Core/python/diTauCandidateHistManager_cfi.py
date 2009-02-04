import FWCore.ParameterSet.Config as cms

diTauCandidateHistManager = cms.PSet(    
  name = cms.string('diTauCandidateHistManager'),
  type = cms.string('DiCandidatePairHistManager'),
      
  diTauCandidateSource = cms.InputTag(''),

  dqmDirectory_store = cms.string('DiTauCandidateQuantities'),

  #requireGenMatch = cms.bool(True)
  requireGenMatch = cms.bool(False)
)
