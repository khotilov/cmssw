import FWCore.ParameterSet.Config as cms

RecoVertexFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices_*_*', 
        'keep  *_offlinePrimaryVerticesWithBS_*_*',
        'keep  *_offlinePrimaryVerticesFromCosmicTracks_*_*',
        'keep  *_nuclearInteractionMaker_*_*',
        'keep *_generalV0Candidates_*_*')
)
#RECO content
RecoVertexRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices_*_*', 
        'keep  *_offlinePrimaryVerticesWithBS_*_*',
        'keep  *_offlinePrimaryVerticesFromCosmicTracks_*_*',
        'keep  *_nuclearInteractionMaker_*_*',
        'keep *_generalV0Candidates_*_*')
)
#AOD content
RecoVertexAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices_*_*', 
        'keep  *_offlinePrimaryVerticesWithBS_*_*',
        'keep  *_offlinePrimaryVerticesFromCosmicTracks_*_*',
        'keep  *_nuclearInteractionMaker_*_*',
        'keep *_generalV0Candidates_*_*')                                           
)

