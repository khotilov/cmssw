import FWCore.ParameterSet.Config as cms

#Full Event content 
RecoTauTagFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_iterativeCone5PFJets_*_*', 
        'keep *_ic5PFJetTracksAssociatorAtVertex_*_*', 
        'keep *_pfRecoTauTagInfoProducer_*_*', 
        'keep *_pfRecoTauProducer*_*_*', 
        'keep *_pfRecoTauDiscrimination*_*_*', 
        'keep *_caloRecoTauTagInfoProducer_*_*', 
        'keep *_caloRecoTauProducer*_*_*', 
        'keep *_caloRecoTauDiscrimination*_*_*')
)
#RECO content
RecoTauTagRECO = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_iterativeCone5PFJets_*_*', 
        'keep *_ic5PFJetTracksAssociatorAtVertex_*_*', 
        'keep *_pfRecoTauTagInfoProducer_*_*', 
        'keep *_pfRecoTauProducer*_*_*', 
        'keep *_pfRecoTauDiscrimination*_*_*', 
        'keep *_caloRecoTauTagInfoProducer_*_*', 
        'keep *_caloRecoTauProducer*_*_*', 
        'keep *_caloRecoTauDiscrimination*_*_*')
)
#AOD content
RecoTauTagAOD = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep *_iterativeCone5PFJets_*_*', 
        'keep *_ic5PFJetTracksAssociatorAtVertex_*_*', 
        'keep *_pfRecoTauTagInfoProducer_*_*', 
        'keep *_pfRecoTauProducer*_*_*', 
        'keep *_pfRecoTauDiscrimination*_*_*', 
        'keep *_caloRecoTauTagInfoProducer_*_*', 
        'keep *_caloRecoTauProducer*_*_*', 
        'keep *_caloRecoTauDiscrimination*_*_*')
)

