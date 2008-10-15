import FWCore.ParameterSet.Config as cms

pfRecoTauDiscriminationAgainstMuon = cms.EDProducer("PFRecoTauDiscriminationAgainstMuon",
    a = cms.double(0.5), ## user definde 2D Cut. Reject tau if calo * a + seg * b < 0 

    b = cms.double(0.5),
    c = cms.double(0.0),
    discriminatorOption = cms.string('noSegMatch'), ## available options are; noSegMatch, twoDCut, merePresence, combined

    PFTauProducer = cms.string('pfRecoTauProducer')
)


