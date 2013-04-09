import FWCore.ParameterSet.Config as cms

muonsCleaned = cms.EDProducer('PFMuonUpdater',
                       pfMuons = cms.InputTag('pfMuons'),
                       muons = cms.InputTag('muons')
)
