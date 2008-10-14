import FWCore.ParameterSet.Config as cms

rpcRecHits = cms.EDProducer("RPCRecHitProducer",
    recAlgoConfig = cms.PSet(

    ),
    recAlgo = cms.string('RPCRecHitStandardAlgo'),
    rpcDigiLabel = cms.InputTag("muonRPCDigis"),
    maskmapfile = cms.FileInPath('RecoLocalMuon/RPCRecHit/data/RPCMaskedStrips.dat')
)


