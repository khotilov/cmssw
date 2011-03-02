import FWCore.ParameterSet.Config as cms

rpcChamberQuality = cms.EDAnalyzer("RPCChamberQuality",
                                   OfflineDQM = cms.untracked.bool(True),
                                   PrescaleFactor  = cms.untracked.int32(5),
                                   MinimumRPCEvents = cms.untracked.int32(10000),
                                   RecHitTypeFolder = cms.untracked.string("Noise")
                                   )


rpcMuonChamberQuality = cms.EDAnalyzer("RPCChamberQuality",
                                       OfflineDQM = cms.untracked.bool(True),
                                       PrescaleFactor  = cms.untracked.int32(5),
                                       MinimumRPCEvents = cms.untracked.int32(10000),
                                       RecHitTypeFolder = cms.untracked.string("Muon")
                                       )
