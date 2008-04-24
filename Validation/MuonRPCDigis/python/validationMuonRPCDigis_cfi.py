import FWCore.ParameterSet.Config as cms

validationMuonRPCDigis = cms.EDFilter("RPCDigiValid",
    # Label to retrieve Digis from the event 
    digiLabel = cms.untracked.string('simMuonRPCDigis'),
    # Name of the root file which will contain the histos
    outputFile = cms.untracked.string('rpcDigiValidPlots.root')
)



