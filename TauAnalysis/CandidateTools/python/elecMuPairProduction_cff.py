import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# produce combinations of electron + muons pairs
#--------------------------------------------------------------------------------

allElecMuPairs = cms.EDProducer("PATElecMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectronsForElecMuTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    dRmin12 = cms.double(-1.),
    srcMET = cms.InputTag('patMETs'),
    srcGenParticles = cms.InputTag('genParticles'),                                  
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                            
    verbosity = cms.untracked.int32(0)
)

produceElecMuPairs = cms.Sequence( allElecMuPairs )

# define additional collections of electron + muon candidates
# with loose track and ECAL isolation applied on electron leg
# (NOTE: to be used for the purpose of factorizing efficiencies
#        of electron isolation from other event selection criteria,
#        in order to avoid problems with limited Monte Carlo statistics)

allElecMuPairsLooseElectronIsolation = cms.EDProducer("PATElecMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatElectronsForElecMuTrkIPlooseIsolationCumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    dRmin12 = cms.double(-1.),
    srcMET = cms.InputTag('patMETs'),
    srcGenParticles = cms.InputTag('genParticles'),                                                        
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                                  
    verbosity = cms.untracked.int32(0)
)

produceElecMuPairsLooseElectronIsolation = cms.Sequence( allElecMuPairsLooseElectronIsolation )

produceElecMuPairsAll = cms.Sequence( produceElecMuPairs * produceElecMuPairsLooseElectronIsolation )


