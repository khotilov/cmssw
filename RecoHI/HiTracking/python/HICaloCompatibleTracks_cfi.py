import FWCore.ParameterSet.Config as cms

hiCaloCompatibleTracks  = cms.EDProducer("HICaloCompatibleTrackSelector",
                                         srcTracks = cms.InputTag("hiSelectedTracks"),
                                         srcTower = cms.InputTag("towerMaker"),
                                         srcPFCands = cms.InputTag("particleFlowTmp"),
                                         usePFCandMatching = cms.untracked.bool(True), 
                                         trkPtMin = cms.untracked.double(5.0),                                              
                                         trkEtaMax = cms.untracked.double(2.4),
                                         towerPtMin = cms.untracked.double(5.0),
                                         matchConeRadius = cms.untracked.double(0.087),
                                         # (calo energy sum/track pt) > caloCut (0 - loose, 1.0 - tight)
                                         caloCut = cms.untracked.double(0.3), 
                                         copyTrajectories = cms.untracked.bool(True), 
                                         copyExtras = cms.untracked.bool(True), ## set to false on AOD
                                         qualityToSet = cms.string("tight"),
                                         # root syntax, pt dependent calo-compatibility cut
                                         funcCaloComp = cms.string("0.2"), 
                                         # root syntax, pt dependent deltaR matching cut
                                         funcDeltaRTowerMatch = cms.string("0.087/(1.0+0.1*exp(-0.28*(x-20.)))") 
)
