import FWCore.ParameterSet.Config as cms

from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import patMuonPFIsolationSelector

#--------------------------------------------------------------------------------
# produce combinations of muon + muon pairs,
# the hypothesis being that the pair of muons results from a Z --> mu+ mu- decay
#--------------------------------------------------------------------------------

selectedPatMuonsForZmumuHypothesesMuonTrack = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon()'),
    filter = cms.bool(False)
)

selectedPatMuonsForZmumuHypothesesEta24 = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsForZmumuHypothesesMuonTrack"),
    cut = cms.string('abs(eta) < 2.4'),
    filter = cms.bool(False)
)

selectedPatMuonsForZmumuHypothesesPt10 = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsForZmumuHypothesesEta24"),
    cut = cms.string('pt > 15.'),
    filter = cms.bool(False)
)

selectedPatMuonsForZmumuHypothesesLoosePFRelIso = cms.EDFilter("PATMuonPFIsolationSelector",
    patMuonPFIsolationSelector.clone(
        sumPtMax = cms.double(0.15)
    ),
    src = cms.InputTag("selectedPatMuonsForZmumuHypothesesPt10"),                                                           
    filter = cms.bool(False)
)

selectedPatMuonsForZmumuHypotheses = cms.Sequence(
    selectedPatMuonsForZmumuHypothesesMuonTrack
   * selectedPatMuonsForZmumuHypothesesEta24
   * selectedPatMuonsForZmumuHypothesesPt10
   * selectedPatMuonsForZmumuHypothesesLoosePFRelIso
)    

allDiMuPairZmumuHypothesesByLooseIsolation = cms.EDProducer("PATDiMuPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('selectedPatMuonsForZmumuHypothesesLoosePFRelIso'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypothesesByLooseIsolation = cms.EDFilter("PATDiMuPairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypothesesByLooseIsolation"),                                   
    cut = cms.string('charge = 0'),
    filter = cms.bool(False)
)

produceDiMuPairsZmumuHypotheses = cms.Sequence(
    selectedPatMuonsForZmumuHypotheses
   * allDiMuPairZmumuHypothesesByLooseIsolation * selectedDiMuPairZmumuHypothesesByLooseIsolation
)
