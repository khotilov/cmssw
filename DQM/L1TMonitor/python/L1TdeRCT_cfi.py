import FWCore.ParameterSet.Config as cms

l1tderct = cms.EDFilter("L1TdeRCT",
    rctSourceData = cms.InputTag("gctDigis"),
    HistFolder = cms.untracked.string('L1TEMU/L1TdeRCT/'),
    outputFile = cms.untracked.string('./L1TDQM.root'),
    verbose = cms.untracked.bool(False),
    DQMStore = cms.untracked.bool(True),
    singlechannelhistos = cms.untracked.bool(False),
    ecalTPGData = cms.InputTag("ecalEBunpacker","EcalTriggerPrimitives"),
    rctSourceEmul = cms.InputTag("valRctDigis"),
    disableROOToutput = cms.untracked.bool(True),
    hcalTPGData = cms.InputTag("hcalDigis"),
    gtDigisLabel = cms.InputTag("gtDigis"),
    gtEGAlgoName = cms.string("L1_SingleEG1"),
    doubleThreshold = cms.int32(3),

)

