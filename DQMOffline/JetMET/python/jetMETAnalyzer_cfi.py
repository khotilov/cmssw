import FWCore.ParameterSet.Config as cms

jetMETAnalyzer = cms.EDAnalyzer("JetMETAnalyzer",
    OutputMEsInRootFile = cms.bool(False),
    PFJetsCollectionLabel = cms.InputTag("iterativeCone5PFJets"),
    caloMETAnalysis = cms.PSet(
    HLTPathsJetMB = cms.vstring("HLT_L1Jet15",
	                    "HLT_Jet30", 
	                    "HLT_Jet50", 
                            "HLT_Jet80", 
	                    "HLT_Jet110", 
	                    "HLT_Jet180",
                            "HLT_DiJetAve15",
                            "HLT_DiJetAve30",
                            "HLT_DiJetAve50",
                            "HLT_DiJetAve70",
                            "HLT_DiJetAve130",
                            "HLT_DiJetAve220",
	                    "HLT_ZeroBias",
                            "HLT_MinBias"),
    ),
    OutputFileName = cms.string('jetMETMonitoring.root'),
    jetAnalysis = cms.PSet(
        eBin    = cms.int32(100),
        phiMin  = cms.double(-3.2),
        ptBin   = cms.int32(100),
        eMin    = cms.double(0.0),
        eMax    = cms.double(500.0),
        pMin    = cms.double(0.0),
        etaBin  = cms.int32(100),
        etaMin  = cms.double(-5.0),
        ptMin   = cms.double(0.0),
        phiBin  = cms.int32(70),
        pBin    = cms.int32(100),
        ptMax   = cms.double(50.0),
        etaMax  = cms.double(5.0),
        pMax    = cms.double(500.0),
        phiMax  = cms.double(3.2)
    ),
    pfJetAnalysis = cms.PSet(

    ),
    DoPFJetAnalysis            = cms.untracked.bool(True),
    SCJetsCollectionLabel      = cms.InputTag("sisCone5CaloJets"),
    DoCaloMETAnalysis          = cms.untracked.bool(True),
    DoJetAnalysis              = cms.untracked.bool(True),
    CaloMETCollectionLabel     = cms.InputTag("met"),
    CaloMETNoHFCollectionLabel = cms.InputTag("metNoHF"),
    ICJetsCollectionLabel      = cms.InputTag("iterativeCone5CaloJets"),
    TriggerResultsLabel        = cms.InputTag("TriggerResults::HLT")
)


