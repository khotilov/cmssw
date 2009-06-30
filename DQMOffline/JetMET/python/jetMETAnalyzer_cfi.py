import FWCore.ParameterSet.Config as cms
#from JetMETCorrections.Configuration.JetPlusTrackCorrections_cff import *
#from JetMETCorrections.Configuration.ZSPJetCorrections152_cff import *

jetMETAnalyzer = cms.EDAnalyzer("JetMETAnalyzer",

    #
    # Output files
    #
    OutputMEsInRootFile = cms.bool(False),
    OutputFileName = cms.string('jetMETMonitoring.root'),

    #
    #
    #
    TriggerResultsLabel        = cms.InputTag("TriggerResults::HLT"),
    JetLo                      = cms.string("HLT_Jet30"),
    JetHi                      = cms.string("HLT_Jet110"),

    #
    # Jet-related
    #                                                                   
    DoPFJetAnalysis            = cms.untracked.bool(True),
    DoJPTJetAnalysis           = cms.untracked.bool(True),
    DoJetAnalysis              = cms.untracked.bool(True),
    PFJetsCollectionLabel      = cms.InputTag("iterativeCone5PFJets"),
    JPTJetsCollectionLabel     = cms.InputTag("JetPlusTrackZSPCorJetIcone5"),
    SCJetsCollectionLabel      = cms.InputTag("sisCone5CaloJets"),
    ICJetsCollectionLabel      = cms.InputTag("iterativeCone5CaloJets"),

    #
    # MET-related
    #                                                                   
    DoCaloMETAnalysis            = cms.untracked.bool(True),
    DoTcMETAnalysis              = cms.untracked.bool(True),
    DoPfMETAnalysis              = cms.untracked.bool(True),
    DoHTMHTAnalysis              = cms.untracked.bool(True),
    CaloMETCollectionLabel       = cms.InputTag("met"),
    CaloMETNoHFCollectionLabel   = cms.InputTag("metNoHF"),
    CaloMETHOCollectionLabel     = cms.InputTag("metHO"),
    CaloMETNoHFHOCollectionLabel = cms.InputTag("metNoHFHO"),
    TcMETCollectionLabel         = cms.InputTag("tcMet"),
    PfMETCollectionLabel         = cms.InputTag("pfMet"),
    JetCollectionForHTMHTLabel   = cms.InputTag("iterativeCone5CaloJets"),

    #
    # For caloMETAnalysis
    #
    caloMETAnalysisMain = cms.PSet(
    HLTPathsJetMB = cms.vstring(),
#   When it is empty, it accepts all the triggers
#    HLTPathsJetMB = cms.vstring("HLT_L1Jet15",
#                            "HLT_Jet30",
#                            "HLT_Jet50",
#                            "HLT_Jet80",
#                            "HLT_Jet110",
#                            "HLT_Jet180",
#                            "HLT_DiJetAve15",
#                            "HLT_DiJetAve30",
#                            "HLT_DiJetAve50",
#                            "HLT_DiJetAve70",
#                            "HLT_DiJetAve130",
#                            "HLT_DiJetAve220",
#                            "HLT_MinBias"),
    etThreshold = cms.double(1.),
    allHist     = cms.bool(True)
    ),

    #
    # For caloMETAnalysis
    #
    caloMETAnalysis = cms.PSet(
    HLTPathsJetMB = cms.vstring(),
    etThreshold = cms.double(1.),
    allHist     = cms.bool(False)
    ),

    #
    # For tcMETAnalysis
    #
    tcMETAnalysis = cms.PSet(
    HLTPathsJetMB = cms.vstring(),
    etThreshold = cms.double(1.)
    ),

    #
    # For pfMETAnalysis
    #
    pfMETAnalysis = cms.PSet(
    HLTPathsJetMB = cms.vstring(),
    etThreshold = cms.double(1.)
    ),

    #
    # For HTMHTAnalysis
    #
    HTMHTAnalysis = cms.PSet(
    HLTPathsJetMB = cms.vstring(),
    ptThreshold = cms.double(30.)
    ),


    #
    # For jetAnalysis
    #
    jetAnalysis = cms.PSet(
        ptThreshold = cms.double(3.),
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

    #
    # For PF jetAnalysis
    #
    pfJetAnalysis = cms.PSet(
        ptThreshold = cms.double(3.),
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

    #
    # For JPT jetAnalysis
    #
    JPTJetAnalysis = cms.PSet(
        ptThreshold = cms.double(3.),
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
    )


)
# no JPT :
jetMETAnalyzerSequence = cms.Sequence(jetMETAnalyzer)
# with JPT take this sequence:
#jetMETAnalyzer = cms.Sequence(ZSPJetCorrections*JetPlusTrackCorrections*jetMETAnalyze)
