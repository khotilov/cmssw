
import FWCore.ParameterSet.Config as cms

dqmElectronGeneralAnalysis = cms.EDAnalyzer("ElectronGeneralAnalyzer",

    Verbosity = cms.untracked.int32(0),
    FinalStep = cms.string("AtRunEnd"),
    InputFile = cms.string(""),
    OutputFile = cms.string(""),
    InputFolderName = cms.string("General"),
    OutputFolderName = cms.string("General"),
    
    ElectronCollection = cms.InputTag("gsfElectrons"),
    MatchingObjectCollection = cms.InputTag("mergedSuperClusters"),
    TrackCollection = cms.InputTag("generalTracks"),
    GsfTrackCollection = cms.InputTag("electronGsfTracks"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    
    TriggerResults = cms.InputTag("TriggerResults::HLT"),
    HltPaths = cms.vstring('HLT_Ele10_SW_L1R','HLT_Ele15_SW_L1R','HLT_Ele15_SW_EleId_L1R','HLT_Ele15_SW_LooseTrackIso_L1R','HLT_Ele15_SC15_SW_LooseTrackIso_L1R','HLT_Ele15_SC15_SW_EleId_L1R','HLT_Ele20_SW_L1R','HLT_Ele20_SC15_SW_L1R','HLT_Ele25_SW_L1R','HLT_Ele25_SW_EleId_LooseTrackIso_L1R','HLT_DoubleEle10_SW_L1R')

)


