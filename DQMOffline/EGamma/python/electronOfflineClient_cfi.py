import FWCore.ParameterSet.Config as cms

from DQMOffline.EGamma.photonAnalyzer_cfi import *


dqmElectronOfflineClient = cms.EDAnalyzer("ElectronOfflineClient",

    Verbosity = cms.untracked.int32(0),
    FinalStep = cms.string("AtRunEnd"),
    InputFile = cms.string(""),
    OutputFile = cms.string(""),
    InputFolderName = cms.string(""),
    OutputFolderName = cms.string(""),
    
    EffHistoTitle = cms.string("fraction of reco ele matching a reco sc")

)
