import FWCore.ParameterSet.Config as cms

from DQMOffline.JetMET.jetMETAnalyzer_cff import *
from DQMOffline.JetMET.caloTowers_cff     import *
from DQMOffline.JetMET.BeamHaloAnalyzer_cfi import *

AnalyzeBeamHalo.StandardDQM = cms.bool(True)

jetMETDQMOfflineSource = cms.Sequence(analyzecaloTowersDQM*AnalyzeBeamHalo*jetMETAnalyzerSequence)
#jetMETDQMOfflineSource = cms.Sequence(analyzecaloTowersDQM*jetMETAnalyzerSequence)
