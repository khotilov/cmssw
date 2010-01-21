import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.MuonTrigRateAnalyzer_cfi import *
from DQMOffline.Trigger.TopTrigRateAnalyzer_cfi import *
#from DQMOffline.Trigger.BPAGTrigRateAnalyzer_cfi import *
from DQMOffline.Trigger.TnPEfficiency_cff import *

muonFullOfflineDQM = cms.Sequence(offlineDQMMuonTrig
								  + topTrigOfflineDQM
								  + TnPEfficiency)

