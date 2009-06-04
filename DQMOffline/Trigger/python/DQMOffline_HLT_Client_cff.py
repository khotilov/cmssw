import FWCore.ParameterSet.Config as cms

from DQM.HLTEvF.HLTEventInfoClient_cfi import *

from DQMOffline.Trigger.FourVectorHLTOfflineClient_cff import *
from DQMOffline.Trigger.EgHLTOfflineClient_cfi import *
from DQMOffline.Trigger.MuonPostProcessor_cfi import *
from DQMOffline.Trigger.JetMETHLTOfflineClient_cfi import *

#    # use include file for dqmEnv 
from DQMServices.Components.DQMEnvironment_cfi import *
dqmEnv.subSystemFolder = 'HLTOffline'

hltOfflineDQMClient = cms.Sequence(hltFourVectorSeqClient*egHLTOffDQMClient*hltEventInfoClient*HLTMuonPostVal*jetMETHLTOfflineClient*dqmEnv)
