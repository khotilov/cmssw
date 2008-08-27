import FWCore.ParameterSet.Config as cms

from DQM.L1TMonitor.L1TGT_readout_cff import *
from DQM.L1TMonitor.L1TGCT_readout_cff import *
l1tmonitor = cms.Sequence(l1tgt*l1tgmt*l1trpctf*l1tcsctf*l1tdttf*l1trct*l1tgct)
l1tgt.gtSource = 'gtDigis::'
l1tgmt.gmtSource = 'gtDigis::'
l1tdttf.dttpgSource = 'dttfDigis::'
l1tcsctf.statusProducer = 'csctfDigis::'
l1tcsctf.lctProducer = 'csctfDigis::'
l1tcsctf.trackProducer = 'csctfDigis::'
l1trpctf.rpctfSource = 'gtDigis::'
l1tgct.gctCentralJetsSource = 'gctDigis:cenJets:'
l1tgct.gctForwardJetsSource = 'gctDigis:forJets:'
l1tgct.gctTauJetsSource = 'gctDigis:tauJets:'
l1tgct.gctEnergySumsSource = 'gctDigis::'
l1tgct.gctIsoEmSource = 'gctDigis:isoEm:'
l1tgct.gctNonIsoEmSource = 'gctDigis:nonIsoEm:'
l1trct.rctSource = 'gctDigis::'

from DQMOffline.Trigger.FourVectorHLTOffline_cfi import *
from DQMOffline.Trigger.Tau.HLTTauDQMOffline_cff import *
from DQMOffline.Trigger.EgammaHLTOffline_cfi import *
from Geometry.CaloEventSetup.CaloTopology_cfi import *
from DQM.L1TMonitor.L1TDEMON_cfi import *
l1temumonitor = cms.Sequence(l1demon) 

triggerOfflineDQMSource = cms.Sequence(l1temumonitor*l1tmonitor*hltResults*egammaHLTDQM*HLTTauDQMOffline)

