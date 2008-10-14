import FWCore.ParameterSet.Config as cms

#Analyzer taken from online dqm
from DQM.TrackingMonitor.MonitorTrackSTAMuons_cfi import *
from DQM.TrackingMonitor.MonitorTrackGLBMuons_cfi import *
from DQMOffline.Muon.dtSegmTask_cfi import *
from DQMOffline.Muon.rpcSource_cfi import *

#dedicated analyzers for offline dqm 
from DQMOffline.Muon.muonAnalyzer_cff import *
from DQMOffline.Muon.CSCMonitor_cfi import *
#dedicated clients for offline dqm 
from DQMOffline.Muon.muonQualityTests_cff import *

muonTrackAnalyzers = cms.Sequence(MonitorTrackSTAMuons*MonitorTrackGLBMuons)
muonMonitors = cms.Sequence(rpcSource*muonTrackAnalyzers*dtSegmentsMonitor*muonAnalyzer)
muonMonitorsAndQualityTests = cms.Sequence(muonMonitors*muonQualityTests)


