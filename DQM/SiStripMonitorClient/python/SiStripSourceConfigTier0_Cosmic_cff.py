import FWCore.ParameterSet.Config as cms

# FED integrity Check
from DQM.SiStripMonitorHardware.siStripFEDCheck_cfi import *
siStripFEDCheck.HistogramUpdateFrequency = 0
siStripFEDCheck.DoPayloadChecks          = True
siStripFEDCheck.CheckChannelLengths      = True
siStripFEDCheck.CheckChannelPacketCodes  = True
siStripFEDCheck.CheckFELengths           = True
siStripFEDCheck.CheckChannelStatus       = True

# FED Monitoring
from DQM.SiStripMonitorHardware.siStripFEDMonitor_Tier0_cff import *

# SiStripMonitorDigi ####
from DQM.SiStripMonitorDigi.SiStripMonitorDigi_cfi import *
SiStripMonitorDigi.Mod_On = False

# SiStripMonitorDigi ####
from DQM.SiStripMonitorCluster.SiStripMonitorCluster_cfi import *
SiStripMonitorCluster.Mod_On = False

# SiStripMonitorTrack ####
# Clone for Cosmic Tracks
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrack_cosmicTk  = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrack_cosmicTk.TrackProducer = 'cosmictrackfinderP5'
SiStripMonitorTrack_cosmicTk.Mod_On        = False
SiStripMonitorTrack_cosmicTk.FolderName    = 'SiStrip/Tracks'
# Clone for CKF Tracks
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrack_ckf = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrack_ckf.TrackProducer      = 'ctfWithMaterialTracksP5'
SiStripMonitorTrack_ckf.Mod_On             = False
SiStripMonitorTrack_ckf.FolderName         = 'SiStrip/Tracks'
# Clone for Road Search  Tracks
import DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi
SiStripMonitorTrack_rs = DQM.SiStripMonitorTrack.SiStripMonitorTrack_cfi.SiStripMonitorTrack.clone()
SiStripMonitorTrack_rs.TrackProducer       = 'rsWithMaterialTracksP5'
SiStripMonitorTrack_rs.Mod_On              = False
SiStripMonitorTrack_rs.FolderName          = 'SiStrip/Tracks'

# TrackerMonitorTrack ####
# Clone for Cosmic Track Finder
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResiduals_cosmicTk = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
MonitorTrackResiduals_cosmicTk.trajectoryInput     = 'cosmictrackfinderP5'
MonitorTrackResiduals_cosmicTk.OutputMEsInRootFile = False
MonitorTrackResiduals_cosmicTk.Mod_On              = False
# Clone for CKF Tracks
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResiduals_ckf = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
MonitorTrackResiduals_ckf.trajectoryInput          = 'ctfWithMaterialTracksP5'
MonitorTrackResiduals_ckf.OutputMEsInRootFile      = False
MonitorTrackResiduals_ckf.Mod_On                   = False
# Clone for Road Search  Tracks
import DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi
MonitorTrackResiduals_rs = DQM.TrackerMonitorTrack.MonitorTrackResiduals_cfi.MonitorTrackResiduals.clone()
MonitorTrackResiduals_rs.trajectoryInput           = 'rsWithMaterialTracksP5'
MonitorTrackResiduals_rs.OutputMEsInRootFile       = False
MonitorTrackResiduals_rs.Mod_On                    = False

# TrackingMonitor ####
# Clone for Cosmic Track Finder
import DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi
TrackMon_cosmicTk = DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi.TrackerCosmicTrackMon.clone()
TrackMon_cosmicTk.TrackProducer                    = 'cosmictrackfinderP5'
TrackMon_cosmicTk.AlgoName                         = 'CosmicTk'
TrackMon_cosmicTk.FolderName                       = 'SiStrip/Tracks'

# Clone for CKF Tracks
import DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi
TrackMon_ckf = DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi.TrackerCosmicTrackMon.clone()
TrackMon_ckf.TrackProducer                         = 'ctfWithMaterialTracksP5'
TrackMon_ckf.AlgoName                              = 'CKFTk'
TrackMon_ckf.FolderName                            = 'SiStrip/Tracks'

# Clone for Road Search  Tracks
import DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi
TrackMon_rs = DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi.TrackerCosmicTrackMon.clone()
TrackMon_rs.TrackProducer                          = 'rsWithMaterialTracksP5'
TrackMon_rs.AlgoName                               = 'RSTk'
TrackMon_rs.FolderName                             = 'SiStrip/Tracks'

# Clone for Beam Halo Muon Tracks
import DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi
TrackMon_bhmuon = DQM.TrackingMonitor.TrackerCosmicsTrackingMonitor_cfi.TrackerCosmicTrackMon.clone()
TrackMon_bhmuon.TrackProducer                      = 'ctfWithMaterialTracksBeamHaloMuon'
TrackMon_bhmuon.AlgoName                           = 'BHMuonTk'
TrackMon_bhmuon.FolderName                         = 'SiStrip/Tracks'

# Tracking Efficiency
# Clone for Cosmic Tracks
import DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi
TrackEffMon_cosmicTk = DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi.TrackEffMon.clone()
TrackEffMon_cosmicTk.TKTrackCollection             = 'cosmictrackfinderP5'
TrackEffMon_cosmicTk.AlgoName                      = 'CosmicTk'
TrackEffMon_cosmicTk.FolderName                    = 'SiStrip/Tracks/Efficiencies'

# Clone for CKF Tracks
import DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi
TrackEffMon_ckf = DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi.TrackEffMon.clone()
TrackEffMon_ckf.TKTrackCollection                  = 'ctfWithMaterialTracksP5'
TrackEffMon_ckf.AlgoName                           = 'CKFTk'
TrackEffMon_ckf.FolderName                         = 'SiStrip/Tracks/Efficiencies'

# Clone for RS Tracks
import DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi
TrackEffMon_rs = DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi.TrackEffMon.clone()
TrackEffMon_rs.TKTrackCollection                   = 'rsWithMaterialTracksP5'
TrackEffMon_rs.AlgoName                            = 'RSTk'
TrackEffMon_rs.FolderName                          = 'SiStrip/Tracks/Efficiencies'

# Clone for Beam Halo  Tracks
import DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi
TrackEffMon_bhmuon = DQM.TrackingMonitor.TrackEfficiencyMonitor_cfi.TrackEffMon.clone()
TrackEffMon_bhmuon.TKTrackCollection               = 'ctfWithMaterialTracksBeamHaloMuon'
TrackEffMon_bhmuon.AlgoName                        = 'BHMuonTk'
TrackEffMon_bhmuon.FolderName                      = 'SiStrip/Tracks/Efficiencies'

# DQM Services
dqmInfoSiStrip = cms.EDFilter("DQMEventInfo",
     subSystemFolder = cms.untracked.string('SiStrip')
)

# Services needed for TkHistoMap
TkDetMap = cms.Service("TkDetMap")
SiStripDetInfoFileReade = cms.Service("SiStripDetInfoFileReader")


# Sequences 
SiStripDQMTier0_cosmicTk = cms.Sequence(SiStripMonitorTrack_cosmicTk*MonitorTrackResiduals_cosmicTk*TrackMon_cosmicTk*TrackEffMon_cosmicTk)

SiStripDQMTier0_ckf = cms.Sequence(SiStripMonitorTrack_ckf*MonitorTrackResiduals_ckf*TrackMon_ckf*TrackEffMon_ckf)

SiStripDQMTier0_rs = cms.Sequence(SiStripMonitorTrack_rs*MonitorTrackResiduals_rs*TrackMon_rs*TrackEffMon_rs)

SiStripDQMTier0 = cms.Sequence(siStripFEDMonitor*SiStripMonitorDigi*SiStripMonitorCluster*SiStripMonitorTrack_ckf*MonitorTrackResiduals_ckf*TrackMon_cosmicTk*TrackMon_ckf*TrackMon_rs*TrackEffMon_ckf*dqmInfoSiStrip)
