import FWCore.ParameterSet.Config as cms

db_client = cms.EDAnalyzer("SiStripCommissioningOfflineDbClient",
  # general parameters
  FilePath         = cms.untracked.string('/tmp'),
  RunNumber        = cms.untracked.uint32(0),
  UseClientFile    = cms.untracked.bool(False),
  UploadHwConfig   = cms.untracked.bool(False),
  UploadAnalyses   = cms.untracked.bool(False),
  DisableDevices   = cms.untracked.bool(False),
  DisableBadStrips = cms.untracked.bool(False),
  SaveClientFile   = cms.untracked.bool(True),
  SummaryXmlFile   = cms.untracked.FileInPath('DQM/SiStripCommissioningClients/data/summary.xml'),
  # individual parameters
  ApvTimingParameters      = cms.PSet(
    TargetDelay = cms.int32(-1)        # -1: latest tick (old default), otherwise target delay for all ticks' rising edge
  ),
  CalibrationParameters    = cms.PSet(),
  DaqScopeModeParameters   = cms.PSet(),
  FastFedCablingParameters = cms.PSet(),
  FedCablingParameters     = cms.PSet(),
  FedTimingParameters      = cms.PSet(),
  FineDelayParameters      = cms.PSet(
    cosmic =  cms.bool(True)
  ),
  LatencyParamameters      = cms.PSet(
    OptimizePerPartition = cms.bool(False)
  ),
  NoiseParameters          = cms.PSet(),
  OptoScanParameters       = cms.PSet(
    TargetGain = cms.double(0.8),     # target gain (0.8 ~ 640ADC for tickmark)
    SkipGainUpdate = cms.bool(False)  # wether to keep the gain the same as already on the db
  ),
  PedestalsParameters      = cms.PSet(
    DeadStripMax  = cms.double(5),   # number times the noise spread below mean noise
    NoisyStripMin = cms.double(5),   # number times the noise spread above mean noise
    HighThreshold = cms.double(5),   # analysis-wide high threshold for the fed zero suppression
    LowThreshold  = cms.double(2)    # analysis-wide low threshold for the fed zero suppression
  ),
<<<<<<< OfflineDbClient_cff.py
  PedsOnlyParameters	   = cms.PSet(),
  PedsFullNoiseParameters  = cms.PSet(
    DeadStripMax  = cms.double(5),  # number times the noise spread below mean noise
    NoisyStripMin = cms.double(5),  # number times the noise spread above mean noise
    HighThreshold = cms.double(5),  # analysis-wide high threshold for the fed zero suppression
    LowThreshold  = cms.double(2)   # analysis-wide low threshold for the fed zero suppression
  ),
  SamplingParameters	   = cms.PSet(),
  VpspScanParameters	   = cms.PSet(),
=======
  PedsOnlyParameters       = cms.PSet(),
  SamplingParameters       = cms.PSet(),
  VpspScanParameters       = cms.PSet(),
>>>>>>> 1.4
)
