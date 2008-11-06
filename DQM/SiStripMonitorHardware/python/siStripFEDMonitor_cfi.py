import FWCore.ParameterSet.Config as cms

siStripFEDMonitor = cms.EDAnalyzer("SiStripFEDMonitorPlugin",
  #Raw data collection
  RawDataTag = cms.untracked.InputTag('source'),
  
  #Dump buffer info and raw data if any error is found
  PrintDebugMessages = cms.untracked.bool(False),
  #Write the DQM store to a file (DQMStore.root) at the end of the run
  WriteDQMStore = cms.untracked.bool(False),
  
  #Do not book expert histograms at global level unless PreBookAllHistos is set
  DisableGlobalExpertHistograms = cms.untracked.bool(False),
  #Disable the FED level histograms
  DisableFEDHistograms = cms.untracked.bool(False),
  #Override previous two option and book and fill all histograms (so that files can be merged)
  FillAllHistograms = cms.untracked.bool(False)
)
