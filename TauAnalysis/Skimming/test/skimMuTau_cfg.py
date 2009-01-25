import FWCore.ParameterSet.Config as cms

process = cms.Process("muTauSkim")

from TauAnalysis.Skimming.EventContent_cff import *

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
#
# Z --> tau tau (all decay modes; simulated with TAUOLA)
# 9k events CMSSW_2_2_3 RelVal sample
#
    '/store/relval/CMSSW_2_2_3/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v4/0003/A4A3988A-BCCB-DD11-A103-001617E30E28.root',
    '/store/relval/CMSSW_2_2_3/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v4/0003/D412FFFC-BCCB-DD11-8B20-000423D952C0.root',
    '/store/relval/CMSSW_2_2_3/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v4/0003/F01E4F34-BDCB-DD11-B87D-001617C3B77C.root',
    '/store/relval/CMSSW_2_2_3/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v4/0004/1CAA08F8-D3CB-DD11-ADF9-000423D6B358.root',
    '/store/relval/CMSSW_2_2_3/RelValZTT/GEN-SIM-RECO/STARTUP_V7_v4/0004/2800478C-08CC-DD11-94BB-0019B9F72BAA.root'
  )
)

#--------------------------------------------------------------------------------
# select muons and tau-jets
#--------------------------------------------------------------------------------

process.selectedMuons = cms.EDFilter("MuonSelector",
  src = cms.InputTag('muons'),
  cut = cms.string("pt > 10 & abs(eta) < 2.5"),
  filter = cms.bool(True)
)

process.selectedTaus = cms.EDFilter("PFTauSelector",
  src = cms.InputTag('pfRecoTauProducer'),
  discriminators = cms.VPSet(
    cms.PSet(
      discriminator = cms.InputTag("pfRecoTauDiscriminationByLeadingPionPtCut"),
      selectionCut = cms.double(0.5)
    )
  ),
  filter = cms.bool(True)
)

#--------------------------------------------------------------------------------
# combine selected muons and tau-jets into pairs
#--------------------------------------------------------------------------------

process.muTauPairs = cms.EDProducer("DiTauProducer",
  hadronicTaus = cms.InputTag('selectedTaus'),
  leptonicTaus = cms.InputTag('selectedMuons'),
  METs = cms.InputTag(''),
  metMode = cms.int32(1),
  useLeadingTaus = cms.bool(False),
  verbose =  cms.untracked.bool(False)
)

#--------------------------------------------------------------------------------
# discard tau-jets that pass muon identification criteria
# (in order to reject events in which there is only one muon and no tau-jet;
#  note that almost all muons get selected as tau-jets !!)
#--------------------------------------------------------------------------------

process.selectedMuTauPairs = cms.EDFilter("DiTauAntiOverlapSelector",
  src = cms.InputTag('muTauPairs'),
  dRmin = cms.double(0.7),
  filter = cms.bool(True)                                     
)

process.muTauSkimPath = cms.Path( (process.selectedTaus + process.selectedMuons) * process.muTauPairs * process.selectedMuTauPairs )

muTauEventSelection = cms.untracked.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('muTauSkimPath')
  )
)
process.muTauSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
#  AODSIMEventContent,
  tauAnalysisEventContent,                                               
  muTauEventSelection,
  fileName = cms.untracked.string('muTauSkim.root')
)
#process.muTauSkimOutputModule.outputCommands.extend(RecoEcalRECO.outputCommands)
#process.muTauSkimOutputModule.outputCommands.extend(RecoParticleFlowRECO.outputCommands)
#
# keep all ECAL + HCAL recHits for computation of IsoDeposits
#
#allEcalRecHits = cms.PSet(
#  outputCommands = cms.untracked.vstring('keep *EcalRecHit_*_*_*',
#                                         'keep *EcalRecHit_*_*_*')
#)
#process.muTauSkimOutputModule.outputCommands.extend(allEcalRecHits.outputCommands)
#allHcalRecHits = cms.PSet(
#  outputCommands = cms.untracked.vstring('keep *HBHERecHit_*_*_*',
#                                         'keep *HFRecHit_*_*_*',
#                                         'keep *HORecHit_*_*_*')
#)
#process.muTauSkimOutputModule.outputCommands.extend(allHcalRecHits.outputCommands)

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath( process.muTauSkimOutputModule )

