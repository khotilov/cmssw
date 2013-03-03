import FWCore.ParameterSet.Config as cms

process = cms.Process('TOPDQM')

## imports of standard configurations
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## --------------------------------------------------------------------
## Frontier Conditions: (adjust accordingly!!!)
##
## For CMSSW_3_8_X MC use             ---> 'START38_V12::All'
## For Data (38X re-processing) use   ---> 'GR_R_38X_V13::All'
## For Data (38X prompt reco) use     ---> 'GR10_P_V10::All'
##
## For more details have a look at: WGuideFrontierConditions
## --------------------------------------------------------------------
##process.GlobalTag.globaltag = 'GR_R_42_V14::All' 
process.GlobalTag.globaltag = 'START61_V8::All'


## input file(s) for testing
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring("file:input.root")
    fileNames = cms.untracked.vstring("/store/relval/CMSSW_6_2_0_pre1-START61_V8/RelValTTbarLepton/GEN-SIM-RECO/v1/00000/C6CC53CC-6E6D-E211-8EAB-003048D3756A.root")
)

## number of events
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

## apply VBTF electronID (needed for the current implementation
## of topSingleElectronDQMLoose and topSingleElectronDQMMedium)
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("DQM.Physics.topElectronID_cff")
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJets")
process.p = cms.Path(process.ak5JetTracksAssociatorAtVertex*process.btagging)


## output
process.output = cms.OutputModule("PoolOutputModule",
  fileName       = cms.untracked.string('topDQM_production.root'),
  outputCommands = cms.untracked.vstring(
    'drop *_*_*_*',
    'keep *_*_*_TOPDQM',
    'drop *_TriggerResults_*_TOPDQM',
    'drop *_simpleEleId70cIso_*_TOPDQM',
    'drop *_ak5JetTracksAssociatorAtVertex_*_TOPDQM',
    'drop *_btagging_*_TOPDQM',
    'drop *_jetProbabilityBJetTags_*_TOPDQM',
    'drop *_ghostTrackBJetTags_*_TOPDQM',
    'drop *_combinedSecondaryVertexMVABJetTags_*_TOPDQM',
    'drop *_trackCountingHighPurBJetTags_*_TOPDQM',
    'drop *_trackCountingHighEffBJetTags_*_TOPDQM',
    'drop *_simpleSecondaryVertexHighEffBJetTags_*_TOPDQM',
    'drop *_simpleSecondaryVertexHighPurBJetTags_*_TOPDQM',
    'drop *_softElectronByIP3dBJetTags_*_TOPDQM',
    'drop *_softElectronByPtBJetTags_*_TOPDQM',
    'drop *_softMuonBJetTags_*_TOPDQM',
    'drop *_softMuonByIP3dBJetTags_*_TOPDQM',
    'drop *_impactParameterTagInfos_*_TOPDQM',
    'drop *_combinedSecondaryVertexBJetTags_*_TOPDQM',
    'drop *_softMuonByPtBJetTags_*_TOPDQM',
    'drop *_ghostTrackVertexTagInfos_*_TOPDQM',
    'drop *_secondaryVertexTagInfos_*_TOPDQM',
    'drop *_softElectronCands_*_TOPDQM',
    'drop *_softMuonTagInfos_*_TOPDQM',
  ),
  splitLevel     = cms.untracked.int32(0),
  dataset = cms.untracked.PSet(
    dataTier   = cms.untracked.string(''),
    filterName = cms.untracked.string('')
  )
)

## load jet corrections
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
#process.prefer("ak5CaloL2L3")
process.prefer("ak5PFL2L3")

## check the event content
process.content = cms.EDAnalyzer("EventContentAnalyzer")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('TopSingleLeptonDQM'   )
process.MessageLogger.cerr.TopSingleLeptonDQM    = cms.untracked.PSet(limit = cms.untracked.int32(1))
process.MessageLogger.categories.append('TopDiLeptonOfflineDQM')
process.MessageLogger.cerr.TopDiLeptonOfflineDQM = cms.untracked.PSet(limit = cms.untracked.int32(1))
process.MessageLogger.categories.append('SingleTopTChannelLeptonDQM'   )
process.MessageLogger.cerr.SingleTopTChannelLeptonDQM    = cms.untracked.PSet(limit = cms.untracked.int32(1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MEtoEDMConverter.deleteAfterCopy = cms.untracked.bool(False)  ## line added to avoid crash when changing run number


## path definitions
process.p      = cms.Path(
    process.ak5JetTracksAssociatorAtVertex *
    process.btagging *
    process.simpleEleId70cIso          *
    process.DiMuonDQM                  +
    process.DiElectronDQM              +
    process.ElecMuonDQM                +
    process.topSingleMuonLooseDQM      +
    process.topSingleMuonMediumDQM     +
    process.topSingleElectronLooseDQM  +
    process.topSingleElectronMediumDQM +
    process.singleTopMuonMediumDQM     +
    process.singleTopElectronMediumDQM
)
process.endjob = cms.Path(
    process.endOfProcess
)
process.fanout = cms.EndPath(
    process.output
)

## schedule definition
process.schedule = cms.Schedule(
    process.p,
    process.endjob,
    process.fanout
)
