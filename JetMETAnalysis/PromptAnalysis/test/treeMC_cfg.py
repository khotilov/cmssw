import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process('TEST')
process.load('JetMETAnalysis.PromptAnalysis.ntuple_cff')

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")

###########
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag ='GR09_R_34X_V5::All'##make sure to check the GT from DBS
###########

# process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
# process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
# process.load("L1Trigger/Configuration/L1RawToDigi_cff")
# process.load("RecoMET/Configuration/RecoMET_BeamHaloId_cff")
# process.load("RecoMET/METProducers/BeamHaloSummary_cfi")
# process.load("RecoMET/METProducers/CSCHaloData_cfi")
# process.load("RecoMET/METProducers/EcalHaloData_cfi")
# process.load("RecoMET/METProducers/HcalHaloData_cfi")
# process.load("RecoMET/METProducers/GlobalHaloData_cfi")

process.load("Configuration/StandardSequences/ReconstructionCosmics_cff")

process.load("RecoMuon/Configuration/RecoMuon_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
fileName = cms.string("MinimumBias__BeamCommissioning09-BSCNOBEAMHALO-Jan23Skim-v1__RAW-RECO.root"), ##give a name to the output file
                           closeFileFast = cms.untracked.bool(True)  ) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring(
    "/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/BSCNOBEAMHALO-Dec14thSkim_v1/0102/BABAF8C3-71EA-DE11-9D8C-0024E8768446.root"
    #'/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8D_900GeV-v1/0005/E4590360-4CD7-DE11-8CB4-002618943896.root'
    #'/store/data/BeamCommissioning09/MinimumBias/RECO/rereco_GR09_P_V7_v1/0099/DABD5D6D-D4E2-DE11-8FFD-00261894387A.root'
    #"/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/123/596/F82DED93-36E2-DE11-9316-000423D9870C.root"
    #"/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/123/151/0E45A7CE-F5DD-DE11-9B2E-001617E30CC8.root"
    ),
    
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    
    secondaryFileNames = cms.untracked.vstring())

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100

# jet corrections
process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_ReReco332_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_900GeV_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_2360GeV_cff")

# summary
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
# from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
# process.bit40OR41 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('40 OR 41'))

# from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
# process.physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])

process.promptanaTree = cms.EDAnalyzer("PromptAnaTree",
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_promptanaevent_*_*',
    'keep *_promptanamet_*_*',
    'keep *_promptanatcmet_*_*',
    'keep *_promptanapfmet_*_*',
    'keep *_promptananohf_*_*',
    'keep *_promptanaic5calojet_*_*',
    #'keep *_promptanasc5calojet_*_*',
    #'keep *_promptanakt4calojet_*_*',
    'keep *_promptanaak5calojet_*_*',
    'keep *_promptanaJPTak5_*_*',
    'keep *_promptanaak5pfjet_*_*',
#    'keep *_promptanahalo_*_*',
    'keep *_promptanacalotowers_*_*',
    'keep *_promptanatrigger_*_*',
    'keep *_promptanavtx_*_*',
    'keep *_promptanatrack_*_*',
    'keep *_promptanacleanup_*_*',
    'keep *_promptanaecalspikes_*_*',
    'keep *_promptanaPMTnoise_*_*'
    ))

process.theBigNtuple = cms.Path(
    #process.physDecl *
    #process.bit40OR41 *
    #process.BeamHaloId *
    (
    process.promptanaevent +
    process.promptanamet   +
    process.promptanatcmet   +
    process.promptanapfmet   +
    process.promptananohf  +
    process.promptanaic5calojet +
    #process.promptanasc5calojet +
    #process.promptanakt4calojet +
    process.promptanaak5calojet +
    process.promptanaJPTak5 +
    process.promptanaak5pfjet +
    #process.promptanahalo +
    process.promptanacalotowers +
    process.promptanatrigger +
    process.promptanavtx +
    process.promptanatrack +
    process.promptanacleanup +
    process.promptanaecalspikes+
    process.promptanaPMTnoise
    )
    * process.promptanaTree )

