import FWCore.ParameterSet.Config as cms

process = cms.Process("elecMuSkim")

from TauAnalysis.Skimming.EventContent_cff import *

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('MC_31X_V2::All')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/A4DD1FAE-B178-DE11-B608-001D09F24EAC.root',
        '/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/9408B54D-CB78-DE11-9AEB-001D09F2503C.root'
    )
)

#--------------------------------------------------------------------------------
# select electrons and muons
#--------------------------------------------------------------------------------

process.selectedElectrons = cms.EDFilter("GsfElectronSelector",
    src = cms.InputTag("gsfElectrons"),
    cut = cms.string("pt > 8 & abs(eta) < 2.5 & eSuperClusterOverP>0.8 & eSuperClusterOverP<1.25"),
    filter = cms.bool(True)
)

process.selectedMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string("pt > 8 & abs(eta) < 2.5"),
    filter = cms.bool(True)
)

process.elecMuPairs = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedElectrons'),
    srcLeg2 = cms.InputTag('selectedMuons'),
    dRmin12 = cms.double(0.),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                 
    verbosity = cms.untracked.int32(0)
)

process.selectedElecMuPairs = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag('elecMuPairs'),
    cut = cms.string("dR12 > 0.3"),
    filter = cms.bool(True)                                     
)

#--------------------------------------------------------------------------------
# keep event in case it passed the muon + electron selection
#--------------------------------------------------------------------------------

process.elecMuSkimPath = cms.Path(
    (process.selectedElectrons + process.selectedMuons)
   * process.elecMuPairs
   * process.selectedElecMuPairs
)

elecMuEventSelection = cms.untracked.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('elecMuSkimPath')
    )
)

#--------------------------------------------------------------------------------
#Dummy filter (always false) necessary to get correct efficiency report.
#--------------------------------------------------------------------------------
process.dummyFilter = cms.EDFilter("HLTHighLevel",
     TriggerResultsTag = cms.InputTag("TriggerResults","",process.name_()),
     HLTPaths = cms.vstring(''),# provide list of HLT paths (or patterns) you want
     eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
     andOr = cms.bool(True),   # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
     throw = cms.bool(True)    # throw exception on unknown path names
)

#--------------------------------------------------------------------------------
# fill validation histograms for events passing the electron + muon selection
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

from TauAnalysis.Skimming.ewkElecMuValHistManager_cfi import *

process.fillElecMuValPlots = cms.EDAnalyzer("EwkTauValidation",

    # list of individual channels                           
    channels = cms.VPSet(
        ewkElecMuValHistManager
    ),

    # disable all warnings
    maxNumWarnings = cms.int32(1)                       
)

process.saveElecMuValPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('elecMuValPlots.root')
)

process.p = cms.Path(
    process.fillElecMuValPlots
   + process.saveElecMuValPlots
   + process.dummyFilter
)

#--------------------------------------------------------------------------------
# save events passing the electron + muon selection
#--------------------------------------------------------------------------------

process.elecMuSkimOutputModule = cms.OutputModule("PoolOutputModule",                                 
    tauAnalysisEventContent,                                               
    elecMuEventSelection,
    fileName = cms.untracked.string('elecMuSkim.root')                                                  
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.o = cms.EndPath(process.elecMuSkimOutputModule)

