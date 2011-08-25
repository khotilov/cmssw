import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('prodBgEstTemplateProductionZtoMuTau')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START311_V2::All')

#--------------------------------------------------------------------------------
# produce templates for W + jets enriched control region without tau id. cuts applied
# (TaNC vloose passed && TaNC medium failed) and events weighted by fake-rates instead
#
# NOTE: fake-rates need to be added to pat::Tau producer before 
#       it gets modified by other functions
#
# import pat::Tau producer module
from TauAnalysis.RecoTools.patPFTauConfig_cfi import patTaus

# configure conditions for accessing fake-rate weight values from k-NearestNeighbour trees
#from TauAnalysis.BgEstimationTools.fakeRateConfiguration_cfi import setupFakeRates
#setupFakeRates(process, patTaus, frType = 'backgroundEstimation')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleZtoMuTauSpecific_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")
process.load("TauAnalysis.RecoTools.filterDataQuality_cfi")

# import sample definitions for running on the grid
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import *

#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveZtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('bgEstTemplatesZtoMuTau.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(1000)    
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        #'file:/data1/veelken/CMSSW_3_6_x/skims/Ztautau_1_1_sXK.root'
        #'file:/data1/veelken/CMSSW_3_8_x/skims/AHtoMuTau/selEvents_AHtoMuTau_woBtag_runs145762to148058_RECO.root'
        #'file:/data1/veelken/CMSSW_3_8_x/skims/test/mcDYttPU156bx_GEN_SIM_RECO_1_1_1VV.root'
        #'/store/relval/CMSSW_4_1_3/RelValZTT/GEN-SIM-RECO/START311_V2-v1/0037/286D4A6C-C651-E011-B6AC-001A92971BBA.root'
		'/store/user/jkolb/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola/skimMuTau_413_v1/83296e62d56aa89124dd9002fb3604c1/muTauSkim_11_1_c8K.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__setattr(process, "genPhaseSpaceCut", copy.deepcopy(#genPhaseSpaceCut#))
#__process.saveZtoMuTauPlots.outputFileName = #plotsOutputFileName#
#__#isBatchMode#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for configuring PAT trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process, hltProcess = 'HLT', outputModule = '')
process.patTrigger.addL1Algos = cms.bool(True)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
#switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)

# comment-out to take new HPS + TaNC combined tau id. algorithm
##switchToPFTauHPSpTaNC(process)
switchToPFTauHPS(process)

# disable preselection on of pat::Taus
# (disabled also in TauAnalysis/RecoTools/python/patPFTauConfig_cfi.py ,
#  but re-enabled after switching tau collection)
process.cleanPatTaus.preselection = cms.string('')

# add "ewkTauId" flag
#setattr(process.patTaus.tauIDSources, "ewkTauId", cms.InputTag('ewkTauId'))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment to replace caloJets by pfJets
switchJetCollection(process, jetCollection = cms.InputTag("ak5PFJets"), doBTagging = True, outputModule = '')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# set Boolean swich to true in order to apply type-1 corrections
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('patMETs'), cms.InputTag('patPFMETs'))
#--------------------------------------------------------------------------------

process.allMuTauPairs.doSVreco = cms.bool(False)
process.allMuTauPairsPFtype1MET.doSVreco = cms.bool(False)
process.allMuTauPairsLooseMuonIsolation.doSVreco = cms.bool(False)

#--------------------------------------------------------------------------------
# define analysis sequences for background enriched selections
process.load("TauAnalysis.Configuration.analyzeZtoMuTau_cff")
from TauAnalysis.Configuration.tools.analysisSequenceTools import addAnalyzer, addSysAnalyzer, removeAnalyzer
from TauAnalysis.CandidateTools.sysErrDefinitions_cfi import *

process.analyzeZtoMuTauEventsOS.name = cms.string('BgEstTemplateAnalyzer_Ztautau')
process.analyzeZtoMuTauEventsOS.eventDumps = cms.VPSet()
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'genPhaseSpaceEventInfoHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'eventWeightHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'muonHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'tauHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'diTauCandidateHistManagerForMuTau')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'diTauCandidateZmumuHypothesisHistManagerForMuTau')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'muPairHistManagerByLooseIsolation')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'jetHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'caloMEtHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'pfMEtHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'particleMultiplicityHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'vertexHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'triggerHistManagerForMuTau')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'dataBinner')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'diTauCandidateNSVfitHistManagerForMuTau')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'diTauCandidateNSVfitVtxMultiplicityBinGridHistManager')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'diTauCandidateNSVfitHistManagerForMuTauPFtype1MET')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'modelBinnerForMuTauGenTauLeptonPairAcc3mZbins')
removeAnalyzer(process.analyzeZtoMuTauEventsOS.analysisSequence, 'modelBinnerForMuTauWrtGenTauLeptonPairAcc3mZbins')
process.diTauCandidateHistManagerForMuTau.diTauCandidateSource = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative')
process.diTauCandidateHistManagerForMuTau.visMassHypothesisSource = cms.InputTag('')
addAnalyzer(process.analyzeZtoMuTauEventsOS, process.diTauCandidateHistManagerForMuTau, 'evtSelDiMuPairZmumuHypothesisVeto')
#process.diTauCandidateNSVfitHistManager.diTauCandidateSource = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative')
#process.diTauCandidateNSVfitHistManager.nSVfitEventHypotheses = cms.PSet(
#    psKine_MEt_logM_fit = cms.string('psKine_MEt_logM_fit')
#)

process.sysUncertaintyHistManager = cms.PSet(
    pluginName = cms.string('sysUncertaintyHistManager'),
    pluginType = cms.string('SysUncertaintyHistManager'),
    histManagers = cms.VPSet(
        cms.PSet(
            config = process.diTauCandidateHistManagerForMuTau,
            systematics = cms.PSet(
                diTauCandidateSource = getSysUncertaintyParameterSets(
                    [ muTauPairSystematics, ]
                )
            )
        )
    ),
    dqmDirectory_store = cms.string('sysUncertaintyHistManagerResults')
)
addSysAnalyzer(process.analyzeZtoMuTauEventsOS, process.sysUncertaintyHistManager, 'evtSelDiMuPairZmumuHypothesisVeto')
    
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauWplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauTTplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauZmumuEnrichedSelection_cff')
#process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauZbbEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauQCDenrichedSelection_cff')

# clone W + jets enriched analysis sequence,
# disable tau id. discriminator and add fake-rate weights
#process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauWplusJetsEnrichedSelection_frWeighted_cff')
#from TauAnalysis.BgEstimationTools.tools.fakeRateAnalysisTools import configureFakeRateWeightProduction
#configureFakeRateWeightProduction(process, method = "simple", 
#	                          preselPFTauJetSource = 'hpsTancTaus', frSet = 'bgEstTemplateTaNCinverted',
#                                  preselPatTauSource = 'selectedPatTausForMuTauLeadTrkPtCumulative',
#		                  addPatTauPreselection = 'tauID("againstMuonTight") > 0.5')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# switch to pfMET in all analysis sequences 
if hasattr(process, "patPFMETs"):
    process.muTauPairsBgEstQCDenriched.srcMET = cms.InputTag('patPFMETs')
    process.muTauPairsBgEstWplusJetsEnriched.srcMET = cms.InputTag('patPFMETs')
    #process.muTauPairsBgEstWplusJetsEnrichedFRweighted.srcMET = cms.InputTag('patPFMETs')
    process.muTauPairsBgEstTTplusJetsEnriched.srcMET = cms.InputTag('patPFMETs')
    process.muTauPairsBgEstZmumuJetMisIdEnriched.srcMET = cms.InputTag('patPFMETs')
    process.muTauPairsBgEstZmumuMuonMisIdEnriched.srcMET = cms.InputTag('patPFMETs')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# change upper limit on tranverse impact parameter of muon track to 2mm
changeCut(process, "selectedPatMuonsTrkIP", 0.2, attribute = "IpMax")

# switch between TaNC and HPS tau id. discriminators
changeCut(process, "selectedPatTausLeadTrk", 'tauID("decayModeFinding") > 0.5')
changeCut(process, "selectedPatTausLeadTrkPt", 'tauID("decayModeFinding") > 0.5')
changeCut(process, "selectedPatTausForMuTauLeadTrk", 'tauID("decayModeFinding") > 0.5')
changeCut(process, "selectedPatTausForMuTauLeadTrkPt", 'tauID("decayModeFinding") > 0.5')
changeCut(process, "selectedPatTausTaNCdiscr", 'tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5')
changeCut(process, "selectedPatTausForMuTauTaNCdiscr", 'tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5')                         ##  <<<<--------------- TAU ID
changeCut(process, "tausBgEstQCDenrichedTaNCdiscr", 'tauID("byVLooseIsolation") > 0.5 & tauID("byMediumIsolation") < 0.5')
changeCut(process, "tausBgEstWplusJetsEnrichedTaNCdiscr", 'tauID("byVLooseIsolation") > 0.5 & tauID("byMediumIsolation") < 0.5')
#changeCut(process, "tausBgEstWplusJetsEnrichedFRweightedTaNCdiscrNotApplied", "tauID('byHPSvloose') > -1000. & tauID('byHPSmedium') < +1000.")
changeCut(process, "tausBgEstZmumuJetMisIdEnrichedTaNCdiscr", 'tauID("byVLooseIsolation") > 0.5 & tauID("byMediumIsolation") < 0.5')

# correct muon veto
changeCut(process, "selectedPatTausForMuTauCaloMuonVeto", 'tauID("againstMuonTight") > -1.')
changeCut(process, "selectedPatTausForMuTauMuonVeto", 'tauID("againstMuonTight") > 0.5')

# change lower limit on separation required between muon and tau-jet to dR > 0.5
changeCut(process, "selectedMuTauPairsAntiOverlapVeto", "dR12 > 0.5")
changeCut(process, "selectedMuTauPairsAntiOverlapVetoLooseMuonIsolation", "dR12 > 0.5")
changeCut(process, "muTauPairsBgEstQCDenriched", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsBgEstWplusJetsEnriched", 0.5, attribute = "dRmin12")
#changeCut(process, "muTauPairsBgEstWplusJetsEnrichedFRweighted", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsBgEstTTplusJetsEnriched", 0.5, attribute = "dRmin12") 
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnriched", 0.5, attribute = "dRmin12") 
changeCut(process, "muTauPairsBgEstZmumuMuonMisIdEnriched", 0.5, attribute = "dRmin12") 

# change upper limit on muon + MET transverse mass to 40 GeV
changeCut(process, "selectedMuTauPairsMt1MET", "mt1MET < 40.")
changeCut(process, "selectedMuTauPairsMt1METlooseMuonIsolation", "mt1MET < 40.")
changeCut(process, "muTauPairsBgEstQCDenrichedMt1MET", "mt1MET < 40.")
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnrichedMt1MET", "mt1MET < 40.")

# disable cut on Pzeta variable
changeCut(process, "selectedMuTauPairsPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "selectedMuTauPairsPzetaDiffLooseMuonIsolation", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "muTauPairsBgEstQCDenrichedPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnrichedPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")

# disable b-tagging for now
#changeCut(process, "jetsBgEstTTplusJetsEnrichedEt40bTag", "bDiscriminator('trackCountingHighEffBJetTags') > -1000.")
#--------------------------------------------------------------------------------

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents
  + process.analyzeZtoMuTauEventsOS
  + process.bgEstWplusJetsEnrichedAnalysisSequence
  #+ process.bgEstWplusJetsEnrichedFRweightedAnalysisSequence
  + process.bgEstTTplusJetsEnrichedAnalysisSequence
  + process.bgEstZmumuEnrichedAnalysisSequence
  #+ process.bgEstZbbEnrichedAnalysisSequence  
  + process.bgEstQCDenrichedAnalysisSequence 
  + process.saveZtoMuTauPlots
)

process.q = cms.Path(process.dataQualityFilters)

process.schedule = cms.Schedule(process.q, process.p)

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True)
)

#--------------------------------------------------------------------------------
# import utility function for applyting Z-recoil corrections to MET
from TauAnalysis.Configuration.tools.mcToDataCorrectionTools import applyZrecoilCorrection_runZtoMuTau_bgEstTemplate
##applyZrecoilCorrection_runZtoMuTau_bgEstTemplate(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# disable accessing generator level information
# if running on data
#from TauAnalysis.Configuration.tools.switchToData import switchToData
#switchToData(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
process.producePatTupleAll = cms.Sequence(process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#
# define "hook" for enabling/disabling production of PAT-tuple event content,
# depending on whether RECO/AOD or PAT-tuples are used as input for analysis
#
#__#patTupleProduction#
if not hasattr(process, "isBatchMode"):
    process.p.replace(process.producePatTupleZtoMuTauSpecific, process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# update InputTag of TriggerResults
from TauAnalysis.Configuration.cfgOptionMethods import _setTriggerProcess
_setTriggerProcess(process, cms.InputTag('TriggerResults::HLT'))
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#print process.dumpPython()

