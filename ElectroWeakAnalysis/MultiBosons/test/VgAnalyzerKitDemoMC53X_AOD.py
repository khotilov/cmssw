import FWCore.ParameterSet.Config as cms

process = cms.Process("VgKit")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('START53_V7A::All')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    #'file:/data3/ncuhep/52X_vgamma/testfiles/WGToLNuG_TuneZ2star_8TeV_PU_S7_START52_V9-v2.root'),
    'file:/data2/poter/MC_CMSSW537_vgamma/testfiles/DYJetsToLL_M-50_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root'),
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets25 = process.kt6PFJets.clone( doAreaFastjet=True, doRhoFastjet=True, rParam=0.6 )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)
process.fjSequence = cms.Sequence( process.kt6PFJets25 )

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
process.fjSequence += process.kt6PFJetsForIsolation

process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.pfJetMETcorr.offsetCorrLabel = cms.string("ak5PFL1Fastjet")
process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
process.preProductionSequence = cms.Sequence( process.producePFMETCorrections)

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Trigger matching
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi")
process.load("ElectroWeakAnalysis.MultiBosons.VgTriggerMatcher_cfi")
process.patTriggerEvent.patTriggerMatches  = cms.VInputTag(
    "electronTriggerMatchHLTEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL",
    "electronTriggerMatchHLTEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVL",
    "electronTriggerMatchHLTEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL",
    "electronTriggerMatchHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass50",
    "electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50",
    "electronTriggerMatchHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Mass50",
    "electronTriggerMatchHLTEle27WP80PFMETMT50",
    "electronTriggerMatchHLTEle27WP80",
    "electronTriggerMatchHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "electronTriggerMatchHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoT",
    "electronTriggerMatchHLTEle8Mass50",
    "electronTriggerMatchHLTSC4Mass50",
    "electronTriggerMatchHLTSC17Mass50",
    "photonTriggerMatchHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8Mass50",
    "photonTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50",
    "photonTriggerMatchHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17Mass50",
    "photonTriggerMatchHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "photonTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT",
    "photonTriggerMatchHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoT",
    "photonTriggerMatchHLTEle8Mass50",
    "photonTriggerMatchHLTSC4Mass50",
    "photonTriggerMatchHLTSC17Mass50",
    "photonTriggerMatchHLTPhoton20CaloIdVLIsoL",
    "photonTriggerMatchHLTPhoton26Photon18",
    "photonTriggerMatchHLTPhoton30CaloIdVL",
    "photonTriggerMatchHLTPhoton30CaloIdVLIsoL",
    "photonTriggerMatchHLTPhoton36Photon22",
    "photonTriggerMatchHLTPhoton36CaloId10Iso50Photon22CaloId10Iso50",
    "photonTriggerMatchHLTPhoton50CaloIdVL",
    "photonTriggerMatchHLTPhoton50CaloIdVLIsoL",
    "photonTriggerMatchHLTPhoton75CaloIdVL",
    "photonTriggerMatchHLTPhoton75CaloIdVLIsoL",
    "photonTriggerMatchHLTPhoton90CaloIdVL",
    "photonTriggerMatchHLTPhoton90CaloIdVLIsoL",
    "photonTriggerMatchHLTPhoton135",
    "photonTriggerMatchHLTPhoton150",
    "muonTriggerMatchHLTIsoMu24",
    "muonTriggerMatchHLTIsoMu30",
    "muonTriggerMatchHLTIsoMu24eta2p1",
    "muonTriggerMatchHLTIsoMu30eta2p1",
    "muonTriggerMatchHLTIsoMu34eta2p1",
    "muonTriggerMatchHLTMu17Mu8",
    "muonTriggerMatchHLTMu17TkMu8",
    "muonTriggerMatchHLTMu17forMu17Mu8",
    "muonTriggerMatchHLTMu17forMu17TkMu8",
    "muonTriggerMatchHLTMu8forMu17Mu8",
    "muonTriggerMatchHLTTkMu8forMu17TkMu8",
    "jetTriggerMatchHLTPFJet40",
    "jetTriggerMatchHLTPFJet80",
    "jetTriggerMatchHLTPFJet140",
    "jetTriggerMatchHLTPFJet200",
    "jetTriggerMatchHLTPFJet260",
    "jetTriggerMatchHLTPFJet320",
    "jetTriggerMatchHLTPFJet400"
    )

process.patTrigger.processName = 'HLT'
process.patTriggerEvent.processName = 'HLT'

process.vgTriggerSequence = cms.Sequence(
     process.patTrigger * process.vgTriggerMatcher * process.patTriggerEvent
    )

#================================================
# ADD OFFICIAL ELECTRON ID FROM EG group!
#================================================
# Electron ID with CIC in 2011 
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
# Electron ID with likelihood in 2011 
process.load("RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi")
# Simple cut based selection in 2010
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence +
                                      process.eidVeryLoose +
                                      process.eidLoose +
                                      process.eidMedium +
                                      process.eidTight +
                                      process.eidSuperTight +
                                      process.eidHyperTight1 +
                                      process.eidHyperTight2 +
                                      process.eidHyperTight3 +
                                      process.eidHyperTight4 +
                                      process.eidVeryLooseMC +
                                      process.eidLooseMC +
                                      process.eidMediumMC +
                                      process.eidTightMC +
                                      process.eidSuperTightMC +
                                      process.eidHyperTight1MC +
                                      process.eidHyperTight2MC +
                                      process.eidHyperTight3MC +
                                      process.eidHyperTight4MC +
                                      process.eidLikelihoodExt)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectronIsolation*process.electronMatch*process.patElectrons)
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso= cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso= cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso= cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso= cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso= cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso= cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso= cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso= cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso= cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso= cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso= cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso= cms.InputTag("simpleEleId60cIso"),
    # For CIC
    eidVeryLoose = cms.InputTag("eidVeryLoose"),
    eidLoose = cms.InputTag("eidLoose"),
    eidMedium = cms.InputTag("eidMedium"),
    eidTight = cms.InputTag("eidTight"),
    eidSuperTight = cms.InputTag("eidSuperTight"),
    eidHyperTight1 = cms.InputTag("eidHyperTight1"),
    eidHyperTight2 = cms.InputTag("eidHyperTight2"),
    eidHyperTight3 = cms.InputTag("eidHyperTight3"),
    eidHyperTight4 = cms.InputTag("eidHyperTight4"),
    # For CIC with MC
    eidVeryLooseMC = cms.InputTag("eidVeryLooseMC"),
    eidLooseMC = cms.InputTag("eidLooseMC"),
    eidMediumMC = cms.InputTag("eidMediumMC"),
    eidTightMC = cms.InputTag("eidTightMC"),
    eidSuperTightMC = cms.InputTag("eidSuperTightMC"),
    eidHyperTight1MC = cms.InputTag("eidHyperTight1MC"),
    eidHyperTight2MC = cms.InputTag("eidHyperTight2MC"),
    eidHyperTight3MC = cms.InputTag("eidHyperTight3MC"),
    eidHyperTight4MC = cms.InputTag("eidHyperTight4MC"),
    # likelihood
    eidLikelihoodExt = cms.InputTag("eidLikelihoodExt")
    )

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

process.patJetCorrFactors.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
process.patJetCorrFactors.rho = cms.InputTag('kt6PFJets','rho')
process.patJetCorrFactors.useRho = True

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 outputModules = []
                 )

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')
process.patMETsPF.metSource = cms.InputTag("pfMet")

process.patPFMETsTypeIcorrected = process.patMETs.clone(
    metSource = cms.InputTag('pfType1CorrectedMet'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
)
process.metAnalysisSequence = cms.Sequence( process.patPFMETsTypeIcorrected )

# apply type I/type I + II PFMEt corrections to pat::MET object 
# and estimate systematic uncertainties on MET
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
runMEtUncertainties(process, doSmearJets = True)

process.smearedPFType1CorrectedMet = process.patMETs.clone(
    metSource = cms.InputTag('patType1CorrectedPFMet'),
    addMuonCorrections = cms.bool(False),
    genMETSource = cms.InputTag('genMetTrue'),
    addGenMET = cms.bool(False)
)
process.metAnalysisSequence += process.smearedPFType1CorrectedMet

process.cleanPatPhotons.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.load("ElectroWeakAnalysis.MultiBosons.VgAnalyzerKit_cfi")
process.VgAnalyzerKit.doStoreJets = cms.bool(True)
process.VgAnalyzerKit.doJetHLTMatch = cms.bool(True)
process.VgAnalyzerKit.doStoreSCs = cms.bool(False)
process.VgAnalyzerKit.triggerResults = cms.InputTag("TriggerResults::HLT")
process.VgAnalyzerKit.rhoLabel = cms.InputTag("kt6PFJets25", "rho")
process.VgAnalyzerKit.sigmaLabel = cms.InputTag("kt6PFJets25", "sigma")
process.VgAnalyzerKit.rho2011Label = cms.InputTag("kt6PFJetsForIsolation", "rho")
process.VgAnalyzerKit.rho2012Label = cms.InputTag("kt6PFJets", "rho")
process.VgAnalyzerKit.TypeIpfMETLabel = cms.InputTag("patPFMETsTypeIcorrected")
process.VgAnalyzerKit.SmearedTypeIpfMETLabel = cms.InputTag("smearedPFType1CorrectedMet")

process.VgAnalyzerKit.doSkim = cms.bool(False)

process.TFileService = cms.Service("TFileService", fileName = cms.string('vgtree.root'))

process.p = cms.Path(process.totalKinematicsFilter*
		     process.pfParticleSelectionSequence*
		     process.eleIsoSequence*
		     process.phoIsoSequence*
		     process.fjSequence*
		     process.preProductionSequence*
                     process.patDefaultSequence*
                     process.metAnalysisSequence*
                     process.vgTriggerSequence*
                     process.VgAnalyzerKit)

