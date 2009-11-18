import FWCore.ParameterSet.Config as cms


import RecoParticleFlow.PFAnalyses.pflowCalibratable_cfi as calibratable 
from RecoParticleFlow.PFAnalyses.pflowFaketracks_cfi import *
from RecoParticleFlow.PFAnalyses.pflowProcessTestbeam_cfi import *
from RecoParticleFlow.PFAnalyses.pflowParticleFiltration_cfi import *
#Need to override clustering to exclude HF components

from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
pfClusteringHCALTB = cms.Sequence(particleFlowRecHitHCAL * particleFlowClusterHCAL)

TFileService = cms.Service("TFileService",
    fileName=cms.string("PFlowTB_Tree.root")
)

finishup = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string("PFlowTB_Events.root"),
    outputCommands=cms.untracked.vstring('keep *'),
    #outputCommands=cms.untracked.vstring('drop *', 'keep *_particleFiltration_*_*', 'keep recoMuons_*_*_*', 'keep *_calibratable_*_*', 'keep *_faketracks_*_*', 'keep recoPFRecTracks_*_*_*', 'keep recoPFRecHits_*_*_*', 'keep recoPFClusters_*_*_*', 'keep recoPFBlocks_*_*_*', 'keep recoPFCandidates_*_*_*'),
	SelectEvents=cms.untracked.PSet(
                SelectEvents=cms.vstring('p1')
    ) 
)

cleanDump = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string("CleanTB_Events.root"),
    outputCommands=cms.untracked.vstring('keep *'),
    SelectEvents=cms.untracked.PSet(
                SelectEvents=cms.vstring('p1')
    ) 
)

#Cleaning only
pflowCleaning = cms.Sequence(particleFiltration)

#Clustering only
pflowClusteringTestbeam = cms.Sequence(pfClusteringECAL * pfClusteringHCALTB)

#Tracking, clustering, pflow - no cleaning
pflowNoCleaning = cms.Sequence(faketracks * pfClusteringECAL * pfClusteringHCALTB * particleFlowBlock * particleFlow)

extraction = cms.EDAnalyzer("ExtractionAnalyzer",
    calibratable.TestbeamDelegate
)

pflowCalibEcalRechits = cms.EDProducer("EcalCalibRechitProducer",
    EEUncalColl=cms.InputTag("ecal2007TBH2WeightUncalibRecHit", "EcalUncalibRecHitsEE"),
    EENoisesFile=cms.untracked.string("/afs/cern.ch/user/b/ballin/scratch0/cmssw/src/RecoParticleFlow/PFAnalyses/macros/EE_noises.txt"),
    EECoeffsFile=cms.untracked.string("/afs/cern.ch/user/b/ballin/scratch0/cmssw/src/RecoParticleFlow/PFAnalyses/macros/ee_calib_test.txt"),
)


towerMakerPF.ecalInputs = cms.VInputTag(
    cms.InputTag("pflowCalibEcalRechits","EcalRecHitsEB"), 
    cms.InputTag("pflowCalibEcalRechits","EcalRecHitsEE"))
towerMakerPF.hfInput=cms.InputTag("");
towerMakerPF.hoInput=cms.InputTag("");
towerMakerPF.UseHO = cms.bool(False)
towerMakerPF.AllowMissingInputs=cms.bool(True)

pflowEndcapRechitMaker = cms.Sequence(particleFiltration * faketracks * pflowCalibEcalRechits * towerMakerPF)


pfAlgoAndExtractionTestbeam = cms.Sequence(faketracks * 
                                    pfClusteringECAL * 
                                    pfClusteringHCALTB * 
                                    particleFlowBlock * 
                                    particleFlow * 
                                    extraction)

pfAlgoAndExtractionTestbeamEndcaps = cms.Sequence(faketracks * 
                                    pflowCalibEcalRechits *
                                    towerMakerPF *
                                    pfClusteringECAL * 
									pfClusteringPS *
                                    pfClusteringHCALTB * 
                                    particleFlowBlock * 
                                    particleFlow * 
                                    extraction)

#The works.
pflowProcessTestbeam = cms.Sequence(particleFiltration * 
                                    pfAlgoAndExtractionTestbeam)

#The works for the endcap.
pflowProcessEndcapTestbeam = cms.Sequence(particleFiltration * 
                                   pfAlgoAndExtractionTestbeamEndcaps)

