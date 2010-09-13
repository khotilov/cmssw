import FWCore.ParameterSet.Config as cms
#import os

process = cms.Process("REPROD")

# General
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
# Global tag for 336patch3
#process.GlobalTag.globaltag = 'GR09_R_V5::All'
#process.GlobalTag.globaltag = 'GR09_R_V6::All'
# Global tag for 341
process.GlobalTag.globaltag = 'GR_R_38X_V9::All'


# Add PF vertices from Maxime
#process.load("RecoParticleFlow.PFTracking.particleFlowDisplacedVertexCandidate_cff")
#process.load("RecoParticleFlow.PFTracking.particleFlowDisplacedVertex_cff")
#process.particleFlowDisplacedVertexCandidate.primaryVertexCut = cms.double(2.0)
#process.particleFlowDisplacedVertex.primaryVertexCut = cms.double(2)
#process.particleFlowDisplacedVertex.tobCut = cms.double(100)
#process.particleFlowDisplacedVertex.tecCut = cms.double(200)

# Other statements

#####################################################################################################
####
####  Top level replaces for handling strange scenarios of early collisions
####

## TRACKING:
## Skip events with HV off
process.newSeedFromTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.newSeedFromPairs.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
process.secTriplets.ClusterCheckPSet.MaxNumberOfPixelClusters=2000
process.fifthSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters = 20000
process.fourthPLSeeds.ClusterCheckPSet.MaxNumberOfCosmicClusters=20000
process.thTripletsA.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000
process.thTripletsB.ClusterCheckPSet.MaxNumberOfPixelClusters = 5000

###### FIXES TRIPLETS FOR LARGE BS DISPLACEMENT ######

### prevent bias in pixel vertex
process.pixelVertices.useBeamConstraint = False

### pixelTracks
#---- new parameters ----
process.pixelTracks.RegionFactoryPSet.RegionPSet.nSigmaZ  = 4.06
process.pixelTracks.RegionFactoryPSet.RegionPSet.originHalfLength = cms.double(40.6)

### 0th step of iterative tracking
#---- new parameters ----
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ   = cms.double(4.06)  
process.newSeedFromTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 40.6

### 2nd step of iterative tracking
#---- new parameters ----
process.secTriplets.RegionFactoryPSet.RegionPSet.nSigmaZ  = cms.double(4.47)  
process.secTriplets.RegionFactoryPSet.RegionPSet.originHalfLength = 44.7

## Primary Vertex
process.offlinePrimaryVerticesWithBS.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVerticesWithBS.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVerticesWithBS.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVerticesWithBS.TkClusParameters.TkGapClusParameters.zSeparation = 1
process.offlinePrimaryVertices.PVSelParameters.maxDistanceToBeam = 2
process.offlinePrimaryVertices.TkFilterParameters.maxNormalizedChi2 = 20
process.offlinePrimaryVertices.TkFilterParameters.maxD0Significance = 100
process.offlinePrimaryVertices.TkFilterParameters.minPixelLayersWithHits = 2
process.offlinePrimaryVertices.TkFilterParameters.minSiliconLayersWithHits = 5
process.offlinePrimaryVertices.TkClusParameters.TkGapClusParameters.zSeparation = 1

## ECAL 
process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]


## HCAL temporary fixes
process.hfreco.samplesToAdd = 4
    
## EGAMMA
process.photons.minSCEtBarrel = 5.
process.photons.minSCEtEndcap =5.
process.photonCore.minSCEt = 5.
process.conversionTrackCandidates.minSCEt =5.
process.conversions.minSCEt =5.
process.trackerOnlyConversions.rCut = 2.
process.trackerOnlyConversions.vtxChi2 = 0.0005

process.hfreco.firstSample=3

## local tracker strip reconstruction
#process.OutOfTime.TOBlateBP=0.071
#process.OutOfTime.TIBlateBP=0.036

## particle flow HF cleaning
process.particleFlowRecHitHCAL.LongShortFibre_Cut = 30.
process.particleFlowRecHitHCAL.ApplyTimeDPG = False
process.particleFlowRecHitHCAL.ApplyPulseDPG = True
process.particleFlowRecHitECAL.timing_Cleaning = True

## HF cleaning for data only
process.hcalRecAlgos.SeverityLevels[3].RecHitFlags.remove("HFDigiTime")
process.hcalRecAlgos.SeverityLevels[4].RecHitFlags.append("HFDigiTime")
    
###
###  end of top level replacements
###
###############################################################################################

# All events with PFMET > 30 GeV in JSON'ed runs (before or after cleaning)
#process.load("PFAnalyses.PFCandidate.METSkim30_ReReco370_JSON_cff")
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
      #'file:Holger_EM_1.root',
      #'file:Holger_EM_2.root',
      #'file:Holger_JM_1.root',
      #'file:Holger_MM_1.root',
      'file:Holger_GM_1.root',
      'file:Holger_GM_2.root',
      'file:Holger_GM_3.root',
      'file:Holger_GM_4.root'
      #'file:roecker_1.root',
      #'file:roecker_2.root',
      #'file:roecker_3.root',
      )
    )
process.source.secondaryFileNames = cms.untracked.vstring()
process.source.noEventSort = cms.untracked.bool(True)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# This is for filtering on L1 technical trigger bit: MB and no beam halo
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(0 AND (36 OR 37 OR 38 OR 39))')

process.scrapping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
#process.HBHENoiseFilter.maxRBXEMF = cms.double(0.01)

#process.tkHVON = cms.EDFilter("PhysDecl",
#                              applyFilter=cms.untracked.bool(True)
#                              )


process.dump = cms.EDAnalyzer("EventContentAnalyzer")


process.load("RecoParticleFlow.Configuration.ReDisplay_EventContent_cff")
process.display = cms.OutputModule("PoolOutputModule",
    process.DisplayEventContent,
    fileName = cms.untracked.string('display_Holger_GM_2.root'),
    #fileName = cms.untracked.string('display_roecker.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
)

process.load("Configuration.EventContent.EventContent_cff")
process.rereco = cms.OutputModule("PoolOutputModule",
    process.RECOSIMEventContent,
    #fileName = cms.untracked.string('NoFilter_METSkimPFClean30.root')
    fileName = cms.untracked.string('reco.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
)

# Maxime !!!@$#^%$^%#@
#process.particleFlowDisplacedVertexCandidate.verbose = False
#process.particleFlowDisplacedVertex.verbose = False

# Local re-reco: Produce tracker rechits, pf rechits and pf clusters
process.towerMakerPF.HcalAcceptSeverityLevel = 9
process.localReReco = cms.Sequence(process.siPixelRecHits+
                                   process.siStripMatchedRecHits+
                                   process.particleFlowCluster)

#Photon re-reco
process.photonReReco = cms.Sequence(process.conversionSequence+
                                    process.trackerOnlyConversionSequence+
                                    process.photonSequence+
                                    process.photonIDSequence)

# Track re-reco
process.globalReReco =  cms.Sequence(process.offlineBeamSpot+
                                     process.recopixelvertexing+
                                     process.ckftracks+
                                     process.ctfTracksPixelLess+
                                     process.offlinePrimaryVertices *
                                     process.offlinePrimaryVerticesWithBS *
                                     process.caloTowersRec+
                                     process.vertexreco+
                                     process.recoJets+
                                     process.muonrecoComplete+
                                     process.electronGsfTracking+
                                     process.photonReReco+
                                     process.metreco)



# Particle Flow re-processing
process.pfReReco = cms.Sequence(process.particleFlowReco+
                                process.recoPFJets+
                                process.recoPFMET+
                                process.PFTau#+
#                                process.particleFlowDisplacedVertexCandidate+
#                                process.particleFlowDisplacedVertex
                                )
                                
# Gen Info re-processing
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.load("RecoMET.Configuration.RecoGenMET_cff")
process.load("RecoParticleFlow.PFProducer.particleFlowSimParticle_cff")
process.load("RecoParticleFlow.Configuration.HepMCCopy_cfi")
process.genReReco = cms.Sequence(process.generator+
                                 process.genParticles+
                                 process.genJetParticles+
                                 process.recoGenJets+
                                 process.genMETParticles+
                                 process.recoGenMET+
                                 process.particleFlowSimParticle)

# The complete reprocessing
process.p = cms.Path(#process.hltLevel1GTSeed+
                     #process.bxSelect+
                     process.scrapping+
                     process.HBHENoiseFilter+
                     #process.tkHVON+
                     process.localReReco+
                     process.globalReReco+
                     process.pfReReco#+
                     #process.genReReco
                     )

# And the output.
# Write out only filtered events
process.display.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
process.rereco.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
#process.outpath = cms.EndPath(process.rereco+process.display)
process.outpath = cms.EndPath(process.display)


# Schedule the paths
process.schedule = cms.Schedule(
    process.p,
    process.outpath
)

# And the logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    #fileMode = cms.untracked.string('NOMERGE'),
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),
    Rethrow = cms.untracked.vstring('Unknown', 
        'ProductNotFound', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ModuleFailure', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError', 
        'FatalRootError', 
        'NotFound')
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1

