import FWCore.ParameterSet.Config as cms

process = cms.Process("EWKDiElectronSkim")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
    'rfio:/tmp/ikesisog/B0E3F649-1B80-DE11-96AF-001D0967E048.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')


# HLT filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.EWK_DiElectronHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

# Uncomment this to access 8E29 menu and filter on it
process.EWK_DiElectronHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT8E29")
process.EWK_DiElectronHLTFilter.HLTPaths = ["HLT_Ele15_LW_L1R"]


#   Make a collection of good SuperClusters.
#
#   Before selection is made, merge the Barrel and EndCap SC's.
process.superClusterMerger =  cms.EDFilter("EgammaSuperClusterMerger",
                                    src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'), cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                  )

#   Get the above merged SC's and select the particle (gamma) to greate SC's Candidates.
process.superClusterCands = cms.EDProducer("ConcreteEcalCandidateProducer",
                                    src = cms.InputTag("superClusterMerger"), particleType = cms.string('gamma')
                                  )

#   Get the above SC's Candidates and place a cut on their Et.
process.goodSuperClusters = cms.EDFilter("CandViewRefSelector",
                                    src = cms.InputTag("superClusterCands"),
                                    cut = cms.string('et > 20.0'),
                                    filter = cms.bool(True)
                                )

process.superClusterFilter = cms.Sequence(process.superClusterMerger + process.superClusterCands + process.goodSuperClusters)


#   Make a collections on good Electrons.
#
process.goodElectrons = cms.EDFilter("CandViewSelector",
                                src = cms.InputTag("gsfElectrons"),
                                cut = cms.string('pt > 20.0'),
                                filter = cms.bool(True)                                
)

#   Filter the above two collections (good SuperClusters and good Electrons)
#
process.electronSuperClusterCombiner = cms.EDFilter("CandViewShallowCloneCombiner",
                                            filter = cms.bool(True),
                                            checkCharge = cms.bool(False),
                                            cut = cms.string('mass > 3.0'),
                                            decay = cms.string('goodElectrons goodSuperClusters')
                                            )

process.electronSuperClusterCounter = cms.EDFilter("CandViewCountFilter",
                                            src = cms.InputTag("electronSuperClusterCombiner"),
                                            minNumber = cms.uint32(1)
                                          )

process.electronSuperClusterFilter = cms.Sequence(process.electronSuperClusterCombiner + process.electronSuperClusterCounter)


# Skim path
process.EWK_DiElectronSkimPath = cms.Path(process.EWK_DiElectronHLTFilter +
                                          process.goodElectrons + 
                                          process.superClusterFilter + 
                                          process.electronSuperClusterFilter
                                         )


# Output module configuration
from Configuration.EventContent.EventContent_cff import *

EWK_DiElectronSkimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring()
)
EWK_DiElectronSkimEventContent.outputCommands.extend(AODEventContent.outputCommands)

EWK_DiElectronSkimEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
           'EWK_DiElectronSkimPath')
    )
)

process.EWK_DiElectronSkimOutputModule = cms.OutputModule("PoolOutputModule",
                                            EWK_DiElectronSkimEventContent,
                                            EWK_DiElectronSkimEventSelection,
                                            dataset = cms.untracked.PSet(
                                                filterName = cms.untracked.string('EWKDiElectronSkim'),
                                                dataTier = cms.untracked.string('USER')
                                            ),
                                            fileName = cms.untracked.string('testEWKDiElectronSkim.root')
)

process.outpath = cms.EndPath(process.EWK_DiElectronSkimOutputModule)


