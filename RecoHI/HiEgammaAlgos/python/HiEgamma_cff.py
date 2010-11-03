import FWCore.ParameterSet.Config as cms

# clustering sequence
from RecoHI.HiEgammaAlgos.HiIslandClusteringSequence_cff import *
from RecoEcal.EgammaClusterProducers.hybridClusteringSequence_cff import *
from RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff import *
from RecoEcal.EgammaClusterProducers.multi5x5PreshowerClusteringSequence_cff import *
from RecoEcal.EgammaClusterProducers.preshowerClusteringSequence_cff import *
from RecoEcal.EgammaClusterProducers.dynamicHybridClusteringSequence_cff import *

hiEcalClusteringSequence = cms.Sequence(islandClusteringSequence*hybridClusteringSequence*dynamicHybridClusteringSequence*multi5x5ClusteringSequence*multi5x5PreshowerClusteringSequence*preshowerClusteringSequence)

# high purity tracks
highPurityTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("hiSelectedTracks"),
    cut = cms.string('quality("highPurity")')
)

# reco photon producer
from RecoEgamma.EgammaPhotonProducers.photonSequence_cff import *

photonCore.scHybridBarrelProducer = cms.InputTag("correctedIslandBarrelSuperClusters") # use island for the moment
photonCore.scIslandEndcapProducer = cms.InputTag("correctedIslandEndcapSuperClusters") # use island for the moment
photonCore.minSCEt = cms.double(8.0)
photons.minSCEtBarrel = cms.double(5.0)
photons.minSCEtEndcap = cms.double(15.0)
photons.minR9Barrel = cms.double(10.)  #0.94
photons.minR9Endcap = cms.double(10.)   #0.95
photons.maxHoverEEndcap = cms.double(0.5)  #0.5
photons.maxHoverEBarrel = cms.double(0.99)  #0.5
photons.primaryVertexProducer = cms.string('hiSelectedVertex') # replace the primary vertex
photons.isolationSumsCalculatorSet.trackProducer = cms.InputTag("highPurityTracks")

hiPhotonSequence = cms.Sequence(highPurityTracks*photonSequence)

# HI Egamma Isolation
from RecoHI.HiEgammaAlgos.HiEgammaIsolation_cff import *

# HI Ecal reconstruction
hiEcalClusters = cms.Sequence(hiEcalClusteringSequence)
hiEgammaSequence = cms.Sequence(hiPhotonSequence)
hiEcalClustersIsolation = cms.Sequence(hiEgammaSequence * hiEgammaIsolationSequence)

# HI Spike Clean Sequence
hiSpikeCleanedSC = cms.EDProducer("HiSpikeCleaner",
                                  recHitProducerBarrel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
                                  recHitProducerEndcap = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
                                  originalSuperClusterProducer = cms.InputTag("correctedIslandBarrelSuperClusters"), 
                                  outputColl  = cms.string( "" ),
                                  etCut          = cms.double(10),
                                  TimingCut    = cms.untracked.double(4.0),
                                  swissCutThr    = cms.untracked.double(0.83)
                                  )

cleanPhotonCore = photonCore.clone()
cleanPhotons = photons.clone()
cleanPhotonCore.scHybridBarrelProducer = cms.InputTag("hiSpikeCleanedSC")
cleanPhotons.photonCoreProducer = cms.InputTag("cleanPhotonCore")

hiPhotonCleaningSequence = cms.Sequence(hiSpikeCleanedSC * cleanPhotonCore * cleanPhotons)
