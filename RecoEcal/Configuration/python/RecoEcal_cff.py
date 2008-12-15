import FWCore.ParameterSet.Config as cms

# Sequence for clustering in ecal barrel & endcap + preshower
# geometry needed for clustering
from RecoEcal.EgammaClusterProducers.geometryForClustering_cff import *
# hybrid clustering sequence
from RecoEcal.EgammaClusterProducers.hybridClusteringSequence_cff import *
# multi5x5 clusters
from RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff import *
# preshower sequence for multi5x5 clusters
from RecoEcal.EgammaClusterProducers.multi5x5PreshowerClusteringSequence_cff import *
#reduced recHit collection producer
from RecoEcal.EgammaClusterProducers.reducedRecHitsSequence_cff import *
# create path with all clustering algos
# NB: preshower MUST be run after multi5x5 clustering in the endcap
ecalClusters = cms.Sequence(hybridClusteringSequence*multi5x5ClusteringSequence*multi5x5PreshowerClusteringSequence)
