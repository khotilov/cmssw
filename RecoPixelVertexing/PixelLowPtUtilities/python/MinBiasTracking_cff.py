import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelLowPtUtilities.common_cff     import *
from RecoPixelVertexing.PixelLowPtUtilities.firstStep_cff  import *
from RecoPixelVertexing.PixelLowPtUtilities.secondStep_cff import *
from RecoPixelVertexing.PixelLowPtUtilities.thirdStep_cff  import *

from RecoVZero.VZeroFinding.VZeros_cff import *

###################################
# First step, triplets, r=0.2 cm
firstStep  = cms.Sequence(pixel3ProtoTracks
                        * pixel3Vertices
                        * pixel3PrimTracks
                        * primSeeds
                        * primTrackCandidates
                        * globalPrimTracks)

###################################
# Second step, triplets, r=3.5 cm
secondStep = cms.Sequence(secondClusters
                        * secondPixelRecHits
                        * secondStripRecHits
                        * pixelSecoTracks
                        * secoSeeds
                        * secoTrackCandidates
                        * globalSecoTracks)

###################################
# Third step, pairs, r=0.2 cm
thirdStep  = cms.Sequence( thirdClusters
                         * thirdPixelRecHits
                         * thirdStripRecHits
                         * pixelTertTracks
                         * tertSeeds
                         * tertTrackCandidates
                         * globalTertTracks)

###################################
# Tracklist combiner
allTracks = cms.EDFilter("TrackListCombiner",
    trackProducers = cms.vstring('globalPrimTracks',
                                 'globalSecoTracks',
                                 'globalTertTracks')
)

###################################
# Minimum bias tracking
minBiasTracking = cms.Sequence(firstStep
                            * secondStep
                             * thirdStep 
                             * allTracks)

