import FWCore.ParameterSet.Config as cms

# Iterative steps
from RecoTracker.IterativeTracking.iterativeTk_cff import *


# RS
from RecoTracker.RoadSearchSeedFinder.RoadSearchSeeds_cff import *
from RecoTracker.RoadSearchCloudMaker.RoadSearchClouds_cff import *
from RecoTracker.RoadSearchTrackCandidateMaker.RoadSearchTrackCandidates_cff import *
from RecoTracker.TrackProducer.RSFinalFitWithMaterial_cff import *

### Not the Tracking uses the 2 seed collections separately. The merged seed collection is produced 
### for backward compatibility with electron reconstruction
newCombinedSeeds = RecoTracker.TkSeedGenerator.GlobalCombinedSeeds_cfi.globalCombinedSeeds.clone()
newCombinedSeeds.seedCollections = cms.VInputTag(
    cms.InputTag('newSeedFromTriplets'),
    cms.InputTag('newSeedFromPairs'),
)


ckftracks = cms.Sequence(iterTracking*trackCollectionMerging*newCombinedSeeds)

rstracks = cms.Sequence(roadSearchSeeds*
                        roadSearchClouds*rsTrackCandidates*
                        rsWithMaterialTracks)



