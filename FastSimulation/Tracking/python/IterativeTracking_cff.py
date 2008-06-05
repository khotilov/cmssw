import FWCore.ParameterSet.Config as cms

from FastSimulation.Tracking.PixelTracksProducer_cff import *
from FastSimulation.Tracking.PixelVerticesProducer_cff import *
from FastSimulation.Tracking.IterativeFirstTracking_cff import *
from FastSimulation.Tracking.IterativeSecondTracking_cff import *
from FastSimulation.Tracking.IterativeThirdTracking_cff import *
from FastSimulation.Tracking.GeneralTracks_cfi import *
#include "RecoTracker/FinalTrackSelectors/data/MergeTrackCollections.cff"
iterativeTracking = cms.Sequence(pixelTracking+pixelVertexing+iterativeFirstTracking+iterativeSecondTracking+iterativeThirdTracking+generalTracks)

