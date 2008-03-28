import FWCore.ParameterSet.Config as cms

# KFTrajectoryFitterESProducer
from TrackingTools.TrackFitters.KFTrajectoryFitterESProducer_cfi import *
# KFTrajectorySmootherESProducer
from TrackingTools.TrackFitters.KFTrajectorySmootherESProducer_cfi import *
# KFFittingSmootherESProducer
from TrackingTools.TrackFitters.KFFittingSmootherESProducer_cfi import *
# TrackerTrajectoryBuilders
#include "RecoTracker/CkfPattern/data/CkfTrajectoryBuilderESProducerP5.cff"
from RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilderESProducerP5_cff import *
import copy
from RecoTracker.CkfPattern.CkfTrackCandidates_cfi import *
# generate CTF track candidates ############
ckfTrackCandidatesP5 = copy.deepcopy(ckfTrackCandidates)
ckfTrackCandidatesP5.NavigationSchool = 'CosmicNavigationSchool'
ckfTrackCandidatesP5.TrajectoryBuilder = 'GroupedCkfTrajectoryBuilderP5'
#replace ckfTrackCandidatesP5.TrajectoryBuilder        = "CkfTrajectoryBuilderP5"
ckfTrackCandidatesP5.SeedProducer = 'combinatorialcosmicseedfinderP5'

