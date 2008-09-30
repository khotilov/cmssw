import FWCore.ParameterSet.Config as cms

from RecoTracker.IterativeTracking.FirstFilter_cfi import *
from RecoTracker.IterativeTracking.SecStep_cff import *
from RecoTracker.IterativeTracking.ThStep_cff import *
from RecoTracker.IterativeTracking.PixelLessStep_cff import *
#iterTracking = cms.Sequence(firstfilter*secondStep*thirdStep*fourthStep)
iterTracking = cms.Sequence(firstfilter*secondStep*thirdStep)

