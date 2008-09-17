import FWCore.ParameterSet.Config as cms

# needed backend
DQMStore = cms.Service("DQMStore",
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0),
    collateHistograms = cms.untracked.bool(False)
)

# actual producer
from DQMServices.Components.MEtoEDMConverter_cfi import *

