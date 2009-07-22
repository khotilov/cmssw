import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.tauCountFilter_cfi import *
layer1Taus = cms.Sequence(allLayer1Taus * selectedLayer1Taus * countLayer1Taus)
