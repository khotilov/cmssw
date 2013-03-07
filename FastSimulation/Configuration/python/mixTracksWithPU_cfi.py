import FWCore.ParameterSet.Config as cms

from FastSimulation.Configuration.mixFastSimObjects_cfi import *

mixSimTracksAndVertices = cms.EDProducer("MixingModule",
                                digitizers = cms.PSet(),
                                LabelPlayback = cms.string(''),
                                maxBunch = cms.int32(0),
                                minBunch = cms.int32(0),
                                bunchspace = cms.int32(25),
                                checktof = cms.bool(False),
                                playback = cms.untracked.bool(False),
                                mixProdStep1 = cms.bool(False),
                                mixProdStep2 = cms.bool(False),
                                input = cms.SecSource("PoolSource",
                                                      type = cms.string('probFunction'),
                                                      nbPileupEvents = cms.PSet(
    probFunctionVariable = cms.vint32(),
    probValue = cms.vdouble(),
    histoFileName = cms.untracked.string('histProbFunction.root'),
    ),
                                                      sequential = cms.untracked.bool(False),
                                                      manage_OOT = cms.untracked.bool(False),  ## manage out-of-time pileup
                                                      ## setting this to True means that the out-of-time pileup
                                                      ## will have a different distribution than in-time, given
                                                      ## by what is described on the next line:
                                                      OOT_type = cms.untracked.string('None'),  ## generate OOT with a Poisson matching the number chosen for in-time
                                                      #OOT_type = cms.untracked.string('fixed'),  ## generate OOT with a fixed distribution
                                                      #intFixed_OOT = cms.untracked.int32(2),
                                                      fileNames = cms.untracked.vstring('file:MinBias_test.root'),
#                                                      fileNames = cms.untracked.vstring('file:SingleMuPt10_forPU.root'),
                                                      ),
                                mixObjects = cms.PSet(
    mixSH = cms.PSet(
    mixSimHits
    ),
    mixVertices = cms.PSet(
    mixSimVertices
    ),
    mixSTracks = cms.PSet(
    mixSimTracks
    )
    )
)
