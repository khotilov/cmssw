import FWCore.ParameterSet.Config as cms

from RecoJets.FFTJetProducers.fftjetcommon_cfi import *

fftjet_pf_pileup_cleaner = cms.EDProducer(
    "FFTJetPFPileupCleaner",
    #
    # Label for the input collection of PFCandidate objects
    PFCandidates = cms.InputTag("particleFlow"),
    #
    # Label for the collection of primary vertices
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    #
    # Find the closest vertex even if the track is not associated
    # with any good vertex?
    checkClosestZVertex = cms.bool(True),
    #
    # Remove the objects associated with the main primary vertex?
    removeMainVertex = cms.bool(False),
    #
    # Remove the objects not associated with any primary vertex?
    removeUnassociated = cms.bool(False),
    #
    # Overall flag to invert the decision
    reverseRemovalDecision = cms.bool(False),
    #
    # Various removal flags by object type. See PFCandidate header
    # for dobject type etails.
    remove_X = cms.bool(False),
    remove_h = cms.bool(True),
    remove_e = cms.bool(True),
    remove_mu = cms.bool(True),
    remove_gamma = cms.bool(False),
    remove_h0 = cms.bool(False),
    remove_h_HF = cms.bool(False),
    remove_egamma_HF = cms.bool(False),
    #
    # Minimum and maximum allowed eta
    etaMin = cms.double(-fftjet_standard_eta_range),
    etaMax = cms.double(fftjet_standard_eta_range),
    #
    # Vertex quality cut
    vertexNdofCut = cms.double(4.0)
)
