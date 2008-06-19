import FWCore.ParameterSet.Config as cms

# $Id: sisCone5JTA_cff.py,v 1.2 2008/04/21 03:27:52 rpw Exp $
from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi import * ##propagator

from RecoJets.JetAssociationProducers.j2tParametersCALO_cfi import *
from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import *
sisCone5JetTracksAssociatorAtVertex = cms.EDFilter("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("sisCone5CaloJets")
)

sisCone5JetTracksAssociatorAtCaloFace = cms.EDFilter("JetTracksAssociatorAtCaloFace",
    j2tParametersCALO,
    jets = cms.InputTag("sisCone5CaloJets")
)

sisCone5JetExtender = cms.EDFilter("JetExtender",
    jets = cms.InputTag("sisCone5CaloJets"),
    jet2TracksAtCALO = cms.InputTag("sisCone5JetTracksAssociatorAtCaloFace"),
    jet2TracksAtVX = cms.InputTag("sisCone5JetTracksAssociatorAtVertex"),
    coneSize = cms.double(0.5)
)

sisCone5JTA = cms.Sequence(sisCone5JetTracksAssociatorAtVertex*sisCone5JetTracksAssociatorAtCaloFace*sisCone5JetExtender)

