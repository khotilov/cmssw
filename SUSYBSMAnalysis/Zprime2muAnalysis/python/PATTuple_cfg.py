#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

# Standard CMSSW configuration (mostly standard in PAT tuple/tools
# use).

process = cms.Process('PAT')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:PlaceHolder.root'))

# Load services needed to run the PAT.
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('PlaceHolder::All')

# Configure the MessageLogger ~sanely. Also direct it to let the PAT
# summary tables be reported -- nice to see how many events had no
# muons, how many had no "selected"/"clean" muons, etc.
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(limit = cms.untracked.int32(-1))

# Define the output file with the output commands defining the
# branches we want to have in our PAT tuple.
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent, patTriggerEventContent
from SUSYBSMAnalysis.Zprime2muAnalysis.EventContent_cff import Zprime2muEventContent
ourEventContent = patEventContent + patTriggerEventContent + Zprime2muEventContent 
process.out = cms.OutputModule('PoolOutputModule',
                               fileName = cms.untracked.string('pat.root'),
                               # If your path in your top-level config is not called 'p', you'll need
                               # to change the next line. In test/PATTuple.py, 'p' is used.
                               SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
                               outputCommands = cms.untracked.vstring('drop *', *ourEventContent)
                               )
process.outpath = cms.EndPath(process.out)

# Load the PAT modules and sequences, and configure them as we
# need. See the individual functions for their documentation.  MC use
# is assumed by default, and should be removed after everything's
# configured in the top-level config using removeMCUse() tool if
# running on data.  (This is due to the design of the PAT: easier to
# do it in this order rather than adding things for MC use later.)
process.load('PhysicsTools.PatAlgos.patSequences_cff')

from PATTools import pruneMCLeptons, addMuonMCClassification, addMuonStations, addMuonHitCount, addHEEPId
pruneMCLeptons(process, use_sim=True) # need to decide whether to move AODOnly() call in here, if so use_sim should just be set False
addMuonMCClassification(process)
addMuonStations(process)
addMuonHitCount(process)
addHEEPId(process)

from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger, switchOnTriggerMatchEmbedding
switchOnTrigger(process)
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.hltTriggerMatch_cfi')
switchOnTriggerMatchEmbedding(process, triggerMatchers=['muonTriggerMatchHLTMuons'])
process.out.outputCommands += ['keep *_cleanPatMuonsTriggerMatch_*_*', 'drop *_cleanPatMuons_*_*']

# Some extra configuration of the PAT.

# Embed the tracker tracks (by default, every other track is already
# embedded).
process.patMuons.embedTrack = True

# Follow VBTF for now in using the beamspot for "correcting" dxy,
# instead of the primary vertex.
process.patMuons.usePV = False

# Define simple quality cuts for muons (analysis cuts to be done at the analysis level).
process.selectedPatMuons.cut = 'isGlobalMuon || isTrackerMuon'

# Filter out events with no selected muons. (countPatMuons counts
# those muons in the cleanPatMuons collection, which by default is
# just a copy of the selectedPatMuons collection.) If, e.g., wanting
# to use the PAT tuple for efficiency studies, or running on just
# electrons, need to change appropriately in your top-level config
# (perhaps using countPatLeptons instead).
process.countPatMuons.minNumber = 1

# Instead of filtering out events at PAT-tupling time based on things
# like GoodVertex and NoScraping, schedule separate paths for all the
# "good data" filters so that the results of them get stored in a
# small TriggerResults::PAT object. This can be read and used to
# filter events in the analyzer process using e.g. the filter in
# Zprime2muAnalysis_cff.py. (This is useful so we don't have to keep
# around the entire generalTracks collection to run the NoScraping
# filter later, for example.)
#
# Make one path for each (a very small storage burden) so they can be
# accessed separately in the TriggerResults object; the "All" path
# isn't necessary because it could be emulated using the AND of all of
# the separate ones, but it's nice for convenience.
process.load('SUSYBSMAnalysis.Zprime2muAnalysis.goodData_cff')
process.goodDataHLTPhysicsDeclared = cms.Path(process.hltPhysicsDeclared)
process.goodDataPrimaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.goodDataNoScraping = cms.Path(process.noscraping)
process.goodDataAll = cms.Path(process.hltPhysicsDeclared * process.primaryVertexFilter * process.noscraping)
