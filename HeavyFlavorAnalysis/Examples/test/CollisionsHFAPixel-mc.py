import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "MC_31X_V5::All"


# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource", 
    fileNames = cms.untracked.vstring(
    "file:/shome/starodumov/out/reco/root/reco-10094.root"
    )
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# ----------------------------------------------------------------------
# -- PixelTree
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.src = 'generalTracks'

try:
    rootFileName = os.environ["JOB"] + "-pixel.root"
except KeyError:
    rootFileName = "collisions-pixel.root"

process.PixelTree = cms.EDAnalyzer(
    "PixelTree",
    verbose                = cms.untracked.int32(0),
    rootFileName           = cms.untracked.string(rootFileName),
    muonCollectionLabel    = cms.untracked.string('muons'),
    trackCollectionLabel   = cms.untracked.string('generalTracks'),
    trajectoryInputLabel   = cms.untracked.string('TrackRefitter'),
    pixelClusterLabel      = cms.untracked.string('siPixelClusters'),
    L1GTReadoutRecordLabel = cms.untracked.string("gtDigis"), 
    hltL1GtObjectMap       = cms.untracked.InputTag("hltL1GtObjectMap"), 
    HLTResultsLabel        = cms.untracked.InputTag("TriggerResults::HLT")
    )


# ----------------------------------------------------------------------
process.HepPDTESSource = cms.ESSource(
    "HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/particle.tbl')
    )
process.genParticles = cms.EDProducer(
    "GenParticleProducer",
    saveBarCodes = cms.untracked.bool(True),
    src = cms.InputTag("generator"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
process.genDump = cms.EDFilter(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )


# ----------------------------------------------------------------------
try:
    rootFileName = os.environ["JOB"] + "-hfa.root"
except KeyError:
    rootFileName = "collisions-hfa.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose  = cms.untracked.int32(0),
    requireCand =  cms.untracked.bool(True),
    fileName = cms.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.trkDump = cms.EDFilter(
    "HFDumpTracks",
    verbose = cms.untracked.int32(0),
    generatorEventLabel = cms.untracked.string('generator'),
    muonsLabel = cms.untracked.InputTag("muons"),
    trackingParticlesLabel = cms.untracked.string('trackingParticles'),
    associatorLabel = cms.untracked.string('TrackAssociatorByHits'),
    doTruthMatching = cms.untracked.int32(3),
    tracksLabel = cms.untracked.string('generalTracks'),
    simTracksLabel = cms.untracked.string('allLayer1TrackCands')
    )

# ----------------------------------------------------------------------
process.muonDump = cms.EDFilter(
    "HFDumpMuons",
    verbose = cms.untracked.int32(0),
    muonsLabel = cms.untracked.InputTag("muons"),
    doTruthMatching = cms.untracked.int32(0),
    )


# ----------------------------------------------------------------------
process.triggerDump = cms.EDFilter(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    L1GTReadoutRecordLabel  = cms.untracked.string("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    hltLabel      = cms.untracked.InputTag("TriggerResults::HLT"), 
    )

# ----------------------------------------------------------------------
process.bmtDump = cms.EDFilter(
    "HFMuonAndTrack",
    verbose = cms.untracked.int32(0), 
    muonsLabel = cms.untracked.InputTag("muons"),
    tracksLabel = cms.untracked.string('generalTracks'),
    muonPt = cms.untracked.double(3.0),
    trackPt = cms.untracked.double(1.0),
    type = cms.untracked.int32(1300), 
    massLow  = cms.untracked.double(2.5), 
    massHigh = cms.untracked.double(11.0)
    )

# ----------------------------------------------------------------------
process.bmmDump = cms.EDFilter(
    "HFDimuons",
    verbose = cms.untracked.int32(0), 
    muonsLabel = cms.untracked.InputTag("muons"),
    tracksLabel = cms.untracked.string('generalTracks'),
    muonPt = cms.untracked.double(3.0),
    type = cms.untracked.int32(531), 
    massLow  = cms.untracked.double(4.0), 
    massHigh = cms.untracked.double(6.0)
    )

# ----------------------------------------------------------------------
process.bupsikpDump = cms.EDFilter(
    "HFBu2JpsiKp",
    verbose = cms.untracked.int32(0), 
    muonsLabel = cms.untracked.InputTag("muons"),
    tracksLabel = cms.untracked.string('generalTracks'),
    muonPt = cms.untracked.double(3.0),
    trackPt = cms.untracked.double(0.1),
    type = cms.untracked.int32(100521) 
    )

# ----------------------------------------------------------------------
process.bdpsikstarDump = cms.EDFilter(
    "HFBd2JpsiKstar",
    verbose = cms.untracked.int32(0), 
    muonsLabel = cms.untracked.InputTag("muons"),
    tracksLabel = cms.untracked.string('generalTracks'),
    muonPt = cms.untracked.double(3.0),
    kaonPt = cms.untracked.double(0.1),
    pionPt = cms.untracked.double(0.1),
    type = cms.untracked.int32(100511) 
    )

# ----------------------------------------------------------------------
process.bspsiphiDump = cms.EDFilter(
    "HFBs2JpsiPhi",
    verbose = cms.untracked.int32(0), 
    muonsLabel = cms.untracked.InputTag("muons"),
    tracksLabel = cms.untracked.string('generalTracks'),
    muonPt = cms.untracked.double(3.0),
    trackPt = cms.untracked.double(0.1),
    type = cms.untracked.int32(100531) 
    )


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.TrackRefitter*
    process.PixelTree*
    process.genParticles* 
    process.genDump*
    process.trkDump*
    process.muonDump*
    process.triggerDump*
    process.bmmDump*
    process.bmtDump*
    process.bupsikpDump*
    process.bdpsikstarDump*
    process.bspsiphiDump*
    process.tree
)




