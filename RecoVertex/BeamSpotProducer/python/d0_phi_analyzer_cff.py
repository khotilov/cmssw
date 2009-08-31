import FWCore.ParameterSet.Config as cms

d0_phi_analyzer = cms.EDAnalyzer("BeamSpotAnalyzer",
    BSAnalyzerParameters = cms.PSet(
        RunAllFitters = cms.bool(False), ## False: run only default fitter
        WriteToDB = cms.bool(False), ## do not write results to DB
    ),
    BeamFitter = cms.PSet(
	Debug = cms.untracked.bool(False),
	TrackCollection = cms.untracked.InputTag('generalTracks'),
        WriteAscii = cms.untracked.bool(True),
	AsciiFileName = cms.untracked.string('BeamFit.txt'),
	MinimumPt = cms.untracked.double(1.2),
	MaximumEta = cms.untracked.double(2.4),
	MaximumImpactParameter = cms.untracked.double(5),
	MaximumZ = cms.untracked.double(60),
	MinimumTotalLayers = cms.untracked.int32(11),
	MinimumPixelLayers = cms.untracked.int32(3),
	MaximumNormChi2 = cms.untracked.double(2.0),
	TrackAlgorithm = cms.untracked.vstring('iter0'), ## ctf,rs,cosmics,iter0,iter1...; for all algos, leave it blank
	TrackQuality = cms.untracked.vstring(), ## loose, tight, highPurity...; for all qualities, leave it blank
        InputBeamWidth = cms.untracked.double(-1.0), ## if -1 use the value calculated by the analyzer
	FractionOfFittedTrks = cms.untracked.double(0.9),
	MinimumInputTracks = cms.untracked.int32(100)
     ),
    OutputFileName = cms.untracked.string('analyze_d0_phi.root')
)

