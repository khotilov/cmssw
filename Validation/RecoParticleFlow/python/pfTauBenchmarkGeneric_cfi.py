import FWCore.ParameterSet.Config as cms

jets = 'iterativeCone5PFJets'

pfTauBenchmarkGeneric = cms.EDAnalyzer("GenericBenchmarkAnalyzer",
    OutputFile = cms.untracked.string('tauBenchmarkGeneric.root'),
    InputTruthLabel = cms.InputTag('tauGenJets'),
    maxEta = cms.double(2.5),
    recPt = cms.double(5.0),
    deltaRMax = cms.double(0.3),
    PlotAgainstRecoQuantities = cms.bool(False),
    BenchmarkLabel = cms.string( jets ),
    InputRecoLabel = cms.InputTag( jets )
)
