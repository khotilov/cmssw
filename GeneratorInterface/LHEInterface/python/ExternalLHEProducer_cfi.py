import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer('ExternalLHEProducer',
    scriptName = cms.FileInPath("GeneratorInterface/LHEInterface/data/run_madgraph_gridpack.sh"),
    outputFile = cms.string("events_final.lhe"),
    args = cms.vstring('slc5_ia32_gcc434/madgraph/V5_1.3.27/test','W1Jet_7TeV_madgraph','100','false','true','wjets','5','20','false','0','99')
)

