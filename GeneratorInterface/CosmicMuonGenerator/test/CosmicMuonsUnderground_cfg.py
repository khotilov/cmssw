import FWCore.ParameterSet.Config as cms

process = cms.Process("runCosMuoGen")
process.load("GeneratorInterface.CosmicMuonGenerator.CMSCGENproducer_cfi")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        engineName = cms.untracked.string('HepJamesRandom')
    )
)

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)
process.CMSCGEN_out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('cosmic.root')
)

process.p = cms.Path(process.generator)
process.outpath = cms.EndPath(process.CMSCGEN_out)
process.generator.MinP = 10.
process.generator.MaxTheta = 80.

# Plug z-position [mm] (default=-14000.)
#process.generator.PlugVz = -33000.;
