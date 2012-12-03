import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("ReggeGribovPartonMCGeneratorFilter",
                    beammomentum = cms.double(1380),
                    targetmomentum = cms.double(-1380),
                    beamid = cms.int32(208),
                    targetid = cms.int32(208),
                    model = cms.int32(7),
                    bmin = cms.double(0),
                    bmax = cms.double(10000),
                    paramFileName = cms.untracked.string("Configuration/Generator/data/ReggeGribovPartonMC.param")
                    )


configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    name = cms.untracked.string('$Source: /local/reps/CMSSW/CMSSW/GeneratorInterface/ReggeGribovPartonMCInterface/python/ReggeGribovPartonMC_cfi.py,v $'),
    annotation = cms.untracked.string('ReggeGribovMC generator')
    )





