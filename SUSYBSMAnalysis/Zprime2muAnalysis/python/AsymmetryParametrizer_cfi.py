import FWCore.ParameterSet.Config as cms

from AsymFitManager_cff import AsymFitManager

AsymmetryParametrizer = cms.EDAnalyzer('AsymmetryParametrizer',
                                       AsymFitManager,
                                       gen_particle_src = cms.InputTag('genParticles'),
                                       lepton_flavor = cms.int32(13),
                                       internal_brem_on = cms.bool(True),
                                       assemble_only = cms.bool(False),
                                       histos_fn = cms.string('AsymmetryParametrizer.root'),
                                       postscript_fn = cms.string('AsymmetryParametrizer.ps'),
                                       )
