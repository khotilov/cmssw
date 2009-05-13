# CSA07 Skim for heavy H+->taunu Events
# Filter events with HLT tau + three other jets Et > 20, abs(eta) < 2.5
# then produces AODSIM selected events
# Created by S.Lehti
# Tested on 2007/08/xx/
# Changed to python 2009/05/12

import FWCore.ParameterSet.Config as cms
process = cms.Process("HChSkim")

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)


process.load("HiggsAnalysis.Skimming.heavyChHiggsToTauNu_SkimPaths_cff")
process.load("HiggsAnalysis.Skimming.heavyChHiggsToTauNu_OutputModule_cff")

process.load("FWCore/MessageService/MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	"file:heavyChHiggsToTauNuSkim.root"
    )
)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/HiggsAnalysis/Skimming/test/Attic/heavyChHiggsToTauNu_Filter_cfg.py,v $'),
    annotation = cms.untracked.string('Skim for heavy H+->tau nu events')
)

process.outpath = cms.EndPath(process.heavyChHiggsToTauNuOutputModuleRECOSIM)


