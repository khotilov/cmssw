import FWCore.ParameterSet.Config as cms

process = cms.Process("Print")

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring('file:COMPHEP_ZP331mumu_M1000_Mcut400.root'),
	skipEvents = cms.untracked.uint32(0)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDFilter("ParticleListDrawer",
	maxEventsToPrint = cms.untracked.int32(1)
)

process.printTree = cms.EDFilter("ParticleTreeDrawer",
	status = cms.untracked.vint32(3),
	src = cms.InputTag("genParticles")##,
## 	printP4 = cms.untracked.bool(False),
## 	printStatus = cms.untracked.bool(True),
## 	printIndex = cms.untracked.bool(True),
## 	printVertex = cms.untracked.bool(True),
## 	printPtEtaPhi = cms.untracked.bool(True)
)

process.p = cms.Path(process.printTree)
