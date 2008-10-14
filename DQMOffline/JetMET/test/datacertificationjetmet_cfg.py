import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.DQM.collectorHost = ''

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1)
)

process.source = cms.Source("EmptySource"
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histos.root")
)

process.demo = cms.EDAnalyzer('DataCertificationJetMET',
                              fileName       = cms.untracked.string("/uscms/home/chlebana/DQM_V0001_R000063463__BeamHalo__BeamCommissioning08-PromptReco-v1__RECO.root"),
                              refFileName    = cms.untracked.string("/uscms/home/chlebana/DQM_V0001_R000063463__BeamHalo__BeamCommissioning08-PromptReco-v1__RECO.root"),
                              OutputFile     = cms.untracked.bool(True),
                              OutputFileName = cms.untracked.string("DQMResult.root")
)

process.p = cms.Path(process.demo)
