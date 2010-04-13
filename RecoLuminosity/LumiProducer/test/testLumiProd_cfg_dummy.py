import FWCore.ParameterSet.Config as cms

process = cms.Process("dbtest")

import FWCore.Framework.test.cmsExceptionsFatalOption_cff
process.options = cms.untracked.PSet(
#  wantSummary = cms.untracked.bool(True),
  Rethrow = FWCore.Framework.test.cmsExceptionsFatalOption_cff.Rethrow
)

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

process.source = cms.Source("EmptySource",
     numberEventsInRun = cms.untracked.uint32(10),
     firstRun = cms.untracked.uint32(124020),
     numberEventsInLuminosityBlock = cms.untracked.uint32(1),
     firstLuminosityBlock = cms.untracked.uint32(1)
)

process.DBService=cms.Service("DBService",
           authPath=cms.untracked.string('/afs/cern.ch/user/x/xiezhen')
)
process.lumiProducer=cms.EDProducer("LumiProducer",
#           connect=cms.string('oracle://devdb10/cms_xiezhen_dev'),
#           connect=cms.string('frontier://cmsfrontier.cern.ch:8000/LumiPrep/CMS_LUMI_DEV_OFFLINE'),
           connect=cms.string('frontier://LumiPrep/CMS_LUMI_DEV_OFFLINE'),                         
#           siteconfpath=cms.untracked.string('/afs/cern.ch/user/x/xiezhen/w1/lumical/CMSSW_3_5_0_pre5/src/RecoLuminosity/LumiProducer'),
           lumiversion=cms.untracked.string('0001') 
)
process.test = cms.EDAnalyzer("TestLumiProducer")
process.out = cms.OutputModule("PoolOutputModule",
           fileName=cms.untracked.string("testLumiProd.root")
)
process.p1 = cms.Path(process.lumiProducer * process.test)
process.e = cms.EndPath(process.out)
