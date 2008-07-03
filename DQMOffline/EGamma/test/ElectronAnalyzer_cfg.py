import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectronValidator")

process.load("RecoEcal.EgammaClusterProducers.geometryForClustering_cff")

process.load("DQMOffline.EGamma.electronAnalyzer_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cff")
process.load("EgammaAnalysis.EgammaIsolationProducers.egammaSuperClusterMerger_cfi")


DQMStore = cms.Service("DQMStore")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'rfio:/castor/cern.ch/cms/store/relval/2008/6/20/RelVal-RelValSingleGammaPt35-1213920853-IDEAL_V2-2nd/0000/28C3DDA9-B03E-DD11-AAA9-000423D9A212.root',
    'rfio:/castor/cern.ch/cms/store/relval/2008/6/20/RelVal-RelValSingleGammaPt35-1213920853-IDEAL_V2-2nd/0000/34A839BD-AF3E-DD11-ADFC-000423D99E46.root',
    'rfio:/castor/cern.ch/cms/store/relval/2008/6/20/RelVal-RelValSingleGammaPt35-1213920853-IDEAL_V2-2nd/0000/4ECD9899-B23E-DD11-90A9-001617C3B6E2.root',
    'rfio:/castor/cern.ch/cms/store/relval/2008/6/20/RelVal-RelValSingleGammaPt35-1213920853-IDEAL_V2-2nd/0000/CC46175D-B03E-DD11-8EE1-000423D9890C.root')

)


process.mergedSuperCluster = cms.EDFilter("SuperClusterMerger",
    src = cms.VInputTag(cms.InputTag("correctedHybridSuperClusters"), cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"))
)

process.FEVT = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring("keep *_MEtoEDMConverter_*_*"),
    #fileName = cms.untracked.string('pippo.root')
)

process.p1 = cms.Path(process.mergedSuperCluster*process.gsfElectronAnalysis*process.MEtoEDMConverter)
process.schedule = cms.Schedule(process.p1)

