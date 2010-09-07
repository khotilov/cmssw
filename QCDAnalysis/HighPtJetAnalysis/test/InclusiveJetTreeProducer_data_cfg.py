import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('RecoJets.JetAssociationProducers.trackExtrapolator_cfi')
process.load('RecoJets.JetAssociationProducers.ak7JTA_cff') 
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2010A/JetMET/RECO/v4/000/144/114/F4E53E34-31B4-DF11-9A2C-0030487CD6D8.root'
    )
)
#############   Monster Event Removal  ###############
process.monster = cms.EDFilter(
"FilterOutScraping",
applyfilter = cms.untracked.bool(True),
debugOn     = cms.untracked.bool(False),
numtrack    = cms.untracked.uint32(10),
thresh      = cms.untracked.double(0.25)
)  
#############  Primary Vertex Filter  ################
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
       vertexCollection = cms.InputTag('offlinePrimaryVertices'),
       minimumNDOF = cms.uint32(4) ,
       maxAbsZ = cms.double(15),
       maxd0 = cms.double(2)
)
#############   inclusive jet tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('InclusiveJetTree.root'))

process.ak5calo = cms.EDAnalyzer('InclusiveJetTreeProducer',
    pfjets                 = cms.string('ak5PFJets'),
    jets                   = cms.string('ak5CaloJets'), 
    jetsID                 = cms.string('ak5JetID'),
    jetExtender            = cms.string('ak5JetExtender'),
    met                    = cms.string('met'),
    metNoHF                = cms.string('metNoHF'),
    hcalNoiseTag           = cms.InputTag('hcalnoise'),
    minJetPt               = cms.double(4.0), 
    minPFJetPt             = cms.double(4.0),
    jetTriggerNames        = cms.vstring('HLT_L1Jet6U','HLT_Jet15U','HLT_Jet30U','HLT_Jet50U'),
    triggerResultsTag      = cms.InputTag('TriggerResults','','HLT'),
    l1TriggerNames         = cms.vstring(),     
    L1GTReadoutRcdSource   = cms.InputTag('gtDigis'),
    L1GTObjectMapRcdSource = cms.InputTag('hltL1GtObjectMap','','HLT')
)

process.ak7calo = process.ak5calo.clone(pfjets = 'ak7PFJets',jets = 'ak7CaloJets', jetsID = 'ak7JetID', jetExtender = 'ak7JetExtender')

process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND NOT (36 OR 37 OR 38 OR 39)')

process.p = cms.Path(process.hltLevel1GTSeed * process.primaryVertexFilter * process.monster * process.trackExtrapolator *  process.ak7JTA 
* process.ak7calo)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10


