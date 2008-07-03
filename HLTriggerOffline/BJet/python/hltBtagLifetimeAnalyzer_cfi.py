import FWCore.ParameterSet.Config as cms

# configuration for HLTBtagLifetimeAnalyzer
hltBtagLifetimeAnalyzer = cms.EDAnalyzer("HLTBtagLifetimeAnalyzer",
    mcRadius = cms.double(0.1),
    outputFile = cms.string('plots.root'),
    computeStepEfficiencies = cms.bool(False),
    computeCumulativeEfficiencies = cms.bool(True),
    offlineRadius = cms.double(0.1), ## matching my pseudo-rapidity cone

    vertex = cms.InputTag("pixelVertices"),
    offlineBJets = cms.InputTag("jetProbabilityBJetTags"), ## match to offline btagged jets

    triggerPath = cms.string('HLT_BTagIP_Jet180'),
    levels = cms.VPSet(cms.PSet(
        filter = cms.InputTag("hltBLifetimeL1seeds","","HLT"),
        jets = cms.InputTag("hltIterativeCone5CaloJets","","HLT"),
        name = cms.string('L1'),
        title = cms.string('L1')
    ), 
        cms.PSet(
            filter = cms.InputTag("hltBLifetime1jetL2filter","","HLT"),
            jets = cms.InputTag("hltBLifetimeL25Jets","","HLT"),
            tracks = cms.InputTag("hltBLifetimeL25Associator","","HLT"),
            name = cms.string('L2'),
            title = cms.string('L2')
        ), 
        cms.PSet(
            filter = cms.InputTag("hltBLifetimeL25filter","","HLT"),
            jets = cms.InputTag("hltBLifetimeL3Jets","","HLT"),
            tracks = cms.InputTag("hltBLifetimeL3Associator","","HLT"),
            name = cms.string('L25'),
            title = cms.string('L2.5')
        ), 
        cms.PSet(
            filter = cms.InputTag("hltBLifetimeL3filter","","HLT"),
            jets = cms.InputTag("hltBLifetimeHLTJets"),
            name = cms.string('L3'),
            title = cms.string('L3')
        )),
    offlineCuts = cms.PSet(
        cut50 = cms.double(0.6), ## roughy 50% eff. in 1.7.0

        cut20 = cms.double(1.1), ## roughy 20% eff. in 1.7.0

        cut80 = cms.double(0.3) ## roughy 80% eff. in 1.7.0

    ),
    jetConfiguration = cms.PSet(
        maxEta = cms.double(5.0), ## pseudorapidity

        maxEnergy = cms.double(300.0) ## GeV

    ),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    mcFlavours = cms.PSet(
        light = cms.vuint32(1, 2, 3, 21), ## udsg

        c = cms.vuint32(4),
        b = cms.vuint32(5),
        g = cms.vuint32(21),
        uds = cms.vuint32(1, 2, 3)
    ),
    mcPartons = cms.InputTag("hltIC5byValAlgo"), ## pick hltIC5byValPhys or hltIC5byValAlgo

    vertexConfiguration = cms.PSet(
        maxZ = cms.double(20.0), ## cm

        maxR = cms.double(0.05) ## cm

    )
)
