# The following comments couldn't be translated into the new config version:

#TCL
#TCM
#TCT
#JBPL
#JBPM
#JBPT
#SMT
#SVM
#SVT
#CSVL 
#CSVM 
#CSVT

import FWCore.ParameterSet.Config as cms

Performance = cms.EDAnalyzer("PerformanceAnalyzer",
    # definition of the Operating Points (L,M,T)
    # cuts estimated either by thomas on 21X, or using old francisco's ones
    # sorted as TCL,TCM,TCT,JPL,JPM,JPT,JBPL, JBPM, JBPT, SMT, SVM, SVT, CSVL, CSVM, CSVT


#    bTagCutList = cms.untracked.vdouble(2.0, 4.6, 4.7, 0.26, 0.5, 
#        0.76, 1.2, 2.3, 3.2, 0.8, 
#        2.0, 3.6, 0.0, 37.0, 0.84, 
#        0.96),


                             bTagCutList = cms.untracked.VPSet(
    cms.PSet(
    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
    cut = cms.untracked.double(2.0),
    name = cms.untracked.string('TCL')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('trackCountingHighEffBJetTags'),
    cut = cms.untracked.double(4.2),
    name = cms.untracked.string('TCM')
    ),
    
    cms.PSet(
    collection = cms.untracked.InputTag('trackCountingHighPurBJetTags'),
    cut = cms.untracked.double(4.1),
    name = cms.untracked.string('TCT')
    ),
    
    cms.PSet(
    collection = cms.untracked.InputTag('jetProbabilityBJetTags'),
    cut = cms.untracked.double(0.24),
    name = cms.untracked.string('JPL')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('jetProbabilityBJetTags'),
    cut = cms.untracked.double(0.49),
    name = cms.untracked.string('JPM')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('jetProbabilityBJetTags'),
    cut = cms.untracked.double(0.74),
    name = cms.untracked.string('JPT')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('jetBProbabilityBJetTags'),
    cut = cms.untracked.double(1.1),
    name = cms.untracked.string('JBPL')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('jetBProbabilityBJetTags'),
    cut = cms.untracked.double(1.4),
    name = cms.untracked.string('JBPM')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('jetBProbabilityBJetTags'),
    cut = cms.untracked.double(1.4),
    name = cms.untracked.string('JBPT')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('softMuonBJetTags'),
    cut = cms.untracked.double(0.8),
    name = cms.untracked.string('SMT')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
    cut = cms.untracked.double(1.3),
    name = cms.untracked.string('SVL')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
    cut = cms.untracked.double(2.1),
    name = cms.untracked.string('SVM')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('simpleSecondaryVertexBJetTags'),
    cut = cms.untracked.double(3.6),
    name = cms.untracked.string('SVT')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('combinedSecondaryVertexBJetTags'),
    cut = cms.untracked.double(0.39),
    name = cms.untracked.string('CSVL')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('combinedSecondaryVertexBJetTags'),
    cut = cms.untracked.double(0.84),
    name = cms.untracked.string('CSVM')
    ),
    cms.PSet(
    collection = cms.untracked.InputTag('combinedSecondaryVertexBJetTags'),
    cut = cms.untracked.double(0.95),
    name = cms.untracked.string('CSVT')
    ),
    ),
                             
                             
                             #
                             # use jet corrections
                             #
                             useJetCorrections = cms.bool (True),
                             jetCorrectionsLabel = cms.string("L2L3JetCorrectorIcone5"),                     
    muoncuts = cms.PSet(
        MinNHits = cms.int32(7),
        MinMuonPt = cms.double(6.0),
        MaxMuonChi2 = cms.double(5.0),
        MaxMuonEta = cms.double(2.5)
    ),
    jetIdParameters = cms.PSet(
        vetoFlavour = cms.vstring(),
        rejectBCSplitting = cms.bool(False),
        physicsDefinition = cms.bool(False),
        coneSizeToAssociate = cms.double(0.3),
        fillLeptons = cms.bool(False),
        fillHeavyHadrons = cms.bool(False),
        fillPartons = cms.bool(True),
        mcSource = cms.string('source')
    ),
    Muons = cms.string('muons'),
    StoreTrackProba = cms.bool(False),
    #PSet jetIdParameters2 = {
    #       string mcSource = "source"
    #       bool fillPartons = true
    #       bool fillHeavyHadrons = false
    #       bool fillLeptons =  false
    #       double coneSizeToAssociate = 0.3
    #       bool physicsDefinition = true
    #       bool rejectBCSplitting = true
    #       vstring vetoFlavour = {  }
    #}
    #
    # prepare pset for TrackHistory
    #
    jetcuts = cms.PSet(
        MaxEta = cms.double(2.5),
        MinDeltaR = cms.double(0.4),
        MinPt = cms.double(20.0),
        MinPtRel = cms.double(-1.0)
    ),
    bTagTrackEventIPtagInfos = cms.string(''),
    # bTagTrackEventIPtagInfos = cms.string('impactParameterTagInfos'),
    flavourMatchOption = cms.string('genParticles'),
    WeightHistograms = cms.bool(False),
    TrackCollection = cms.untracked.string('generalTracks'),
    GenJets = cms.string('iterativeCone5GenJets'),
    Jets = cms.string('iterativeCone5CaloJets'),
    bTagTrackEvent = cms.bool(False),
    #InputTag simG4 = g4SimHits
    outputFile = cms.untracked.string('results.root'),
    StorePtHat = cms.bool(False),
    SimTracks = cms.string('g4SimHits'),
    AwayJetTagger = cms.string('TCL'),
    flavourSource = cms.InputTag("IC5byValAlgo"),
    StoreWeightsInNtuple = cms.bool(False),
    PrimaryVertexCollection = cms.untracked.string('offlinePrimaryVerticesFromCTFTracks'),
    WritePerformancePlots = cms.bool(True)
)



