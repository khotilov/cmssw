import FWCore.ParameterSet.Config as cms

# load jet energy correction parameters
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

# load config file for producing collection of ak5PFJets
# using pfNoPileUp collection of PFCandidates as input
from JetMETCorrections.Type1MET.ak5PFchsJets_cff import ak5PFchsJetsSequence

#--------------------------------------------------------------------------------
# select PFCandidates ("unclustered energy") not within jets
# for Type 2 MET correction
from CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi import pfNoJet
pfchsCandsNotInJet = pfNoJet.clone(
    topCollection = cms.InputTag('ak5PFchsJets'),
    bottomCollection = cms.InputTag('particleFlow')
)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 1 + 2 MET corrections for PFJets
pfchsJetMETcorr = cms.EDProducer("PFJetMETcorrInputProducer",
    src = cms.InputTag('ak5PFchsJets'),
    offsetCorrLabel = cms.string("ak5PFchsL1Fastjet"),
    jetCorrLabel = cms.string("ak5PFchsL1FastL2L3"), # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)                                         
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 2 MET corrections for selected PFCandidates
pfchsCandMETcorr = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag('pfchsCandsNotInJet')                         
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce Type 0 MET corrections for selected vertices
pfchsMETcorrCHS = cms.EDProducer("PFchsMETcorrInputProducer",
    src = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof = cms.uint32(4),
    goodVtxZ = cms.double(24)
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# use MET corrections to produce Type 1 / Type 1 + 2 corrected PFMET objects
pfchsType1CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(True),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfchsJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)   

pfchsType1p2CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(True),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfchsJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(True),
    srcUnclEnergySums = cms.VInputTag(
        cms.InputTag('pfchsJetMETcorr', 'type2'),
        cms.InputTag('pfchsJetMETcorr', 'offset'),
        cms.InputTag('pfchsCandMETcorr')                                    
    ),                              
    type2CorrFormula = cms.string("A"),
    type2CorrParameter = cms.PSet(
        A = cms.double(1.4)
    )
)   
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define sequence to run all modules
producePFMETchsCorrections = cms.Sequence(
    ak5PFchsJetsSequence
   + pfchsCandsNotInJet
   + pfchsJetMETcorr
   + pfchsCandMETcorr
   + pfchsMETcorr
   + pfchsType1CorrectedMet
   + pfchsType1p2CorrectedMet
)
#--------------------------------------------------------------------------------
