import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts

'''

Configuration for 'shrinkingCone' PFTau Producer

See PFT-08-001 for a description of the algorithm.

'''

_shrinkingConeRecoTausConfig = cms.PSet(
    name = cms.string("shrinkingCone"),
    primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
    qualityCuts = PFTauQualityCuts.signalQualityCuts,
    # If true, consider PFLeptons (e/mu) as charged hadrons.
    usePFLeptons = cms.bool(True),
    pfCandSrc = cms.InputTag("particleFlow"),
    plugin = cms.string("RecoTauBuilderConePlugin"),
    leadObjectPt = cms.double(5.0),
    matchingCone = cms.string('0.1'),
    signalConeChargedHadrons = cms.string('min(max(5.0/et(), 0.07), 0.15)'),
    isoConeChargedHadrons = cms.string('0.5'),
    signalConePiZeros = cms.string('0.15'),
    isoConePiZeros = cms.string('0.5'),
    signalConeNeutralHadrons = cms.string('0.15'),
    isoConeNeutralHadrons = cms.string('0.5'),
)

_fixedConeRecoTausConfig = cms.PSet(
    name = cms.string("fixedCone"),
    primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
    qualityCuts = PFTauQualityCuts.signalQualityCuts,
    # If true, consider PFLeptons (e/mu) as charged hadrons.
    usePFLeptons = cms.bool(True),
    pfCandSrc = cms.InputTag("particleFlow"),
    plugin = cms.string("RecoTauBuilderConePlugin"),
    leadObjectPt = cms.double(5.0),
    matchingCone = cms.string('0.1'),
    signalConeChargedHadrons = cms.string('0.07'),
    isoConeChargedHadrons = cms.string('0.5'),
    signalConePiZeros = cms.string('0.15'),
    isoConePiZeros = cms.string('0.5'),
    signalConeNeutralHadrons = cms.string('0.10'),
    isoConeNeutralHadrons = cms.string('0.5'),
)

shrinkingConeRecoTaus = cms.EDProducer(
    "RecoTauProducer",
    jetSrc = cms.InputTag("ak5PFJets"),
    piZeroSrc = cms.InputTag("ak5PFJetsRecoTauPiZeros"),
    builders = cms.VPSet(
        _shrinkingConeRecoTausConfig
    ),
    # Build an empty tau in the case that a jet does not have any tracks
    buildNullTaus = cms.bool(True),
    modifiers = cms.VPSet(
        #cms.PSet(
            #name = cms.string('twoprong'),
            #plugin = cms.string("RecoTauTwoProngFilter"),
            #minPtFractionForSecondProng = cms.double(0.1),
        #),
        #cms.PSet(
            #name = cms.string('filterphotons'),
            #plugin = cms.string("RecoTauPhotonFilter"),
            #minPtFractionSinglePhotons  = cms.double(0.10),
            #minPtFractionPiZeroes       = cms.double(0.15),
        #),
        # Electron rejection
        cms.PSet(
            name = cms.string("shrinkingConeElectronRej"),
            plugin = cms.string("RecoTauElectronRejectionPlugin"),
            #Electron rejection parameters
            ElectronPreIDProducer                = cms.InputTag("elecpreid"),
            EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
            EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
            EcalStripSumE_minClusEnergy          = cms.double(0.1),
            EcalStripSumE_deltaEta               = cms.double(0.03),
            ElecPreIDLeadTkMatch_maxDR           = cms.double(0.01),
            maximumForElectrionPreIDOutput       = cms.double(-0.1),
            DataType = cms.string("AOD"),
        ),
    ),
)

fixedConeRecoTaus = shrinkingConeRecoTaus.clone(
    builders = cms.VPSet(
        _fixedConeRecoTausConfig
    )
)
