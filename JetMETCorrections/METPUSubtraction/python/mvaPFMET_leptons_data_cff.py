
import FWCore.ParameterSet.Config as cms

#from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff  import *
from JetMETCorrections.Configuration.DefaultJEC_cff                     import *
from JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cfi            import *
from RecoJets.JetProducers.PileupJetIDParams_cfi                        import JetIdParams

calibratedAK5PFJetsForPFMEtMVA = cms.EDProducer('PFJetCorrectionProducer',
    src = cms.InputTag('ak5PFJets'),
    correctors = cms.vstring("ak5PFL1FastL2L3Residual") for Data
)

pfMEtMVA = cms.EDProducer("PFMETProducerMVA",
    srcCorrJets = cms.InputTag('calibratedAK5PFJetsForPFMEtMVA'),
    srcUncorrJets = cms.InputTag('ak5PFJets'),
    srcPFCandidates = cms.InputTag('particleFlow'),
    srcVertices = cms.InputTag('offlinePrimaryVertices'),
    srcLeptons = cms.VInputTag("isomuons","isoelectrons","isotaus"),#"muons","hpsPFTauProducer"), # NOTE: you need to set this to collections of electrons, muons and tau-jets
                          #       passing the lepton reconstruction & identification criteria applied in your analysis
    srcRho = cms.InputTag('kt6PFJets','rho'),
    globalThreshold = cms.double(-1.),#pfMet.globalThreshold,
    minCorrJetPt = cms.double(-1.),
    inputFileNames = cms.PSet(
        U     = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmet_52.root'),
        DPhi  = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetphi_52.root'),
        CovU1 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetu1cov_52.root'),
        CovU2 = cms.FileInPath('JetMETCorrections/METPUSubtraction/data/gbrmetu2cov_52.root')
    ),
    corrector = cms.string("ak5PFL1Fastjet"),
    useType1 = cms.bool(False),
    dZcut = cms.double(0.1),
    impactParTkThreshold = cms.double(0.),
    tmvaWeights = cms.string("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml"),
    tmvaMethod = cms.string("JetID"),
    version = cms.int32(-1),
    cutBased = cms.bool(False),                      
    tmvaVariables = cms.vstring(
        "nvtx",
        "jetPt",
        "jetEta",
        "jetPhi",
        "dZ",
        "d0",
        "beta",
        "betaStar",
        "nCharged",
        "nNeutrals",
        "dRMean",
        "frac01",
        "frac02",
        "frac03",
        "frac04",
        "frac05",
    ),
    tmvaSpectators = cms.vstring(),
    JetIdParams = JetIdParams,
    label = cms.string("PhilV1"),
    verbosity = cms.int32(0)
)

pfMEtMVAsequence  = cms.Sequence(
    (isomuonseq+isotauseq+isoelectronseq)*
    calibratedAK5PFJetsForPFMEtMVA*
    pfMEtMVA
    )

