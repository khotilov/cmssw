# Auto generated configuration file
# using: 
# Revision: 1.99.2.3 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/PYTHIA6_Bs2JpsiPhi_10TeV_cff.py -s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions FrontierConditions_GlobalTag,IDEAL_V11::All -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/MixingNoPileUp_cff')
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('Configuration/StandardSequences/VtxSmearedEarly10TeVCollision_cff')
process.load('Configuration/StandardSequences/Sim_cff')
process.load('Configuration/StandardSequences/Digi_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/L1TriggerDefaultMenu_cff')
process.load('Configuration/StandardSequences/DigiToRaw_cff')
process.load('HLTrigger/Configuration/HLT_2E30_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1. $'),
    annotation = cms.untracked.string('BstoJpsiPhi_mumuKK'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/GenProduction/python/PYTHIA6_Bs2JpsiPhi_10TeV_cff.py,v $')
)

process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring('BsToJpsiPhi_uptoHLT_1.log'),
                                    simul = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
                                    )

################################################
# grab all starting seeds from /dev/urandom
################################################

process.load("IOMC.RandomEngine.IOMC_cff")
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randHelper = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PythiaSource",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    comEnergy = cms.untracked.double(10000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTJ(11)=3     ! Choice of the fragmentation function', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(51)=10042     ! CTEQ6L1 structure function chosen', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model', 
            'MSTU(21)=1     ! Check on possible errors during program execution', 
            'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
            'PARP(89)=1960. ! sqrts for which PARP82 is set', 
            'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
            'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
            'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
            'PARP(67)=2.5    ! amount of initial-state radiation', 
            'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
            'PARP(62)=1.25   ! ', 
            'PARP(64)=0.2    ! ', 
            'MSTP(91)=1     !', 
            'PARP(91)=2.1   ! kt distribution', 
            'PARP(93)=15.0  ! '),
        processParameters = cms.vstring('MDCY(134,1) = 0', 
            'MDCY(137,1) = 0', 
            'MDCY(138,1) = 0', 
            'MDCY(135,1) = 0', 
            'MDCY(141,1) = 0', 
            'MDCY(140,1) = 0', 
            'MDCY(15,1) = 0', 
            'MDCY(123,1) = 0', 
            'MDCY(126,1) = 0', 
            'MDCY(129,1) = 0', 
            'MDCY(122,1) = 0', 
            'MDCY(125,1) = 0', 
            'MDCY(128,1) = 0', 
            'MDCY(262,1) = 0', 
            'MDCY(264,1) = 0', 
            'MDCY(263,1) = 0', 
            'MDCY(265,1) = 0', 
            'MDCY(286,1) = 0', 
            'MDCY(287,1) = 0', 
            'MDCY(124,1) = 0', 
            'MDCY(127,1) = 0', 
            'MDCY(266,1) = 0', 
            'MDCY(288,1) = 0', 
            'MDCY(267,1) = 0', 
            'MDCY(130,1) = 0', 
            'MDCY(112,1) = 0', 
            'MDCY(113,1) = 0', 
            'MDCY(114,1) = 0', 
            'MDCY(117,1) = 0', 
            'MDCY(258,1) = 0', 
            'MDCY(256,1) = 0', 
            'MDCY(257,1) = 0', 
            'MDCY(259,1) = 0', 
            'MDCY(284,1) = 0', 
            'MDCY(283,1) = 0', 
            'MDCY(118,1) = 0', 
            'MDCY(115,1) = 0', 
            'MDCY(102,1) = 0', 
            'MDCY(109,1) = 0', 
            'MDCY(103,1) = 0', 
            'MDCY(107,1) = 0', 
            'MDCY(110,1) = 0', 
            'MDCY(119,1) = 0', 
            'MDCY(120,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(280,1) = 0', 
            'MDCY(281,1) = 0', 
            'MDCY(108,1) = 0', 
            'MDCY(104,1) = 0', 
            'MDCY(253,1) = 0', 
            'MDCY(251,1) = 0', 
            'MDCY(250,1) = 0', 
            'MDCY(252,1) = 0', 
            'MDCY(254,1) = 0', 
            'MDCY(282,1) = 0', 
            'MDCY(285,1) = 0', 
            'MDCY(111,1) = 0', 
            'MDCY(121,1) = 0', 
            'MDCY(255,1) = 0', 
            'MDCY(261,1) = 0', 
            'MDCY(131,1) = 0', 
            'MDCY(132,1) = 0', 
            'MDCY(295,1) = 0', 
            'MDCY(268,1) = 0', 
            'MDCY(289,1) = 0', 
            'MDCY(133,1) = 0', 
            'MDCY(146,1) = 0', 
            'MDCY(147,1) = 0', 
            'MDCY(296,1) = 0', 
            'MDCY(278,1) = 0', 
            'MDCY(294,1) = 0', 
            'MDCY(148,1) = 0', 
            'MDCY(279,1) = 0', 
            'MDCY(181,1) = 0', 
            'MDCY(182,1) = 0', 
            'MDCY(84,1) = 0', 
            'MDCY(179,1) = 0', 
            'MDCY(185,1) = 0', 
            'MDCY(189,1) = 0', 
            'MDCY(187,1) = 0', 
            'MDCY(194,1) = 0', 
            'MDCY(192,1) = 0', 
            'MDCY(164,1) = 0', 
            'MDCY(169,1) = 0', 
            'MDCY(158,1) = 0', 
            'MDCY(159,1) = 0', 
            'MDCY(175,1) = 0', 
            'MDCY(155,1) = 0', 
            'MDCY(151,1) = 0', 
            'MDCY(162,1) = 0', 
            'MDCY(167,1) = 0', 
            'MDCY(163,1) = 0', 
            'MDCY(170,1) = 0', 
            'MDCY(168,1) = 0', 
            'MDCY(174,1) = 0', 
            'MDCY(172,1) = 0', 
            'MDCY(173,1) = 0', 
            'MDCY(176,1) = 0', 
            'MDCY(180,1) = 0', 
            'MDCY(186,1) = 0', 
            'MDCY(188,1) = 0', 
            'MDCY(193,1) = 0', 
            'MDCY(195,1) = 0', 
            'MDCY(196,1) = 0', 
            'MDCY(197,1) = 0', 
            'MDCY(43,1) = 0', 
            'MDCY(44,1) = 0', 
            'MDCY(269,1) = 0', 
            'MDCY(210,1) = 0', 
            'MDCY(211,1) = 0', 
            'MDCY(219,1) = 0', 
            'MDCY(227,1) = 0', 
            'MDCY(217,1) = 0', 
            'MDCY(208,1) = 0', 
            'MDCY(215,1) = 0', 
            'MDCY(143,1) = 0', 
            'MDCY(223,1) = 0', 
            'MDCY(225,1) = 0', 
            'MDCY(272,1) = 0', 
            'MDCY(291,1) = 0', 
            'MDCY(273,1) = 0', 
            'MDCY(139,1) = 0', 
            'MDCY(270,1) = 0', 
            'MDCY(290,1) = 0', 
            'MDCY(271,1) = 0', 
            'MDCY(136,1) = 0', 
            'MDCY(274,1) = 0', 
            'MDCY(292,1) = 0', 
            'MDCY(275,1) = 0', 
            'MDCY(142,1) = 0', 
            'MDCY(144,1) = 0', 
            'MDCY(145,1) = 0', 
            'MDCY(209,1) = 0', 
            'MDCY(218,1) = 0', 
            'MDCY(216,1) = 0', 
            'MDCY(224,1) = 0', 
            'MDCY(226,1) = 0', 
            'MDCY(228,1) = 0', 
            'MDCY(276,1) = 0', 
            'MDCY(277,1) = 0', 
            'MDCY(293,1) = 0', 
            'MDCY(105,1) = 0', 
            'MSEL=1         ! All process'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('BsToJpsiPhi_uptoHLT_1.root'),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

process.evtgenproducer = cms.EDProducer("EvtGenProducer",
    use_default_decay = cms.untracked.bool(False),
    decay_table = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/DECAY.DEC'),
    particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
    user_decay_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/Bs_Jpsiphi_mumuKK.dec'),
    list_forced_decays = cms.vstring('MyB_s0', 
        'Myanti-B_s0'),
    processParameters = cms.vstring('MDCY(134,1) = 0', 
        'MDCY(137,1) = 0', 
        'MDCY(138,1) = 0', 
        'MDCY(135,1) = 0', 
        'MDCY(141,1) = 0', 
        'MDCY(140,1) = 0', 
        'MDCY(15,1) = 0', 
        'MDCY(123,1) = 0', 
        'MDCY(126,1) = 0', 
        'MDCY(129,1) = 0', 
        'MDCY(122,1) = 0', 
        'MDCY(125,1) = 0', 
        'MDCY(128,1) = 0', 
        'MDCY(262,1) = 0', 
        'MDCY(264,1) = 0', 
        'MDCY(263,1) = 0', 
        'MDCY(265,1) = 0', 
        'MDCY(286,1) = 0', 
        'MDCY(287,1) = 0', 
        'MDCY(124,1) = 0', 
        'MDCY(127,1) = 0', 
        'MDCY(266,1) = 0', 
        'MDCY(288,1) = 0', 
        'MDCY(267,1) = 0', 
        'MDCY(130,1) = 0', 
        'MDCY(112,1) = 0', 
        'MDCY(113,1) = 0', 
        'MDCY(114,1) = 0', 
        'MDCY(117,1) = 0', 
        'MDCY(258,1) = 0', 
        'MDCY(256,1) = 0', 
        'MDCY(257,1) = 0', 
        'MDCY(259,1) = 0', 
        'MDCY(284,1) = 0', 
        'MDCY(283,1) = 0', 
        'MDCY(118,1) = 0', 
        'MDCY(115,1) = 0', 
        'MDCY(102,1) = 0', 
        'MDCY(109,1) = 0', 
        'MDCY(103,1) = 0', 
        'MDCY(107,1) = 0', 
        'MDCY(110,1) = 0', 
        'MDCY(119,1) = 0', 
        'MDCY(120,1) = 0', 
        'MDCY(281,1) = 0', 
        'MDCY(280,1) = 0', 
        'MDCY(281,1) = 0', 
        'MDCY(108,1) = 0', 
        'MDCY(104,1) = 0', 
        'MDCY(253,1) = 0', 
        'MDCY(251,1) = 0', 
        'MDCY(250,1) = 0', 
        'MDCY(252,1) = 0', 
        'MDCY(254,1) = 0', 
        'MDCY(282,1) = 0', 
        'MDCY(285,1) = 0', 
        'MDCY(111,1) = 0', 
        'MDCY(121,1) = 0', 
        'MDCY(255,1) = 0', 
        'MDCY(261,1) = 0', 
        'MDCY(131,1) = 0', 
        'MDCY(132,1) = 0', 
        'MDCY(295,1) = 0', 
        'MDCY(268,1) = 0', 
        'MDCY(289,1) = 0', 
        'MDCY(133,1) = 0', 
        'MDCY(146,1) = 0', 
        'MDCY(147,1) = 0', 
        'MDCY(296,1) = 0', 
        'MDCY(278,1) = 0', 
        'MDCY(294,1) = 0', 
        'MDCY(148,1) = 0', 
        'MDCY(279,1) = 0', 
        'MDCY(181,1) = 0', 
        'MDCY(182,1) = 0', 
        'MDCY(84,1) = 0', 
        'MDCY(179,1) = 0', 
        'MDCY(185,1) = 0', 
        'MDCY(189,1) = 0', 
        'MDCY(187,1) = 0', 
        'MDCY(194,1) = 0', 
        'MDCY(192,1) = 0', 
        'MDCY(164,1) = 0', 
        'MDCY(169,1) = 0', 
        'MDCY(158,1) = 0', 
        'MDCY(159,1) = 0', 
        'MDCY(175,1) = 0', 
        'MDCY(155,1) = 0', 
        'MDCY(151,1) = 0', 
        'MDCY(162,1) = 0', 
        'MDCY(167,1) = 0', 
        'MDCY(163,1) = 0', 
        'MDCY(170,1) = 0', 
        'MDCY(168,1) = 0', 
        'MDCY(174,1) = 0', 
        'MDCY(172,1) = 0', 
        'MDCY(173,1) = 0', 
        'MDCY(176,1) = 0', 
        'MDCY(180,1) = 0', 
        'MDCY(186,1) = 0', 
        'MDCY(188,1) = 0', 
        'MDCY(193,1) = 0', 
        'MDCY(195,1) = 0', 
        'MDCY(196,1) = 0', 
        'MDCY(197,1) = 0', 
        'MDCY(43,1) = 0', 
        'MDCY(44,1) = 0', 
        'MDCY(269,1) = 0', 
        'MDCY(210,1) = 0', 
        'MDCY(211,1) = 0', 
        'MDCY(219,1) = 0', 
        'MDCY(227,1) = 0', 
        'MDCY(217,1) = 0', 
        'MDCY(208,1) = 0', 
        'MDCY(215,1) = 0', 
        'MDCY(143,1) = 0', 
        'MDCY(223,1) = 0', 
        'MDCY(225,1) = 0', 
        'MDCY(272,1) = 0', 
        'MDCY(291,1) = 0', 
        'MDCY(273,1) = 0', 
        'MDCY(139,1) = 0', 
        'MDCY(270,1) = 0', 
        'MDCY(290,1) = 0', 
        'MDCY(271,1) = 0', 
        'MDCY(136,1) = 0', 
        'MDCY(274,1) = 0', 
        'MDCY(292,1) = 0', 
        'MDCY(275,1) = 0', 
        'MDCY(142,1) = 0', 
        'MDCY(144,1) = 0', 
        'MDCY(145,1) = 0', 
        'MDCY(209,1) = 0', 
        'MDCY(218,1) = 0', 
        'MDCY(216,1) = 0', 
        'MDCY(224,1) = 0', 
        'MDCY(226,1) = 0', 
        'MDCY(228,1) = 0', 
        'MDCY(276,1) = 0', 
        'MDCY(277,1) = 0', 
        'MDCY(293,1) = 0', 
        'MDCY(105,1) = 0', 
        'MSEL=1         ! All Process')
)

# Other statements
process.GlobalTag.globaltag = 'IDEAL_V11::All'

process.bfilter = cms.EDFilter("PythiaFilter",
    moduleLabel = cms.untracked.string('source'),
    ParticleID = cms.untracked.int32(5)
)
process.decayfilter = cms.EDFilter("BdecayFilter",
    secondDaughterDecayPtMin = cms.double(0.0),
    moduleLabel = cms.untracked.string('evtgenproducer'),
    secondDaughterDecayEtaMax = cms.double(10.0),
    firstDaughterDecayPtMin = cms.double(2.5),
    motherParticle = cms.int32(531),
    secondDaughterDecayEtaMin = cms.double(-10.0),
    firstDaughterDecayEtaMin = cms.double(-2.5),
    secondDaughter = cms.int32(333),
    secondDaughterDecay = cms.vint32(321, 321),
    firstDaughterDecayEtaMax = cms.double(2.5),
    firstDaughter = cms.int32(443),
    firstDaughterDecay = cms.vint32(13, 13),
)

process.VtxSmeared.src = 'evtgenproducer'
process.genParticleCandidates.src = 'evtgenproducer'
process.g4SimHits.Generator.HepMCProductLabel = 'evtgenproducer'
process.genParticles.src = 'evtgenproducer'
process.genParticles.abortOnUnknownPDGCode = False

process.ProductionFilterSequence = cms.Sequence(process.bfilter*process.evtgenproducer*process.decayfilter)

# Path and EndPath definitions
process.generation_step = cms.Path(process.ProductionFilterSequence*process.pgen)
process.simulation_step = cms.Path(process.ProductionFilterSequence*process.psim)
process.digitisation_step = cms.Path(process.ProductionFilterSequence*process.pdigi)
process.L1simulation_step = cms.Path(process.ProductionFilterSequence*process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.ProductionFilterSequence*process.DigiToRaw)
process.HLT_TripleJet60_MET60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sTripleJet60MET60+process.hltPreTripleJet60MET60+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt3jet60+process.HLTEndSequence))
process.HLT_BTagMu_QuadJet30_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuQuadJet30Relaxed+process.HLTBCommonL2recoSequence+process.hltBSoftmuon4jetL2filter30+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filterRelaxed+process.HLTEndSequence))
process.HLT_MinBiasPixel_Trk5 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMinBiasPixel+process.hltPreMinBiasPixelTrk5+process.HLTDoLocalPixelSequence+process.HLTPixelTrackingForMinBiasSequence+process.hltPixelCands+process.hltPixelMBForAlignment+process.HLTEndSequence))
process.HLT_IsoEle12_Jet40 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEJet+process.hltPreIsoEle12Jet40+process.HLTEJetElectronSequence+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltej1jet40+process.HLTEndSequence))
process.HLT_Jet250 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet250+process.hltPreJet250+process.HLTRecoJetRegionalSequence+process.hlt1jet250+process.HLTEndSequence))
process.HLT_IsoMu15 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuIso10+process.hltPreIsoMu15+process.hltSingleMuIsoL1Filtered10+process.HLTL2muonrecoSequence+process.hltSingleMuIsoL2PreFiltered12+process.HLTL2muonisorecoSequence+process.hltSingleMuIsoL2IsoFiltered12+process.HLTL3muonrecoSequence+process.hltSingleMuIsoL3PreFiltered15+process.HLTL3muonisorecoSequence+process.hltSingleMuIsoL3IsoFiltered15+process.HLTEndSequence))
process.HLT_DoubleEle5_SW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedDoubleEgammaEt5+process.hltPreDoubleEle5SWL1R+process.HLTDoubleElectronEt5L1NonIsoHLTnonIsoSequence+process.HLTEndSequence))
process.HLT_IsoMu11 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuIso7+process.hltPreIsoMu11+process.hltSingleMuIsoL1Filtered+process.HLTL2muonrecoSequence+process.hltSingleMuIsoL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltSingleMuIsoL2IsoFiltered+process.HLTL3muonrecoSequence+process.hltSingleMuIsoL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltSingleMuIsoL3IsoFiltered+process.HLTEndSequence))
process.HLT_MET75 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET75+process.hltPreMET75+process.HLTRecoJetMETSequence+process.hlt1MET75+process.HLTEndSequence))
process.AlCa_HcalPhiSym = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sAlCaHcalPhiSym+process.hltPreAlCaHcalPhiSym+process.hltAlCaHcalFEDSelector+process.HLTEndSequence))
process.HLT_Ele15_SW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt12+process.hltPreEle15SWL1R+process.HLTSingleElectronEt15L1NonIsoHLTNonIsoSequence+process.HLTEndSequence))
process.HLT_IsoMu7_BTagIP_Jet35 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuB+process.hltPreIsoMu7BTagIPJet35+process.hltMuBLifetimeL1Filtered+process.HLTL2muonrecoSequence+process.hltMuBLifetimeIsoL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltMuBLifetimeIsoL2IsoFiltered+process.HLTBCommonL2recoSequence+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTL3muonrecoSequence+process.hltMuBLifetimeIsoL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltMuBLifetimeIsoL3IsoFiltered+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.HLT_DiJetAve50 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve50+process.hltPreDiJetAve50+process.HLTRecoJetMETSequence+process.hltdijetave50+process.HLTEndSequence))
process.HLT_IsoEle12_IsoTau_Trk3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sElectronTau+process.hltPreIsoEle12IsoTauTrk3+process.HLTETauSingleElectronL1IsolatedHOneOEMinusOneOPFilterSequence+process.HLTL2TauJetsElectronTauSequnce+process.hltL2ElectronTauIsolationProducer+process.hltL2ElectronTauIsolationSelector+process.hltFilterEcalIsolatedTauJetsElectronTau+process.HLTRecopixelvertexingSequence+process.hltJetTracksAssociatorAtVertexL25ElectronTau+process.hltConeIsolationL25ElectronTau+process.hltIsolatedTauJetsSelectorL25ElectronTauPtLeadTk+process.hltFilterIsolatedTauJetsL25ElectronTauPtLeadTk+process.hltIsolatedTauJetsSelectorL25ElectronTau+process.hltFilterIsolatedTauJetsL25ElectronTau+process.HLTEndSequence))
process.HLT_DiJetAve220 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve220+process.hltPreDiJetAve220+process.HLTRecoJetMETSequence+process.hltdijetave220+process.HLTEndSequence))
process.HLT_IsoTau_MET35_Trk15_L1MET = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleTauMET+process.hltPreIsoTauMET35Trk15L1MET+process.HLTCaloTausCreatorSequence+process.hltMet+process.hlt1METSingleTauMET+process.hltL2SingleTauMETJets+process.hltL2SingleTauMETIsolationProducer+process.hltL2SingleTauMETIsolationSelector+process.hltFilterSingleTauMETEcalIsolation+process.HLTDoLocalPixelSequence+process.HLTRecopixelvertexingSequence+process.hltAssociatorL25SingleTauMET+process.hltConeIsolationL25SingleTauMET+process.hltIsolatedL25SingleTauMET+process.hltFilterL25SingleTauMET+process.HLTDoLocalStripSequence+process.hltL3SingleTauMETPixelSeeds+process.hltCkfTrackCandidatesL3SingleTauMET+process.hltCtfWithMaterialTracksL3SingleTauMET+process.hltAssociatorL3SingleTauMET+process.hltConeIsolationL3SingleTauMET+process.hltIsolatedL3SingleTauMET+process.hltFilterL3SingleTauMET+process.HLTEndSequence))
process.HLT_DoubleIsoEle12_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedDoubleEgamma+process.hltPreDoubleIsoEle12L1R+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoDoubleElectronL1MatchFilterRegional+process.hltL1NonIsoDoubleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsolatedElectronHcalIsol+process.hltL1NonIsoDoubleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.HLTPixelMatchElectronL1NonIsoSequence+process.hltL1NonIsoDoubleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.HLTPixelMatchElectronL1NonIsoTrackingSequence+process.hltL1NonIsoDoubleElectronEoverpFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.HLTL1NonIsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1NonIsoElectronTrackIsol+process.hltL1NonIsoDoubleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_MET25 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET25+process.hltPreMET25+process.HLTRecoJetMETSequence+process.hlt1MET25+process.HLTEndSequence))
process.HLT_DoubleMu4_BJPsi = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJpsitoMumu+process.hltPreDoubleMu4BJPsi+process.hltJpsitoMumuL1Filtered+process.HLTL2muonrecoSequence+process.HLTL3displacedMumurecoSequence+process.hltDisplacedJpsitoMumuFilter+process.HLTEndSequence))
process.HLTriggerFinalPath = cms.Path(process.ProductionFilterSequence*(process.hltTriggerSummaryAOD+process.hltPreTriggerSummaryRAW+process.hltTriggerSummaryRAW+process.hltBoolFinalPath))
process.HLT_Mu11 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso7+process.hltPreMu11+process.hltSingleMuNoIsoL1Filtered7+process.HLTL2muonrecoSequence+process.hltSingleMuNoIsoL2PreFiltered9+process.HLTL3muonrecoSequence+process.hltSingleMuNoIsoL3PreFiltered11+process.HLTEndSequence))
process.HLT_BTagMu_QuadJet40 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuQuadJet40+process.HLTBCommonL2recoSequence+process.hltBSoftmuon4jetL2filter+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filter+process.HLTEndSequence))
process.HLT_IsoEle18_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgamma+process.hltPreIsoEle18L1R+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoSingleElectronL1MatchFilterRegional+process.hltL1NonIsoSingleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsolatedElectronHcalIsol+process.hltL1NonIsoSingleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.HLTPixelMatchElectronL1NonIsoSequence+process.hltL1NonIsoSingleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.HLTPixelMatchElectronL1NonIsoTrackingSequence+process.hltL1NonIsoSingleElectronHOneOEMinusOneOPFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.HLTL1NonIsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1NonIsoElectronTrackIsol+process.hltL1NonIsoSingleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_Mu5_TripleJet30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuNoIsoJets30+process.hltPreMu5TripleJet30+process.hltMuNoIsoJetsMinPt4L1Filtered+process.HLTL2muonrecoSequence+process.hltMuNoIsoJetsMinPt4L2PreFiltered+process.HLTL3muonrecoSequence+process.hltMuNoIsoJetsMinPtL3PreFiltered+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltMuNoIsoHLTJets3jet30+process.HLTEndSequence))
process.HLT_IsoEle10_BTagIP_Jet35 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sElectronB+process.hltPreIsoEle10BTagIPJet35+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltElBElectronL1MatchFilter+process.hltElBElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsolatedElectronHcalIsol+process.hltElBElectronHcalIsolFilter+process.HLTBCommonL2recoSequence+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1NonIsoSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltElBElectronPixelMatchFilter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.HLTPixelMatchElectronL1NonIsoTrackingSequence+process.hltElBElectronEoverpFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.HLTL1NonIsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1NonIsoElectronTrackIsol+process.hltElBElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_MET65 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET65+process.hltPreMET65+process.HLTRecoJetMETSequence+process.hlt1MET65+process.HLTEndSequence))
process.HLT_Ele15_LW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt10+process.hltPreEle15LWL1R+process.HLTSingleElectronLWEt15L1NonIsoHLTNonIsoSequence+process.HLTEndSequence))
process.HLT_BTagMu_TripleJet40_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuTripleJet40Relaxed+process.HLTBCommonL2recoSequence+process.hltBSoftmuon3jetL2filter40+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filterRelaxed+process.HLTEndSequence))
process.AlCa_IsoTrack = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sAlCaIsoTrack+process.hltPreAlCaIsoTrack+process.HLTL3PixelIsolFilterSequence+process.HLTIsoTrRegFEDSelection+process.HLTEndSequence))
process.HLT_EM200 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgamma+process.hltPreEM200+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoSingleEMVeryHighEtL1MatchFilterRegional+process.hltL1NonIsoSinglePhotonEMVeryHighEtEtFilter+process.HLTEndSequence))
process.HLT_DoubleJet60_MET60_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet60MET60Aco+process.hltPreDoubleJet60MET60Aco+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hltPhi2metAco+process.HLTEndSequence))
process.HLT_DoubleIsoMu3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiMuonIso+process.hltPreDoubleIsoMu3+process.hltDiMuonIsoL1Filtered+process.HLTL2muonrecoSequence+process.hltDiMuonIsoL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltDiMuonIsoL2IsoFiltered+process.HLTL3muonrecoSequence+process.hltDiMuonIsoL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltDiMuonIsoL3IsoFiltered+process.HLTEndSequence))
process.HLT_Jet30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet30+process.hltPreJet30+process.HLTRecoJetMETSequence+process.hlt1jet30+process.HLTEndSequence))
process.HLT_DoubleMu3_Psi2S = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJpsiMM+process.hltPreDoubleMu3Psi2S+process.hltJpsiMML1Filtered+process.HLTL2muonrecoSequence+process.hltPsi2SMML2Filtered+process.HLTL3muonrecoSequence+process.hltPsi2SMML3Filtered+process.HLTEndSequence))
process.HLT_CSCBeamHaloOverlapRing1 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sCSCBeamHaloOverlapRing1+process.hltPreCSCBeamHaloOverlapRing1+process.hltMuonCSCDigis+process.hltCsc2DRecHits+process.hltOverlapsHLTCSCBeamHaloOverlapRing1+process.HLTEndSequence))
process.HLT_CSCBeamHaloOverlapRing2 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sCSCBeamHaloOverlapRing2+process.hltPreCSCBeamHaloOverlapRing2+process.hltMuonCSCDigis+process.hltCsc2DRecHits+process.hltOverlapsHLTCSCBeamHaloOverlapRing2+process.HLTEndSequence))
process.HLT_QuadJet35_MET60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sQuadJet35MET60+process.hltPreQuadJet35MET60+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt4jet35+process.HLTEndSequence))
process.HLT_TripleMu3_TauTo3Mu = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuMuk+process.hltPreTripleMu3TauTo3Mu+process.hltMuMukL1Filtered+process.HLTL2muonrecoSequence+process.HLTL3displacedMumurecoSequence+process.hltDisplacedMuMukFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTRecopixelvertexingSequence+process.hltMumukPixelSeedFromL2Candidate+process.hltCkfTrackCandidatesMumuk+process.hltCtfWithMaterialTracksMumuk+process.hltMumukAllConeTracks+process.hltmmkFilter+process.HLTEndSequence))
process.HLT_DoubleIsoEle12_LW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedDoubleEgamma+process.hltPreDoubleIsoEle12LWL1R+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoLargeWindowDoubleElectronL1MatchFilterRegional+process.hltL1NonIsoLargeWindowDoubleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsolatedElectronHcalIsol+process.hltL1NonIsoLargeWindowDoubleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoLargeWindowSequence+process.HLTPixelMatchElectronL1NonIsoLargeWindowSequence+process.hltL1NonIsoLargeWindowDoubleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoLargeWindowTrackingSequence+process.HLTPixelMatchElectronL1NonIsoLargeWindowTrackingSequence+process.hltL1NonIsoLargeWindowDoubleElectronEoverpFilter+process.HLTL1IsoLargeWindowElectronsRegionalRecoTrackerSequence+process.HLTL1NonIsoLargeWindowElectronsRegionalRecoTrackerSequence+process.hltL1IsoLargeWindowElectronTrackIsol+process.hltL1NonIsoLargeWindowElectronTrackIsol+process.hltL1NonIsoLargeWindowDoubleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_BTagIP_Jet180 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetime+process.hltPreBTagIPJet180+process.HLTBCommonL2recoSequence+process.hltBLifetime1jetL2filter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.AlCa_EcalPhiSym = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sAlCaEcalPhiSym+process.hltPreAlCaEcalPhiSym+process.hltEcalDigis+process.hltEcalWeightUncalibRecHit+process.hltEcalRecHit+process.hltAlCaPhiSymStream+process.HLTEndSequence))
process.HLT_DoubleEle10_Z = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleEgamma+process.hltPreDoubleEle10Z+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoDoubleElectronZeeL1MatchFilterRegional+process.hltL1IsoDoubleElectronZeeEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoDoubleElectronZeeHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltL1IsoDoubleElectronZeePixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.hltL1IsoDoubleElectronZeeEoverpFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1IsoDoubleElectronZeeTrackIsolFilter+process.hltL1IsoDoubleElectronZeePMMassFilter+process.HLTEndSequence))
process.HLT_Photon15_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt10+process.hltPrePhoton15L1R+process.HLTSinglePhoton15L1NonIsolatedHLTNonIsoSequence+process.HLTEndSequence))
process.HLT_IsoPhoton25_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt15+process.hltPreIsoPhoton25L1R+process.HLTSinglePhoton25L1NonIsolatedHLTIsoSequence+process.HLTEndSequence))
process.HLT_TrackerCosmics = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sTrackerCosmics+process.hltPreTrackerCosmics+process.HLTEndSequence))
process.HLT_IsoEle12_DoubleJet80 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEJet+process.hltPreIsoEle12DoubleJet80+process.HLTEJetElectronSequence+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltej2jet80+process.HLTEndSequence))
process.HLT_DiJetAve70 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve70+process.hltPreDiJetAve70+process.HLTRecoJetMETSequence+process.hltdijetave70+process.HLTEndSequence))
process.HLT_QuadJet60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sQuadJet60+process.hltPreQuadJet60+process.HLTRecoJetMETSequence+process.hlt4jet60+process.HLTEndSequence))
process.HLT_DoubleJet50_MET70_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet50MET70Aco+process.hltPreDoubleJet50MET70Aco+process.HLTRecoJetMETSequence+process.hlt1MET70+process.hltPhiJet2metAco+process.HLTEndSequence))
process.HLT_L2Mu9 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso7+process.hltPreL2Mu9+process.hltSingleMuNoIsoL1Filtered7+process.HLTL2muonrecoSequence+process.hltSingleMuLevel2NoIsoL2PreFiltered+process.HLTEndSequence))
process.HLT_DoubleEle6_Exclusive = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sExclusiveDoubleEgamma+process.hltPreDoubleEle6Exclusive+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoDoubleExclElectronL1MatchFilterRegional+process.hltL1IsoDoubleExclElectronEtPhiFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoDoubleExclElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltL1IsoDoubleExclElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.hltL1IsoDoubleExclElectronEoverpFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1IsoDoubleExclElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_Photon25_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt15+process.hltPrePhoton25L1R+process.HLTSinglePhoton25L1NonIsolatedHLTNonIsoSequence+process.HLTEndSequence))
process.HLT_DoubleMu3_SameSign = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSameSignMu+process.hltPreDoubleMu3SameSign+process.hltSameSignMuL1Filtered+process.HLTL2muonrecoSequence+process.hltSameSignMuL2PreFiltered+process.HLTL3muonrecoSequence+process.hltSameSignMuL3PreFiltered+process.HLTEndSequence))
process.HLT_IsoMu13 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuIso10+process.hltPreIsoMu13+process.hltSingleMuIsoL1Filtered10+process.HLTL2muonrecoSequence+process.hltSingleMuIsoL2PreFiltered11+process.HLTL2muonisorecoSequence+process.hltSingleMuIsoL2IsoFiltered11+process.HLTL3muonrecoSequence+process.hltSingleMuIsoL3PreFiltered13+process.HLTL3muonisorecoSequence+process.hltSingleMuIsoL3IsoFiltered13+process.HLTEndSequence))
process.HLT_IsoMu14_IsoTau_Trk3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuonTau+process.hltPreIsoMu14IsoTauTrk3+process.hltMuonTauL1Filtered+process.HLTL2muonrecoSequence+process.hltMuonTauIsoL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltMuonTauIsoL2IsoFiltered+process.HLTDoLocalStripSequence+process.HLTL3muonrecoSequence+process.HLTL3muonisorecoSequence+process.hltMuonTauIsoL3PreFiltered+process.hltMuonTauIsoL3IsoFiltered+process.HLTCaloTausCreatorRegionalSequence+process.hltL2TauJetsProviderMuonTau+process.hltL2MuonTauIsolationProducer+process.hltL2MuonTauIsolationSelector+process.hltFilterEcalIsolatedTauJetsMuonTau+process.HLTDoLocalPixelSequence+process.HLTRecopixelvertexingSequence+process.hltJetsPixelTracksAssociatorMuonTau+process.hltPixelTrackConeIsolationMuonTau+process.hltIsolatedTauJetsSelectorL25MuonTauPtLeadTk+process.hltFilterL25MuonTauPtLeadTk+process.hltPixelTrackIsolatedTauJetsSelectorMuonTau+process.hltFilterPixelTrackIsolatedTauJetsMuonTau+process.HLTEndSequence))
process.HLT_DoubleMu7_Z = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sZMM+process.hltPreDoubleMu7Z+process.hltZMML1Filtered+process.HLTL2muonrecoSequence+process.hltZMML2Filtered+process.HLTL3muonrecoSequence+process.hltZMML3Filtered+process.HLTEndSequence))
process.HLT_Jet100_MET60_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet100MET60Aco+process.hltPreJet100MET60Aco+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt1jet100+process.hlt1jet1METAco+process.HLTEndSequence))
process.HLT_BTagIP_HT320_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetimeLowEnergy+process.hltPreBTagIPHT320Relaxed+process.HLTBCommonL2recoSequence+process.hltBLifetimeHTL2filter320+process.HLTBLifetimeL25recoSequenceRelaxed+process.hltBLifetimeL25filterRelaxed+process.HLTBLifetimeL3recoSequenceRelaxed+process.hltBLifetimeL3filterRelaxed+process.HLTEndSequence))
process.HLT_BTagMu_DoubleJet60_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBtagMuDoubleJet60Relaxed+process.HLTBCommonL2recoSequence+process.hltBSoftmuon2jetL2filter60+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filterRelaxed+process.HLTEndSequence))
process.HLT_DoubleFwdJet40_MET60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleFwdJet40MET60+process.hltPreDoubleFwdJet40MET60+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt2jetvbf+process.HLTEndSequence))
process.HLT_DiJetAve15 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve15+process.hltPreDiJetAve15+process.HLTRecoJetMETSequence+process.hltdijetave15+process.HLTEndSequence))
process.HLT_DoubleEle10_LW_OnlyPixelM_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedDoubleEgammaEt5+process.hltPreDoubleEle10LWOnlyPixelML1R+process.HLTDoubleElectronLWonlyPMEt10L1NonIsoHLTNonIsoSequence+process.HLTEndSequence))
process.HLT_MET65_HT350 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET35HT350+process.hltPreMET35HT350+process.HLTRecoJetMETSequence+process.hlt1MET65+process.hlt1HT350+process.HLTEndSequence))
process.HLT_MinBiasHcal = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMinBiasHcal+process.hltPreMinBiasHcal+process.HLTEndSequence))
process.HLT_L1Jet15 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sL1Jet15+process.hltPreL1Jet15+process.HLTEndSequence))
process.HLT_BTagIP_Jet120_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetimeLowEnergy+process.hltPreBTagIPJet120Relaxed+process.HLTBCommonL2recoSequence+process.hltBLifetime1jetL2filter120+process.HLTBLifetimeL25recoSequenceRelaxed+process.hltBLifetimeL25filterRelaxed+process.HLTBLifetimeL3recoSequenceRelaxed+process.hltBLifetimeL3filterRelaxed+process.HLTEndSequence))
process.HLT_DoubleIsoTau_Trk3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleTau+process.hltPreDoubleIsoTauTrk3+process.HLTCaloTausCreatorRegionalSequence+process.hltL2DoubleTauJets+process.hltL2DoubleTauIsolationProducer+process.hltL2DoubleTauIsolationSelector+process.hltFilterDoubleTauEcalIsolation+process.HLTDoLocalPixelSequence+process.HLTRecopixelvertexingSequence+process.hltAssociatorL25PixelTauIsolated+process.hltConeIsolationL25PixelTauIsolated+process.hltIsolatedL25PixelTauPtLeadTk+process.hltFilterL25PixelTauPtLeadTk+process.hltIsolatedL25PixelTau+process.hltFilterL25PixelTau+process.HLTEndSequence))
process.HLT_L1Mu = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sL1Mu+process.hltPreL1Mu+process.hltMuLevel1PathL1Filtered+process.HLTEndSequence))
process.HLT_BTagIP_QuadJet40 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetime+process.hltPreBTagIPQuadJet40+process.HLTBCommonL2recoSequence+process.hltBLifetime4jetL2filter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.HLT_BTagIP_DoubleJet120 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetime+process.hltPreBTagIPDoubleJet120+process.HLTBCommonL2recoSequence+process.hltBLifetime2jetL2filter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.HLT_Jet50 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet50+process.hltPreJet50+process.HLTRecoJetMETSequence+process.hlt1jet50+process.HLTEndSequence))
process.HLT_LooseIsoEle15_LW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt12+process.hltPreLooseIsoEle15LWL1R+process.HLTSingleElectronLWEt15L1NonIsoHLTLooseIsoSequence+process.HLTEndSequence))
process.HLTriggerFirstPath = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltGetRaw+process.hltPreFirstPath+process.hltBoolFirstPath+process.HLTEndSequence))
process.HLT_BTagIP_TripleJet40_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetimeLowEnergy+process.hltPreBTagIPTripleJet40Relaxed+process.HLTBCommonL2recoSequence+process.hltBLifetime3jetL2filter40+process.HLTBLifetimeL25recoSequenceRelaxed+process.hltBLifetimeL25filterRelaxed+process.HLTBLifetimeL3recoSequenceRelaxed+process.hltBLifetimeL3filterRelaxed+process.HLTEndSequence))
process.HLT_Ele10_SW_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt8+process.hltPreEle10SWL1R+process.HLTSingleElectronEt10L1NonIsoHLTnonIsoSequence+process.HLTEndSequence))
process.HLT_BTagIP_QuadJet30_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetimeLowEnergy+process.hltPreBTagIPQuadJet30Relaxed+process.HLTBCommonL2recoSequence+process.hltBLifetime4jetL2filter30+process.HLTBLifetimeL25recoSequenceRelaxed+process.hltBLifetimeL25filterRelaxed+process.HLTBLifetimeL3recoSequenceRelaxed+process.hltBLifetimeL3filterRelaxed+process.HLTEndSequence))
process.HLT_Mu13 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso10+process.hltPreMu13+process.hltSingleMuNoIsoL1Filtered10+process.HLTL2muonrecoSequence+process.hltSingleMuNoIsoL2PreFiltered11+process.HLTL3muonrecoSequence+process.hltSingleMuNoIsoL3PreFiltered13+process.HLTEndSequence))
process.HLT_IsoPhoton40_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgamma+process.hltPreIsoPhoton40L1R+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoSinglePhotonL1MatchFilterRegional+process.hltL1NonIsoSinglePhotonEtFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1NonIsolatedPhotonEcalIsol+process.hltL1NonIsoSinglePhotonEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedPhotonHcalIsol+process.hltL1NonIsolatedPhotonHcalIsol+process.hltL1NonIsoSinglePhotonHcalIsolFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.HLTL1NonIsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1NonIsoPhotonTrackIsol+process.hltL1NonIsoSinglePhotonTrackIsolFilter+process.HLTEndSequence))
process.HLT_DiJetAve130 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve130+process.hltPreDiJetAve130+process.HLTRecoJetMETSequence+process.hltdijetave130+process.HLTEndSequence))
process.HLT_Mu15 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso10+process.hltPreMu15+process.hltSingleMuNoIsoL1Filtered10+process.HLTL2muonrecoSequence+process.hltSingleMuNoIsoL2PreFiltered12+process.HLTL3muonrecoSequence+process.hltSingleMuNoIsoL3PreFiltered15+process.HLTEndSequence))
process.HLT_L1MuOpen = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sL1MuOpen+process.hltPreL1MuOpen+process.hltMuLevel1PathL1OpenFiltered+process.HLTEndSequence))
process.HLT_CSCBeamHaloRing2or3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sCSCBeamHaloRing2or3+process.hltPreCSCBeamHaloRing2or3+process.hltMuonCSCDigis+process.hltCsc2DRecHits+process.hltFilter23HLTCSCBeamHaloRing2or3+process.HLTEndSequence))
process.HLT_DoubleIsoEle10_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleEgamma+process.hltPreDoubleIsoEle10L1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoDoubleElectronL1MatchFilterRegional+process.hltL1IsoDoubleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoDoubleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltL1IsoDoubleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.hltL1IsoDoubleElectronEoverpFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1IsoDoubleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_MinBiasPixel = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMinBiasPixel+process.hltPreMinBiasPixel+process.HLTDoLocalPixelSequence+process.HLTPixelTrackingForMinBiasSequence+process.hltPixelCands+process.hltMinBiasPixelFilter+process.HLTEndSequence))
process.HLT_DoubleMu3_BJPsi = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJpsitoMumuRelaxed+process.hltPreDoubleMu3BJPsi+process.hltJpsitoMumuL1FilteredRelaxed+process.HLTL2muonrecoSequence+process.HLTL3displacedMumurecoSequence+process.hltDisplacedJpsitoMumuFilterRelaxed+process.HLTEndSequence))
process.HLT_IsoEle15_LW_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleEgamma+process.hltPreIsoEle15LWL1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoLargeWindowSingleL1MatchFilter+process.hltL1IsoLargeWindowSingleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoLargeWindowSingleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoLargeWindowSequence+process.hltL1IsoLargeWindowSingleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoLargeWindowTrackingSequence+process.hltL1IsoLargeWindowSingleElectronHOneOEMinusOneOPFilter+process.HLTL1IsoLargeWindowElectronsRegionalRecoTrackerSequence+process.hltL1IsoLargeWindowElectronTrackIsol+process.hltL1IsoLargeWindowSingleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_MET35 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET35+process.hltPreMET35+process.HLTRecoJetMETSequence+process.hlt1MET35+process.HLTEndSequence))
process.HLT_DoubleMu3_JPsi = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJpsiMM+process.hltPreDoubleMu3JPsi+process.hltJpsiMML1Filtered+process.HLTL2muonrecoSequence+process.hltJpsiMML2Filtered+process.HLTL3muonrecoSequence+process.hltJpsiMML3Filtered+process.HLTEndSequence))
process.HLT_ForwardBSC = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sForwardBSC+process.hltPreForwardBSC+process.HLTEndSequence))
process.HLT_CSCBeamHalo = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sCSCBeamHalo+process.hltPreCSCBeamHalo+process.HLTEndSequence))
process.HLT_IsoEle15_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleEgamma+process.hltPreIsoEle15L1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoSingleL1MatchFilter+process.hltL1IsoSingleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoSingleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltL1IsoSingleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.hltL1IsoSingleElectronHOneOEMinusOneOPFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1IsoSingleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_DoubleMu3_Vtx2mm = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiMuonNoIso+process.hltPreDoubleMu3Vtx2mm+process.hltDiMuonNoIsoL1Filtered+process.HLTL2muonrecoSequence+process.hltDiMuonNoIsoL2PreFiltered+process.HLTL3muonrecoSequence+process.hltDiMuonNoIsoL3PreFilteredRelaxedVtx2mm+process.HLTEndSequence))
process.HLT_Jet110 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet110+process.hltPreJet110+process.HLTRecoJetRegionalSequence+process.hlt1jet110+process.HLTEndSequence))
process.HLT_DoubleMu3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiMuonNoIso+process.hltPreDoubleMu3+process.hltDiMuonNoIsoL1Filtered+process.HLTL2muonrecoSequence+process.hltDiMuonNoIsoL2PreFiltered+process.HLTL3muonrecoSequence+process.hltDiMuonNoIsoL3PreFiltered+process.HLTEndSequence))
process.HLT_DoubleJet40_MET70_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet40MET70Aco+process.hltPreDoubleJet40MET70Aco+process.HLTRecoJetMETSequence+process.hlt1MET70+process.hltPhiJet1Jet2Aco+process.HLTEndSequence))
process.HLT_DoubleJet125_MET60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet125MET60+process.hltPreDoubleJet125MET60+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt2jet125New+process.HLTEndSequence))
process.HLT_L1MET20 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sL1MET20+process.hltPreL1MET20+process.HLTEndSequence))
process.HLT_Jet180_MET60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet180MET60+process.hltPreJet180MET60+process.HLTRecoJetMETSequence+process.hlt1MET60+process.hlt1jet180+process.HLTEndSequence))
process.HLT_IsoEle5_TripleJet30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEJet30+process.hltPreIsoEle5TripleJet30+process.HLTE3Jet30ElectronSequence+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltej3jet30+process.HLTEndSequence))
process.HLT_NoL2IsoMu8_Jet40 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuNoL2IsoJets+process.hltPreNoL2IsoMu8Jet40+process.hltMuNoL2IsoJetsL1Filtered+process.HLTL2muonrecoSequence+process.hltMuNoL2IsoJetsL2PreFiltered+process.HLTL3muonrecoSequence+process.hltMuNoL2IsoJetsL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltMuNoL2IsoJetsL3IsoFiltered+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltMuNoL2IsoJetsHLT1jet40+process.HLTEndSequence))
process.HLT_IsoEle10_Mu10_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEgMuNonIso+process.hltPreIsoEle10Mu10L1R+process.hltNonIsoEMuL1MuonFilter+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltemuNonIsoL1MatchFilterRegional+process.hltemuNonIsoL1IsoEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsolatedElectronHcalIsol+process.hltemuNonIsoL1HcalIsolFilter+process.HLTL2muonrecoSequence+process.hltNonIsoEMuL2MuonPreFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.HLTPixelMatchElectronL1NonIsoSequence+process.hltemuNonIsoL1IsoPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.HLTPixelMatchElectronL1NonIsoTrackingSequence+process.hltemuNonIsoL1IsoEoverpFilter+process.HLTL3muonrecoSequence+process.hltNonIsoEMuL3MuonPreFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.HLTL1NonIsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltL1NonIsoElectronTrackIsol+process.hltemuNonIsoL1IsoTrackIsolFilter+process.HLTEndSequence))
process.HLT_IsoPhoton20_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt15+process.hltPreIsoPhoton20L1R+process.HLTSinglePhoton20L1NonIsolatedHLTIsoSequence+process.HLTEndSequence))
process.HLT_QuadJet30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sQuadJet30+process.hltPreQuadJet30+process.HLTRecoJetMETSequence+process.hlt4jet30+process.HLTEndSequence))
process.HLT_BackwardBSC = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBackwardBSC+process.hltPreBackwardBSC+process.HLTEndSequence))
process.HLT_IsoMu7_Jet40 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuJets+process.hltPreIsoMu7Jet40+process.hltMuJetsL1Filtered+process.HLTL2muonrecoSequence+process.hltMuJetsL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltMuJetsL2IsoFiltered+process.HLTL3muonrecoSequence+process.hltMuJetsL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltMuJetsL3IsoFiltered+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltMuJetsHLT1jet40+process.HLTEndSequence))
process.HLT_FwdJet20 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sFwdJet20+process.hltPreFwdJet20+process.HLTRecoJetMETSequence+process.hltRapGap+process.HLTEndSequence))
process.HLT_DoublePhoton10_Exclusive = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sExclusiveDoubleEgamma+process.hltPreDoublePhoton10Exclusive+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoDoubleExclPhotonL1MatchFilterRegional+process.hltL1IsoDoubleExclPhotonEtPhiFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1IsoDoubleExclPhotonEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedPhotonHcalIsol+process.hltL1IsoDoubleExclPhotonHcalIsolFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1IsoDoubleExclPhotonTrackIsolFilter+process.HLTEndSequence))
process.HLT_MinBiasEcal = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMinBiasEcal+process.hltPreMinBiasEcal+process.HLTEndSequence))
process.HLT_MinBias = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMinBias+process.hltPreMinBias+process.HLTEndSequence))
process.AlCa_EcalPi0 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sAlCaEcalPi0+process.hltPreAlCaEcalPi0+process.HLTDoRegionalEgammaEcalSequence+process.hltAlCaPi0RegRecHits+process.HLTEndSequence))
process.HLT_SumET120 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSumET120+process.hltPreSumET120+process.HLTRecoJetMETSequence+process.hlt1SumET120+process.HLTEndSequence))
process.HLT_DoubleMu3_Upsilon = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sUpsilonMM+process.hltPreDoubleMu3Upsilon+process.hltUpsilonMML1Filtered+process.HLTL2muonrecoSequence+process.hltUpsilonMML2Filtered+process.HLTL3muonrecoSequence+process.hltUpsilonMML3Filtered+process.HLTEndSequence))
process.HLT_DoubleIsoPhoton20_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedDoubleEgamma+process.hltPreDoubleIsoPhoton20L1R+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoDoublePhotonL1MatchFilterRegional+process.hltL1NonIsoDoublePhotonEtFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1NonIsolatedPhotonEcalIsol+process.hltL1NonIsoDoublePhotonEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedPhotonHcalIsol+process.hltL1NonIsolatedPhotonHcalIsol+process.hltL1NonIsoDoublePhotonHcalIsolFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.HLTL1NonIsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1NonIsoPhotonTrackIsol+process.hltL1NonIsoDoublePhotonTrackIsolFilter+process.hltL1NonIsoDoublePhotonDoubleEtFilter+process.HLTEndSequence))
process.HLT_BTagMu_HT250_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonHTLowEnergy+process.hltPreBTagMuHT250Relaxed+process.HLTBCommonL2recoSequence+process.hltBSoftmuonHTL2filter250+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filterRelaxed+process.HLTEndSequence))
process.HLT_Jet80 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet80+process.hltPreJet80+process.HLTRecoJetRegionalSequence+process.hlt1jet80+process.HLTEndSequence))
process.HLT_IsoEle12_QuadJet35 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEJet+process.hltPreIsoEle12QuadJet35+process.HLTEJetElectronSequence+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltej4jet35+process.HLTEndSequence))
process.HLT_Jet180 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet180+process.hltPreJet180+process.HLTRecoJetRegionalSequence+process.hlt1jet180regional+process.HLTEndSequence))
process.HLT_DiJetAve30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDiJetAve30+process.hltPrediJetAve30+process.HLTRecoJetMETSequence+process.hltdijetave30+process.HLTEndSequence))
process.HLT_BTagMu_Jet20_Calib = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuJet20Calib+process.HLTBCommonL2recoSequence+process.hltBSoftmuon1jetL2filter+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonByDRL3filter+process.HLTEndSequence))
process.HLT_BTagMu_DoubleJet120 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuDoubleJet120+process.HLTBCommonL2recoSequence+process.hltBSoftmuon2jetL2filter+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filter+process.HLTEndSequence))
process.HLT_IsoMu9 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuIso7+process.hltPreIsoMu9+process.hltSingleMuIsoL1Filtered+process.HLTL2muonrecoSequence+process.hltSingleMuIsoL2PreFiltered7+process.HLTL2muonisorecoSequence+process.hltSingleMuIsoL2IsoFiltered7+process.HLTL3muonrecoSequence+process.hltSingleMuIsoL3PreFiltered9+process.HLTL3muonisorecoSequence+process.hltSingleMuIsoL3IsoFiltered9+process.HLTEndSequence))
process.HLT_IsoTau_MET65_Trk20 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleTau+process.hltPreIsoTauMET65Trk20+process.HLTCaloTausCreatorSequence+process.hltMet+process.hlt1METSingleTau+process.hltL2SingleTauJets+process.hltL2SingleTauIsolationProducer+process.hltL2SingleTauIsolationSelector+process.hltFilterSingleTauEcalIsolation+process.HLTDoLocalPixelSequence+process.HLTRecopixelvertexingSequence+process.hltAssociatorL25SingleTau+process.hltConeIsolationL25SingleTau+process.hltIsolatedL25SingleTau+process.hltFilterL25SingleTau+process.HLTDoLocalStripSequence+process.hltL3SingleTauPixelSeeds+process.hltCkfTrackCandidatesL3SingleTau+process.hltCtfWithMaterialTracksL3SingleTau+process.hltAssociatorL3SingleTau+process.hltConeIsolationL3SingleTau+process.hltIsolatedL3SingleTau+process.hltFilterL3SingleTau+process.HLTEndSequence))
process.HLT_BTagIP_DoubleJet60_Relaxed = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetimeLowEnergy+process.hltPreBTagIPDoubleJet60Relaxed+process.HLTBCommonL2recoSequence+process.hltBLifetime2jetL2filter60+process.HLTBLifetimeL25recoSequenceRelaxed+process.hltBLifetimeL25filterRelaxed+process.HLTBLifetimeL3recoSequenceRelaxed+process.hltBLifetimeL3filterRelaxed+process.HLTEndSequence))
process.HLT_MET50 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMET50+process.hltPreMET50+process.HLTRecoJetMETSequence+process.hlt1MET50+process.HLTEndSequence))
process.HLT_DoubleLooseIsoTau = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleTauRelaxed+process.hltPreDoubleLooseIsoTau+process.HLTCaloTausCreatorRegionalSequence+process.hltL2DoubleTauJetsRelaxed+process.hltL2DoubleTauIsolationProducerRelaxed+process.hltL2DoubleTauIsolationSelectorRelaxed+process.hltFilterDoubleTauEcalIsolationRelaxed+process.HLTEndSequence))
process.HLT_EM80 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgamma+process.hltPreEM80+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.HLTL1NonIsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1NonIsoRecoEcalCandidate+process.hltL1NonIsoSingleEMHighEtL1MatchFilterRegional+process.hltL1NonIsoSinglePhotonEMHighEtEtFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1NonIsolatedPhotonEcalIsol+process.hltL1NonIsoSingleEMHighEtEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1NonIsolatedElectronHcalIsol+process.hltL1IsolatedElectronHcalIsol+process.hltL1NonIsoSingleEMHighEtHOEFilter+process.hltHcalDoubleCone+process.hltL1NonIsoEMHcalDoubleCone+process.hltL1NonIsoSingleEMHighEtHcalDBCFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.HLTL1NonIsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1NonIsoPhotonTrackIsol+process.hltL1NonIsoSingleEMHighEtTrackIsolFilter+process.HLTEndSequence))
process.HLT_IsoEle12_TripleJet60 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sEJet+process.hltPreIsoEle12TripleJet60+process.HLTEJetElectronSequence+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltej3jet60+process.HLTEndSequence))
process.HLT_Jet60_MET70_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sJet60MET70Aco+process.hltPreJet60MET70Aco+process.HLTRecoJetMETSequence+process.hlt1MET70+process.hltPhiJet1metAco+process.HLTEndSequence))
process.HLT_DoubleJet150 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet150+process.hltPreDoubleJet150+process.HLTRecoJetRegionalSequence+process.hlt2jet150+process.HLTEndSequence))
process.HLT_IsoMu7_BTagMu_Jet20 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuB+process.hltPreIsoMu7BTagMuJet20+process.hltMuBSoftL1Filtered+process.HLTL2muonrecoSequence+process.hltMuBSoftIsoL2PreFiltered+process.HLTL2muonisorecoSequence+process.hltMuBSoftIsoL2IsoFiltered+process.HLTBCommonL2recoSequence+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTL3muonrecoSequence+process.hltMuBSoftIsoL3PreFiltered+process.HLTL3muonisorecoSequence+process.hltMuBSoftIsoL3IsoFiltered+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filter+process.HLTEndSequence))
process.HLT_IsoPhoton10_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt8+process.hltPreIsoPhoton10L1R+process.HLTSinglePhotonEt10L1NonIsolatedSequence+process.HLTEndSequence))
process.HLT_DoubleJet125_Aco = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleJet125Aco+process.hltPreDoubleJet125Aco+process.HLTRecoJetRegionalSequence+process.hlt2jet125+process.hlt2jetAco+process.HLTEndSequence))
process.HLT_BTagMu_TripleJet70 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftmuonNjet+process.hltPreBTagMuTripleJet70+process.HLTBCommonL2recoSequence+process.hltBSoftmuon3jetL2filter+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filter+process.HLTEndSequence))
process.HLT_IsoEle8_IsoMu7 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sIsoEgMu+process.hltPreIsoEle8IsoMu7+process.hltEMuL1MuonFilter+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltemuL1IsoSingleL1MatchFilter+process.hltemuL1IsoSingleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltemuL1IsoSingleElectronHcalIsolFilter+process.HLTL2muonrecoSequence+process.hltEMuL2MuonPreFilter+process.HLTL2muonisorecoSequence+process.hltEMuL2MuonIsoFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoSequence+process.hltemuL1IsoSingleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoTrackingSequence+process.hltemuL1IsoSingleElectronEoverpFilter+process.HLTL3muonrecoSequence+process.hltEMuL3MuonPreFilter+process.HLTL3muonisorecoSequence+process.hltEMuL3MuonIsoFilter+process.HLTL1IsoElectronsRegionalRecoTrackerSequence+process.hltL1IsoElectronTrackIsol+process.hltemuL1IsoSingleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_Mu14_Jet50 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sMuNoIsoJets+process.hltPreMu14Jet50+process.hltMuNoIsoJetsL1Filtered+process.HLTL2muonrecoSequence+process.hltMuNoIsoJetsL2PreFiltered+process.HLTL3muonrecoSequence+process.hltMuNoIsoJetsL3PreFiltered+process.HLTDoCaloSequence+process.HLTDoJetRecoSequence+process.hltMuNoIsoJetsHLT1jet50+process.HLTEndSequence))
process.HLT_Mu15_Vtx2mm = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso7+process.hltPreMu15Vtx2mm+process.hltSingleMuNoIsoL1Filtered7+process.HLTL2muonrecoSequence+process.hltSingleMuNoIsoL2PreFiltered12L1pre7+process.HLTL3muonrecoSequence+process.hltSingleMuNoIsoL3PreFilteredRelaxedVtx2mm+process.HLTEndSequence))
process.HLT_LooseIsoTau_MET30_L1MET = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleTauMET+process.hltPreLooseIsoTauMET30L1MET+process.HLTCaloTausCreatorSequence+process.hltMet+process.hlt1METSingleTauMETRelaxed+process.hltL2SingleTauMETJets+process.hltL2SingleTauMETIsolationProducer+process.hltL2SingleTauMETIsolationSelectorRelaxed+process.hltFilterSingleTauMETEcalIsolationRelaxed+process.HLTEndSequence))
process.HLT_BTagIP_HT470 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetime+process.hltPreBTagIPHT470+process.HLTBCommonL2recoSequence+process.hltBLifetimeHTL2filter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.HLT_DoubleIsoEle10_LW_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleEgamma+process.hltPreDoubleIsoEle10LWL1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoLargeWindowDoubleElectronL1MatchFilterRegional+process.hltL1IsoLargeWindowDoubleElectronEtFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedElectronHcalIsol+process.hltL1IsoLargeWindowDoubleElectronHcalIsolFilter+process.HLTDoLocalPixelSequence+process.HLTDoLocalStripSequence+process.HLTPixelMatchElectronL1IsoLargeWindowSequence+process.hltL1IsoLargeWindowDoubleElectronPixelMatchFilter+process.HLTPixelMatchElectronL1IsoLargeWindowTrackingSequence+process.hltL1IsoLargeWindowDoubleElectronEoverpFilter+process.HLTL1IsoLargeWindowElectronsRegionalRecoTrackerSequence+process.hltL1IsoLargeWindowElectronTrackIsol+process.hltL1IsoLargeWindowDoubleElectronTrackIsolFilter+process.HLTEndSequence))
process.HLT_LooseIsoTau_MET30 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleTau+process.hltPreLooseIsoTauMET30+process.HLTCaloTausCreatorSequence+process.hltMet+process.hlt1METSingleTauRelaxed+process.hltL2SingleTauJets+process.hltL2SingleTauIsolationProducer+process.hltL2SingleTauIsolationSelectorRelaxed+process.hltFilterSingleTauEcalIsolationRelaxed+process.HLTEndSequence))
process.HLT_TripleJet85 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sTripleJet85+process.hltPreTripleJet85+process.HLTRecoJetRegionalSequence+process.hlt3jet85+process.HLTEndSequence))
process.HLT_BTagIP_TripleJet70 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBLifetime+process.hltPreBTagIPTripleJet70+process.HLTBCommonL2recoSequence+process.hltBLifetime3jetL2filter+process.HLTBLifetimeL25recoSequence+process.hltBLifetimeL25filter+process.HLTBLifetimeL3recoSequence+process.hltBLifetimeL3filter+process.HLTEndSequence))
process.HLT_Mu3 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuPrescale3+process.hltPreMu3+process.hltSingleMuPrescale3L1Filtered+process.HLTL2muonrecoSequence+process.hltSingleMuPrescale3L2PreFiltered+process.HLTL3muonrecoSequence+process.hltSingleMuPrescale3L3PreFiltered+process.HLTEndSequence))
process.HLT_Mu5 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuPrescale5+process.hltPreMu5+process.hltSingleMuPrescale5L1Filtered+process.HLTL2muonrecoSequence+process.hltSingleMuPrescale5L2PreFiltered+process.HLTL3muonrecoSequence+process.hltSingleMuPrescale5L3PreFiltered+process.HLTEndSequence))
process.HLT_Mu7 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuPrescale77+process.hltPreMu7+process.hltSingleMuPrescale77L1Filtered+process.HLTL2muonrecoSequence+process.hltSingleMuPrescale77L2PreFiltered+process.HLTL3muonrecoSequence+process.hltSingleMuPrescale77L3PreFiltered+process.HLTEndSequence))
process.HLT_Mu9 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleMuNoIso7+process.hltPreMu9+process.hltSingleMuNoIsoL1Filtered7+process.HLTL2muonrecoSequence+process.hltSingleMuNoIsoL2PreFiltered7+process.HLTL3muonrecoSequence+process.hltSingleMuNoIsoL3PreFiltered9+process.HLTEndSequence))
process.HLT_DoubleFwdJet50 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleFwdJet50+process.hltPreDoubleFwdJet50+process.HLTRecoJetMETSequence+process.hlt2jetGapFilter+process.HLTEndSequence))
process.HLT_IsoPhoton30_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sSingleEgamma+process.hltPreIsoPhoton30L1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoSinglePhotonL1MatchFilter+process.hltL1IsoSinglePhotonEtFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1IsoSinglePhotonEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedPhotonHcalIsol+process.hltL1IsoSinglePhotonHcalIsolFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1IsoSinglePhotonTrackIsolFilter+process.HLTEndSequence))
process.HLT_ZeroBias = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sZeroBias+process.hltPreZeroBias+process.HLTEndSequence))
process.HLT_BTagMu_HT370 = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sBSoftMuonHT+process.hltPreBTagMuHT370+process.HLTBCommonL2recoSequence+process.hltBSoftmuonHTL2filter+process.HLTBSoftmuonL25recoSequence+process.hltBSoftmuonL25filter+process.HLTBSoftmuonL3recoSequence+process.hltBSoftmuonL3filter+process.HLTEndSequence))
process.HLT_DoubleIsoPhoton20_L1I = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sDoubleEgamma+process.hltPreDoubleIsoPhoton20L1I+process.HLTDoRegionalEgammaEcalSequence+process.HLTL1IsolatedEcalClustersSequence+process.hltL1IsoRecoEcalCandidate+process.hltL1IsoDoublePhotonL1MatchFilterRegional+process.hltL1IsoDoublePhotonEtFilter+process.hltL1IsolatedPhotonEcalIsol+process.hltL1IsoDoublePhotonEcalIsolFilter+process.HLTDoLocalHcalWithoutHOSequence+process.hltL1IsolatedPhotonHcalIsol+process.hltL1IsoDoublePhotonHcalIsolFilter+process.HLTDoLocalTrackerSequence+process.HLTL1IsoEgammaRegionalRecoTrackerSequence+process.hltL1IsoPhotonTrackIsol+process.hltL1IsoDoublePhotonTrackIsolFilter+process.hltL1IsoDoublePhotonDoubleEtFilter+process.HLTEndSequence))
process.HLT_IsoPhoton15_L1R = cms.Path(process.ProductionFilterSequence*(process.HLTBeginSequence+process.hltL1sRelaxedSingleEgammaEt12+process.hltPreIsoPhoton15L1R+process.HLTSinglePhoton15L1NonIsolatedHLTIsoSequence+process.HLTEndSequence))
process.endjob_step = cms.Path(process.ProductionFilterSequence*process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.HLTriggerFirstPath,process.HLT_L1Jet15,process.HLT_Jet30,process.HLT_Jet50,process.HLT_Jet80,process.HLT_Jet110,process.HLT_Jet180,process.HLT_Jet250,process.HLT_FwdJet20,process.HLT_DoubleJet150,process.HLT_DoubleJet125_Aco,process.HLT_DoubleFwdJet50,process.HLT_DiJetAve15,process.HLT_DiJetAve30,process.HLT_DiJetAve50,process.HLT_DiJetAve70,process.HLT_DiJetAve130,process.HLT_DiJetAve220,process.HLT_TripleJet85,process.HLT_QuadJet30,process.HLT_QuadJet60,process.HLT_SumET120,process.HLT_L1MET20,process.HLT_MET25,process.HLT_MET35,process.HLT_MET50,process.HLT_MET65,process.HLT_MET75,process.HLT_MET65_HT350,process.HLT_Jet180_MET60,process.HLT_Jet60_MET70_Aco,process.HLT_Jet100_MET60_Aco,process.HLT_DoubleJet125_MET60,process.HLT_DoubleFwdJet40_MET60,process.HLT_DoubleJet60_MET60_Aco,process.HLT_DoubleJet50_MET70_Aco,process.HLT_DoubleJet40_MET70_Aco,process.HLT_TripleJet60_MET60,process.HLT_QuadJet35_MET60,process.HLT_IsoEle15_L1I,process.HLT_IsoEle18_L1R,process.HLT_IsoEle15_LW_L1I,process.HLT_LooseIsoEle15_LW_L1R,process.HLT_Ele10_SW_L1R,process.HLT_Ele15_SW_L1R,process.HLT_Ele15_LW_L1R,process.HLT_EM80,process.HLT_EM200,process.HLT_DoubleIsoEle10_L1I,process.HLT_DoubleIsoEle12_L1R,process.HLT_DoubleIsoEle10_LW_L1I,process.HLT_DoubleIsoEle12_LW_L1R,process.HLT_DoubleEle5_SW_L1R,process.HLT_DoubleEle10_LW_OnlyPixelM_L1R,process.HLT_DoubleEle10_Z,process.HLT_DoubleEle6_Exclusive,process.HLT_IsoPhoton30_L1I,process.HLT_IsoPhoton10_L1R,process.HLT_IsoPhoton15_L1R,process.HLT_IsoPhoton20_L1R,process.HLT_IsoPhoton25_L1R,process.HLT_IsoPhoton40_L1R,process.HLT_Photon15_L1R,process.HLT_Photon25_L1R,process.HLT_DoubleIsoPhoton20_L1I,process.HLT_DoubleIsoPhoton20_L1R,process.HLT_DoublePhoton10_Exclusive,process.HLT_L1Mu,process.HLT_L1MuOpen,process.HLT_L2Mu9,process.HLT_IsoMu9,process.HLT_IsoMu11,process.HLT_IsoMu13,process.HLT_IsoMu15,process.HLT_Mu3,process.HLT_Mu5,process.HLT_Mu7,process.HLT_Mu9,process.HLT_Mu11,process.HLT_Mu13,process.HLT_Mu15,process.HLT_Mu15_Vtx2mm,process.HLT_DoubleIsoMu3,process.HLT_DoubleMu3,process.HLT_DoubleMu3_Vtx2mm,process.HLT_DoubleMu3_JPsi,process.HLT_DoubleMu3_Upsilon,process.HLT_DoubleMu7_Z,process.HLT_DoubleMu3_SameSign,process.HLT_DoubleMu3_Psi2S,process.HLT_BTagIP_Jet180,process.HLT_BTagIP_Jet120_Relaxed,process.HLT_BTagIP_DoubleJet120,process.HLT_BTagIP_DoubleJet60_Relaxed,process.HLT_BTagIP_TripleJet70,process.HLT_BTagIP_TripleJet40_Relaxed,process.HLT_BTagIP_QuadJet40,process.HLT_BTagIP_QuadJet30_Relaxed,process.HLT_BTagIP_HT470,process.HLT_BTagIP_HT320_Relaxed,process.HLT_BTagMu_DoubleJet120,process.HLT_BTagMu_DoubleJet60_Relaxed,process.HLT_BTagMu_TripleJet70,process.HLT_BTagMu_TripleJet40_Relaxed,process.HLT_BTagMu_QuadJet40,process.HLT_BTagMu_QuadJet30_Relaxed,process.HLT_BTagMu_HT370,process.HLT_BTagMu_HT250_Relaxed,process.HLT_DoubleMu3_BJPsi,process.HLT_DoubleMu4_BJPsi,process.HLT_TripleMu3_TauTo3Mu,process.HLT_IsoTau_MET65_Trk20,process.HLT_IsoTau_MET35_Trk15_L1MET,process.HLT_LooseIsoTau_MET30,process.HLT_LooseIsoTau_MET30_L1MET,process.HLT_DoubleIsoTau_Trk3,process.HLT_DoubleLooseIsoTau,process.HLT_IsoEle8_IsoMu7,process.HLT_IsoEle10_Mu10_L1R,process.HLT_IsoEle12_IsoTau_Trk3,process.HLT_IsoEle10_BTagIP_Jet35,process.HLT_IsoEle12_Jet40,process.HLT_IsoEle12_DoubleJet80,process.HLT_IsoEle5_TripleJet30,process.HLT_IsoEle12_TripleJet60,process.HLT_IsoEle12_QuadJet35,process.HLT_IsoMu14_IsoTau_Trk3,process.HLT_IsoMu7_BTagIP_Jet35,process.HLT_IsoMu7_BTagMu_Jet20,process.HLT_IsoMu7_Jet40,process.HLT_NoL2IsoMu8_Jet40,process.HLT_Mu14_Jet50,process.HLT_Mu5_TripleJet30,process.HLT_BTagMu_Jet20_Calib,process.HLT_ZeroBias,process.HLT_MinBias,process.HLT_MinBiasHcal,process.HLT_MinBiasEcal,process.HLT_MinBiasPixel,process.HLT_MinBiasPixel_Trk5,process.HLT_BackwardBSC,process.HLT_ForwardBSC,process.HLT_CSCBeamHalo,process.HLT_CSCBeamHaloOverlapRing1,process.HLT_CSCBeamHaloOverlapRing2,process.HLT_CSCBeamHaloRing2or3,process.HLT_TrackerCosmics,process.AlCa_IsoTrack,process.AlCa_EcalPhiSym,process.AlCa_HcalPhiSym,process.AlCa_EcalPi0,process.HLTriggerFinalPath,process.endjob_step,process.out_step)
