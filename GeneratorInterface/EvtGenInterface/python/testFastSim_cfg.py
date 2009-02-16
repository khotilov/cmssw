import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

from Configuration.Generator.PythiaUESettings_cfi import *
process.load('Configuration/StandardSequences/MagneticField_38T_cff')

process.load('Configuration/StandardSequences/Generator_cff')

process.load('FastSimulation/Configuration/FamosSequences_cff')

process.load('Configuration.StandardSequences.L1TriggerDefaultMenu_cff')

process.load('FastSimulation.PileUpProducer.PileUpSimulator10TeV_cfi')

process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')

process.load('FastSimulation/Configuration/CommonInputs_cff')

process.load('FastSimulation/Configuration/EventContent_cff')

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
       threshold = cms.untracked.string('ERROR')
    ),
    destinations = cms.untracked.vstring('pythiaevtgen.log', 
        'cout')
)
process.Timing = cms.Service("Timing")

process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
     saveFileName = cms.untracked.string(''), 
     theSource = cms.PSet(
                    initialSeed = cms.untracked.uint32(123456789),
                    engineName = cms.untracked.string('TRandom3')
     ),
     VtxSmeared = cms.PSet(
                    initialSeed = cms.untracked.uint32(123456789),
                    engineName = cms.untracked.string('TRandom3')
     ),
     famosPileUp = cms.PSet(
                    initialSeed = cms.untracked.uint32(918273),
                    engineName = cms.untracked.string('TRandom3')
     ),
     famosSimHits = cms.PSet(
                    initialSeed = cms.untracked.uint32(13579),
                    engineName = cms.untracked.string('TRandom3')
     ),
     siTrackerGaussianSmearingRecHits = cms.PSet(
                    initialSeed = cms.untracked.uint32(24680),
                    engineName = cms.untracked.string('TRandom3')
     ),
     caloRecHits = cms.PSet(
                    initialSeed = cms.untracked.uint32(654321),
                    engineName = cms.untracked.string('TRandom3')
     ),
     paramMuons = cms.PSet(
                    initialSeed = cms.untracked.uint32(54525),
                    engineName = cms.untracked.string('TRandom3')
     ),
     l1ParamMuons = cms.PSet(
                    initialSeed = cms.untracked.uint32(6453209),
                    engineName = cms.untracked.string('TRandom3')
     ),
     MuonSimHits = cms.PSet(
                    initialSeed = cms.untracked.uint32(987346),
                    engineName = cms.untracked.string('TRandom3')
     ),
     simMuonRPCDigis = cms.PSet(
                    initialSeed = cms.untracked.uint32(524964),
                    engineName = cms.untracked.string('TRandom3')
     ),
     simMuonCSCDigis = cms.PSet(
                    initialSeed = cms.untracked.uint32(525432),
                    engineName = cms.untracked.string('TRandom3')
     ),
     simMuonDTDigis = cms.PSet(
                    initialSeed = cms.untracked.uint32(67673876),
                    engineName = cms.untracked.string('TRandom3')
     ),
     evtgenproducer = cms.PSet(
                   initialSeed = cms.untracked.uint32(67676)
     )
   
)

process.load("Configuration.Generator.PythiaUESettings_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)
process.source = cms.Source("PythiaSource",
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    maxEventsToPrint = cms.untracked.int32(4),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        process.pythiaUESettingsBlock,
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
            'MSEL=5         ! b-bbar'),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

process.evtgenproducer = cms.EDProducer("EvtGenProducer",
     use_default_decay = cms.untracked.bool(True),
     decay_table = cms.FileInPath('GeneratorInterface/EvtGenInterface//data/DECAY.DEC'),
     particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
     user_decay_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/Bd_Kstarmumu_Kpi.dec'),
     list_forced_decays = cms.vstring(),
     processParameters = cms.vstring('MSEL=5     ! bbbar', 
            'MDCY(134,1) = 0', 
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
            'MDCY(105,1) = 0')
 )

process.myout = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(

            'drop *',
            "keep *_siTrackerGaussianSmearingRecHits_*_*",
	    "keep *_evtgenproducer_*_*",
	    "keep SimTracks_*_*_*",
	    "keep *_*_TrackerHits_*",
	    "keep *_generalTracks_*_*",
	    "drop *_*_MuonSimTracks_*",
	    "drop TrajectorysToOnerecoTracksAssociation_*_*_*",
            'drop *_source_*_*'
         )
)

process.p1 = cms.Path(process.evtgenproducer*process.famosWithTracks*process.randomEngineStateProducer)

process.outpath = cms.EndPath(process.myout)
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.famosSimHits.SourceLabel = "evtgenproducer"
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True
process.famosPileUp.PileUpSimulator.averageNumber = 5.0

