import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

process.load("IOMC.EventVertexGenerators.VtxSmearedBeamProfile_cfi")

process.load("SimG4CMS.HcalTestBeam.test.TB2004GeometryXML_cfi")

process.load("Configuration.EventContent.EventContent_cff")

process.load("SimG4Core.Application.g4SimHits_cfi")

process.load("DQMServices.Core.DQM_cfg")

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalGeom = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HcalTBSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        CaloSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        DEBUG = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        VertexGenerator = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HCalGeom = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HcalSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        threshold = cms.untracked.string('INFO'),
        SimHCalData = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        EcalSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    categories = cms.untracked.vstring('CaloSim', 
        'EcalGeom', 
        'EcalSim', 
        'HCalGeom', 
        'HcalSim', 
        'HcalTBSim', 
        'SimHCalData', 
        'VertexGenerator'),
    destinations = cms.untracked.vstring('cout')
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(9876),
        VtxSmeared = cms.untracked.uint32(123456789)
    ),
    sourceSeed = cms.untracked.uint32(135799753)
)

process.common_beam_direction_parameters = cms.PSet(
    MaxEta = cms.untracked.double(0.5655),
    MaxPhi = cms.untracked.double(-0.1309),
    MinEta = cms.untracked.double(0.5655),
    MinPhi = cms.untracked.double(-0.1309),
    BeamPosition = cms.untracked.double(-521.5)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("FlatRandomEGunSource",
    PGunParameters = cms.untracked.PSet(
        process.common_beam_direction_parameters,
        MaxE = cms.untracked.double(10.01),
        MinE = cms.untracked.double(9.99),
        PartID = cms.untracked.vint32(211)
    ),
    Verbosity = cms.untracked.int32(0)
)

process.o1 = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    fileName = cms.untracked.string('sim2004.root')
)

process.Timing = cms.Service("Timing")

process.p1 = cms.Path(process.VtxSmeared*process.g4SimHits)
process.outpath = cms.EndPath(process.o1)
process.g4SimHits.UseMagneticField = False
process.g4SimHits.Physics.type = 'SimG4Core/Physics/QGSP'
process.g4SimHits.CaloSD = cms.PSet(
    process.common_beam_direction_parameters,
    process.common_heavy_suppression,
    SuppressHeavy = cms.bool(False),
    DetailedTiming = cms.untracked.bool(False),
    Verbosity = cms.untracked.int32(0),
    CheckHits = cms.untracked.int32(25),
    CorrectTOFBeam = cms.untracked.bool(False),
    UseMap = cms.untracked.bool(True),
    EminTrack = cms.double(1.0)
)
process.g4SimHits.ECalSD.UseBirkLaw = False
process.g4SimHits.ECalSD.BirkC1 = 0.33333
process.g4SimHits.ECalSD.BirkC2 = 0.0
process.g4SimHits.ECalSD.BirkL3Parametrization = True
process.g4SimHits.ECalSD.SlopeLightYield = 0.05
process.g4SimHits.HCalSD.UseBirkLaw = False
process.g4SimHits.HCalSD.BirkC1 = 0.013
process.g4SimHits.HCalSD.BirkC2 = '9.6e-6'
process.g4SimHits.HCalSD.UseShowerLibrary = False
process.g4SimHits.HCalSD.TestNumberingScheme = True
process.g4SimHits.HCalSD.UseHF = False
process.g4SimHits.HCalSD.ForTBH2 = True
process.g4SimHits.CaloTrkProcessing.TestBeam = True
process.g4SimHits.Watchers = cms.VPSet(cms.PSet(
    HcalTB04Analysis = cms.PSet(
        process.common_beam_direction_parameters,
        ScaleHB16 = cms.double(0.5),
        HcalOnly = cms.bool(False),
        ScaleHO = cms.double(0.4),
        Type = cms.int32(2),
        Verbose = cms.untracked.bool(True),
        FileName = cms.untracked.string('HcalTB04.root'),
        Names = cms.vstring('HcalHits', 
            'EcalHitsEB'),
        ETtotMax = cms.untracked.double(20.0),
        ScaleHE0 = cms.double(0.5),
        Mode = cms.int32(1),
        ScaleHB0 = cms.double(0.5),
        EcalNoise = cms.double(0.13),
        EHCalMax = cms.untracked.double(2.0)
    ),
    HcalQie = cms.PSet(
        NumOfBuckets = cms.int32(10),
        BaseLine = cms.int32(4),
        BinOfMax = cms.int32(6),
        PreSamples = cms.int32(0),
        EDepPerPE = cms.double(0.0005),
        SignalBuckets = cms.int32(2),
        SigmaNoise = cms.double(0.5),
        qToPE = cms.double(4.0)
    ),
    type = cms.string('HcalTB04Analysis')
))
process.DQM.collectorHost = ''

