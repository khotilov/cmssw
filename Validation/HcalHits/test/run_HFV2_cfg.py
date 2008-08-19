import FWCore.ParameterSet.Config as cms

process = cms.Process("CaloTest")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")

#Magnetic Field 		
process.load("Configuration.StandardSequences.MagneticField_cff")

#Geometry
process.load("Validation.HcalHits.testGeometryPMTXML_cfi")

# Calo geometry service model
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")

process.load("SimG4Core.Application.g4SimHits_cfi")

process.load("Validation.HcalHits.HcalHitValidation_cfi")

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout'),
    categories = cms.untracked.vstring('HFShower', 
        'HcalSim', 
        'CaloSim', 
        'G4cout', 
        'G4cerr', 
        'ValidHcal'),
    cout = cms.untracked.PSet(
        G4cerr = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        ValidHcal = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        CaloSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HFShower = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        G4cout = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        HcalSim = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    )
)

process.Timing = cms.Service("Timing")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(9876),
        VtxSmeared = cms.untracked.uint32(123456789)
    ),
    sourceSeed = cms.untracked.uint32(135799753)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("FlatRandomEGunSource",
    PGunParameters = cms.untracked.PSet(
        # you can request more than 1 particle
        PartID = cms.untracked.vint32(11),
        MinEta = cms.untracked.double(3.070),
        MaxEta = cms.untracked.double(3.071),
        MinPhi = cms.untracked.double(0.0872),
        MaxPhi = cms.untracked.double(0.0873),
        #untracked double MinEta = 3.0
        #untracked double MaxEta = 4.5
        MinE = cms.untracked.double(99.9),
        MaxE = cms.untracked.double(100.1)
    ),
    firstRun = cms.untracked.uint32(1)
)

process.USER = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('100_pi_bs.root')
)

process.DQMStore = cms.Service("DQMStore",
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0)
)

process.DQM = cms.Service("DQM",
    debug = cms.untracked.bool(False),
    publishFrequency = cms.untracked.double(5.0),
    collectorPort = cms.untracked.int32(9090),
    collectorHost = cms.untracked.string('')
)

process.p1 = cms.Path(process.VtxSmeared*process.g4SimHits*process.hcalHitValid)
process.outpath = cms.EndPath(process.USER)
process.VtxSmeared.SigmaX = 0.00001
process.VtxSmeared.SigmaY = 0.00001
process.VtxSmeared.SigmaZ = 0.00001
process.g4SimHits.UseMagneticField = False
process.g4SimHits.Physics.type = 'SimG4Core/Physics/QGSP_EMV'
process.g4SimHits.HCalSD.UseParametrize = True
process.g4SimHits.HCalSD.UsePMTHits = True
process.g4SimHits.HCalSD.BetaThreshold = 0.70
process.g4SimHits.HFShower.TrackEM = True
process.g4SimHits.G4Commands = ['/physics_engine/neutron/energyLimit 0 keV', '/physics_engine/neutron/timeLimit 0.01 ms']
process.g4SimHits.Watchers = cms.VPSet(cms.PSet(
    SimG4HcalValidation = cms.PSet(
        TimeLowLimit = cms.double(0.0),
        LabelNxNInfo = cms.untracked.string('HcalInfoNxN'),
        LabelLayerInfo = cms.untracked.string('HcalInfoLayer'),
        HcalHitThreshold = cms.double(1e-20),
        Phi0 = cms.double(0.3054),
        ConeSize = cms.double(0.5),
        InfoLevel = cms.int32(2),
        JetThreshold = cms.double(5.0),
        EcalHitThreshold = cms.double(1e-20),
        TimeUpLimit = cms.double(999.0),
        HcalClusterOnly = cms.bool(False),
        Eta0 = cms.double(0.3045),
        LabelJetsInfo = cms.untracked.string('HcalInfoJets'),
        Names = cms.vstring('HcalHits', 
            'EcalHitsEB', 
            'EcalHitsEE', 
            'EcalHitsES'),
        HcalSampling = cms.bool(True)
    ),
    type = cms.string('SimG4HcalValidation')
))
process.hcalHitValid.outputFile = '100_pi_bs_plots.root'


