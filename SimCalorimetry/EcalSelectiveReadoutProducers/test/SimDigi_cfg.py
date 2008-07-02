import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# event vertex smearing - applies only once (internal check)
# Note : all internal generators will always do (0,0,0) vertex
#
process.load("IOMC.EventVertexGenerators.VtxSmearedGauss_cfi")

#Geometry
#
#include "Geometry/CMSCommonData/data/cmsSimIdealGeometryXML.cfi"
process.load("Geometry.EcalCommonData.EcalOnly_cfi")

# Calo geometry service model
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")

process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")

process.load("Geometry.EcalMapping.EcalMapping_cfi")

process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

#Magnetic Field
#
process.load("Configuration.StandardSequences.MagneticField_cff")

# Step 2 : CMS Detector Simulation
process.load("SimG4Core.Application.g4SimHits_cfi")

# Pile-up processing:
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

#
# Ecal digi production:
process.load("SimCalorimetry.EcalSimProducers.ecaldigi_cfi")

#
#
# Get hardcoded conditions the same used for standard digitization
process.load("CalibCalorimetry.Configuration.Ecal_FakeConditions_cff")

# TPG: defines the ecalTriggerPrimitiveDigis module
process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTriggerPrimitiveDigis_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.o1 = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(9),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_simEcalUnsuppressedDigis_*_*', 
        'keep *_simEcalUnsuppressedDigis_*_*', 
        'keep *_simEcalTriggerPrimitiveDigis_*_*'),
    fileName = cms.untracked.string('srp_validation_in.root')
)

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('cout')
)

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    moduleSeeds = cms.PSet(
        g4SimHits = cms.untracked.uint32(123456782),
        simEcalUnsuppressedDigis = cms.untracked.uint32(123456784),
        ecalUnsuppressedDigis = cms.untracked.uint32(123456783),
        VtxSmeared = cms.untracked.uint32(123456781)
    ),
    sourceSeed = cms.untracked.uint32(135799753)
)

process.source = cms.Source("FlatRandomEGunSource",
    PGunParameters = cms.untracked.PSet(
        # you can request more than 1 particle
        #vint32  PartID = {211,11}
        #vint32 PartID = { 13 } 
        PartID = cms.untracked.vint32(11),
        MaxEta = cms.untracked.double(3.0),
        MaxPhi = cms.untracked.double(3.141592653589793),
        MinEta = cms.untracked.double(-3.0),
        MinE = cms.untracked.double(99.99),
        MinPhi = cms.untracked.double(-3.141592653589793), ## must be in radians

        MaxE = cms.untracked.double(100.01)
    ),
    Verbosity = cms.untracked.int32(0) ## set to 1 (or greater)  for printouts

)

process.detSim = cms.Sequence(process.VtxSmeared*process.g4SimHits)
process.p = cms.Path(process.detSim*process.mix*process.simEcalUnsuppressedDigis*process.simEcalTriggerPrimitiveDigis)
process.outpath = cms.EndPath(process.o1)

