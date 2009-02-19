import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.untracked.PSet(
        PartID = cms.untracked.vint32(111),
        MaxEta = cms.untracked.double(3.0),
        MaxPhi = cms.untracked.double(3.14159265359),
        MinEta = cms.untracked.double(-3.0),
        MinE = cms.untracked.double(9.99),
        MinPhi = cms.untracked.double(-3.14159265359), ## in radians

        MaxE = cms.untracked.double(10.01)
    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('single pi0 E 10'),
    firstRun = cms.untracked.uint32(1)
)



ProductionFilterSequence = cms.Sequence(generator)
