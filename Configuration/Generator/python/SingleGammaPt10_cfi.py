import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.untracked.PSet(
        MaxPt = cms.untracked.double(10.01),
        MinPt = cms.untracked.double(9.99),
        PartID = cms.untracked.vint32(22),
        MaxEta = cms.untracked.double(2.5),
        MaxPhi = cms.untracked.double(3.14159265359),
        MinEta = cms.untracked.double(-2.5),
        MinPhi = cms.untracked.double(-3.14159265359) ## in radians

    ),
    Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

    psethack = cms.string('single gamma pt 10'),
    firstRun = cms.untracked.uint32(1)
)



ProductionFilterSequence = cms.Sequence(generator)
