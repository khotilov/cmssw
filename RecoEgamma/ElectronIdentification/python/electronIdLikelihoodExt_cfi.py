import FWCore.ParameterSet.Config as cms

from RecoEcal.EgammaClusterProducers.geometryForClustering_cff import *
from RecoEgamma.ElectronIdentification.likelihoodPdfsDB_cfi import *
from RecoEgamma.ElectronIdentification.likelihoodESetup_cfi import *
eidLikelihoodExt = cms.EDProducer("EleIdLikelihoodExtProducer",
    src = cms.InputTag("pixelMatchGsfElectrons"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    doLikelihood = cms.bool(True),
    filter = cms.bool(False),
    threshold = cms.double(0.5),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB")
)


