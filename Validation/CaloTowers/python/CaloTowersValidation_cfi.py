import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff import *
from Configuration.StandardSequences.Simulation_cff import *
from Configuration.StandardSequences.MixingNoPileUp_cff import *
from Configuration.StandardSequences.Reconstruction_cff import *
from Configuration.StandardSequences.GeometryECALHCAL_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
GlobalTag.globaltag = 'IDEAL_30X::All'

from DQMServices.Core.DQM_cfg import *
maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(True),
    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring('file:/')
)

MessageLogger = cms.Service("MessageLogger")

myanalyzer = cms.EDFilter("CaloTowersValidation",
    outputFile = cms.untracked.string('CaloTowersValidationHB.root'),
    CaloTowerCollectionLabel = cms.untracked.string('towerMaker'),
    hcalselector = cms.untracked.string('HB')
)

DQM.collectorHost = ''

hbhereco.digiLabel = 'simHcalDigis'
horeco.digiLabel = 'simHcalDigis'
hfreco.digiLabel = 'simHcalDigis'

ecalPreshowerRecHit.ESdigiCollection = 'simEcalPreshowerDigis'
ecalWeightUncalibRecHit.EBdigiCollection = 'simEcalDigis:ebDigis'
ecalWeightUncalibRecHit.EEdigiCollection = 'simEcalDigis:eeDigis'

p = cms.Path(mix*calDigi*calolocalreco*caloTowersRec*myanalyzer)

