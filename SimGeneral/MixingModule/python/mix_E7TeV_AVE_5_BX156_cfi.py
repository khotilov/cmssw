import FWCore.ParameterSet.Config as cms

# configuration to model pileup for initial physics phase
from SimGeneral.MixingModule.mixObjects_cfi import *
mix = cms.EDProducer("MixingModule",
    LabelPlayback = cms.string(''),
    maxBunch = cms.int32(3),
    minBunch = cms.int32(-5), ## in terms of 25 nsec

    bunchspace = cms.int32(450), ##ns
    mixProdStep1 = cms.bool(False),
    mixProdStep2 = cms.bool(False),

    playback = cms.untracked.bool(False),
    useCurrentProcessOnly = cms.bool(False),
                   
    input = cms.SecSource("PoolSource",
        nbPileupEvents = cms.PSet(
            averageNumber = cms.double(5.0)
        ),
        seed = cms.int32(1234567),
        type = cms.string('poisson'),
	sequential = cms.untracked.bool(False),
        fileNames = cms.untracked.vstring('/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/18890F4C-FD99-DD11-BFF9-000423D996C8.root',
        '/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/1A423C16-6099-DD11-9320-000423D9853C.root',
        '/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/5EDA8A7F-5D99-DD11-B1CD-001617C3B706.root',
        '/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/600D1E6A-5F99-DD11-A7D5-000423D9890C.root',
        '/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/68ECAE92-5F99-DD11-ACAB-000423D98E6C.root',
        '/store/relval/CMSSW_2_1_10/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP_V7_v2/0000/8802D325-5E99-DD11-B858-000423D98A44.root')
    ),
    mixObjects = cms.PSet(
        mixCH = cms.PSet(
            mixCaloHits
        ),
        mixTracks = cms.PSet(
            mixSimTracks
        ),
        mixVertices = cms.PSet(
            mixSimVertices
        ),
        mixSH = cms.PSet(
            mixSimHits
        ),
        mixHepMC = cms.PSet(
            mixHepMCProducts
        )
    )
)


