import FWCore.ParameterSet.Config as cms

hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)

from Geometry.HcalEventSetup.HcalRelabel_cfi import HcalReLabel
#------------------------- HARDCODED conditions

es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
toGet = cms.untracked.vstring('GainWidths'),
    iLumi = cms.double(20.),                 # for Upgrade
    HcalReLabel =  HcalReLabel,              # for Upgrade
    HERecalibration = cms.bool(False),       # for Upgrade   
    HFRecalibration = cms.bool(False)        # for Upgrade   
)


es_prefer_hcalHardcode = cms.ESPrefer("HcalHardcodeCalibrations", "es_hardcode")

import Geometry.HcalEventSetup.hcalTopologyConstants_cfi as hcalTopologyConstants_cfi
es_hardcode.hcalTopologyConstants = cms.PSet(hcalTopologyConstants_cfi.hcalTopologyConstants)
