import FWCore.ParameterSet.Config as cms

from CommonTools.SiStripClusterization.SiStripClusterization_cfi import *
SiStripRawToClustersFacility = cms.EDProducer("SiStripRawToClusters",
    SiStripClusterization,
    FedRawData = cms.InputTag('rawDataCollector')
)


