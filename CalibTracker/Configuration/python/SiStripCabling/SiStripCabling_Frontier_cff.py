import FWCore.ParameterSet.Config as cms

import copy
from CalibTracker.Configuration.Common.PoolDBESSource_cfi import *
siStripFedCabling = copy.deepcopy(poolDBESSource)
sistripconn = cms.ESProducer("SiStripConnectivity")

siStripFedCabling.connect = 'frontier://FrontierProd/CMS_COND_20X_STRIP'
siStripFedCabling.toGet = cms.VPSet(cms.PSet(
    record = cms.string('SiStripFedCablingRcd'),
    tag = cms.string('SiStripFedCabling_20X')
))

