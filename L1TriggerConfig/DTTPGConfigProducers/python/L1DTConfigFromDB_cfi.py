import FWCore.ParameterSet.Config as cms

#Include configuration ParameterSets
from L1TriggerConfig.DTTPGConfigProducers.L1DTConfigParams_cff import *
from L1TriggerConfig.DTTPGConfigProducers.L1DTConfigMap_cff import *
L1DTConfigFromDB = cms.ESProducer("DTConfigDBProducer",
    DTTPGMapBlock,
    DTTPGParametersBlock,
    cfgConfig  = cms.bool(False),
    TracoLutsFromDB = cms.bool(False),
    UseBtiAcceptParam = cms.bool(False),
    debugDB    = cms.bool(False),
    debugBti   = cms.int32(0),
    debugTraco = cms.int32(0),
    debugTSP   = cms.bool(False),
    debugTST   = cms.bool(False),
    debugTU    = cms.bool(False),
    debugSC    = cms.bool(False),
    debugLUTs  = cms.bool(False),             
    debug      = cms.bool(False),
    catalog = cms.string("file:testcatalog.xml"),
    authPath = cms.string("/afs/cern.ch/cms/DB/conddb"),
    siteLocalConfig = cms.bool(False),
##### CHANGED  for 3XY #####
###    contact = cms.string("oracle://cms_orcoff_prod/CMS_COND_21X_DT"),
    contact = cms.string("oracle://cms_orcoff_prod/CMS_COND_31X_DT"),
#    connect = cms.string("frontier://FrontierProd/CMS_COND_31X_DT"),
##### CHANGED  for 3XY #####
###    tag = cms.string("conf_ccb_V01"),
    tag = cms.string("DT_config_V02"),
##### CHANGED  for 3XY #####
###    token = cms.string("[DB=00000000-0000-0000-0000-000000000000][CNT=DTConfigList][CLID=9CB14BE8-30A2-DB11-9935-000E0C5CE283][TECH=00000B01][OID=00000004-00000000]")
    token = cms.string("DB=00000000-0000-0000-0000-000000000000][CNT=DTConfigList][CLID=9CB14BE8-30A2-DB11-9935-000E0C5CE283][TECH=00000B01][OID=0000000B-00000000]")
)



