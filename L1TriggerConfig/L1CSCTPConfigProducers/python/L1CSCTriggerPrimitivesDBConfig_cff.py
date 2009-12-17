import FWCore.ParameterSet.Config as cms

#from L1TriggerConfig.L1CSCTPConfigProducers.L1CSCTriggerPrimitivesConfig_cfi import *

from CondCore.DBCommon.CondDBSetup_cfi import *

# Read constants from DB.
l1csctpdbconfsrc = cms.ESSource("PoolDBESSource",
    CondDBSetup,
    timetype = cms.string('runnumber'),
    connect = cms.string('frontier://FrontierDev/CMS_COND_CSC'),
    #connect = cms.string('frontier://cmsfrontier.cern.ch:8000/FrontierProd/CMS_COND_21X_CSC'),
    authenticationMethod = cms.untracked.uint32(1),
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('CSCL1TPParametersRcd'),
            tag = cms.string('CSCL1TPParameters')
            #type = cms.string('L1CSCTPParameters')
        )
    )
 
)

l1csctpdbconfsrc.DBParameters.authenticationPath = cms.untracked.string('/afs/cern.ch/cms/DB/conddb')


# Reading from DB has precedence over dummy producers (which use constants
# defined in cfi files).
es_prefer_l1csctpdbconfsrc = cms.ESPrefer("PoolDBESSource","l1csctpdbconfsrc")
