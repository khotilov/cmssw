import FWCore.ParameterSet.Config as cms


process = cms.Process("readDB")
process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.load("CondCore.DBCommon.CondDBSetup_cfi")



process.BeamSpotDBSource = cms.ESSource("PoolDBESSource",
                                        process.CondDBSetup,
                                        toGet = cms.VPSet(cms.PSet(
    record = cms.string('BeamSpotObjectsRcd'),
    tag = cms.string('Early7TeVCollision_4p2cm_v1_mc_IDEAL')
    )),
                                        connect = cms.string('sqlite_file:EarlyCollision_IDEAL.db')
                                        #connect = cms.string('oracle://cms_orcon_prod/CMS_COND_21X_BEAMSPOT')
                                        )

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(1)
                    )
process.beamspot = cms.EDFilter("BeamSpotFromDB")


process.p = cms.Path(process.beamspot)

