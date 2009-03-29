def initCondDBSource( process,
                      inputDBConnect = 'frontier://FrontierPrep/CMS_COND_L1T',
                      inputDBAuth = '.',
                      tagBase = 'CRAFT_hlt',
                      use30XTagList = False ):
    import FWCore.ParameterSet.Config as cms
    from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup

    process.l1conddb = cms.ESSource("PoolDBESSource",
                            CondDBSetup,
                            connect = cms.string(inputDBConnect),
                            toGet = cms.VPSet(cms.PSet(
        record = cms.string('L1TriggerKeyListRcd'),
        tag = cms.string('L1TriggerKeyList_' + tagBase)
        ),
                                              cms.PSet(
        record = cms.string('L1TriggerKeyRcd'),
        tag = cms.string('L1TriggerKey_' + tagBase)
        )),
                            BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService')
                            )
    process.l1conddb.DBParameters.authenticationPath = inputDBAuth

    if use30XTagList == True:
        from CondTools.L1Trigger.L1SubsystemParams30X_cfi import initL1Subsystems
    else:
        from CondTools.L1Trigger.L1SubsystemParams_cfi import initL1Subsystems
    initL1Subsystems( tagBase = tagBase )
    process.l1conddb.toGet.extend(initL1Subsystems.params.recordInfo)

    process.es_prefer_l1conddb = cms.ESPrefer("PoolDBESSource","l1conddb")
