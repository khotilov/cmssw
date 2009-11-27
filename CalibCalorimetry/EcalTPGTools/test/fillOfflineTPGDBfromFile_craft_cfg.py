import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('*'),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('DEBUG')
    ),
    destinations = cms.untracked.vstring('cout')
)

process.source = cms.Source("EmptyIOVSource",
    lastValue = cms.uint64(1),
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    interval = cms.uint64(1)
)

process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_craft_cff")

process.load("CondCore.DBCommon.CondDBCommon_cfi")

process.CondDBCommon.connect = 'sqlite_file:DB_craft.db'

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDBCommon,
    messagelevel = cms.untracked.uint32(3),
    catalog = cms.untracked.string('file:PoolFileCatalog_DB.xml'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('EcalTPGPedestalsRcd'),
        tag = cms.string('EcalTPGPedestals_craft')
    ), 
        cms.PSet(
            record = cms.string('EcalTPGLinearizationConstRcd'),
            tag = cms.string('EcalTPGLinearizationConst_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGSlidingWindowRcd'),
            tag = cms.string('EcalTPGSlidingWindow_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainEBIdMapRcd'),
            tag = cms.string('EcalTPGFineGrainEBIdMap_craft')
        ),
        cms.PSet(
            record = cms.string('EcalTPGFineGrainStripEERcd'),
            tag = cms.string('EcalTPGFineGrainStripEE_craft')
        ),
        cms.PSet(
            record = cms.string('EcalTPGFineGrainTowerEERcd'),
            tag = cms.string('EcalTPGFineGrainTowerEE_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGLutIdMapRcd'),
            tag = cms.string('EcalTPGLutIdMap_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGWeightIdMapRcd'),
            tag = cms.string('EcalTPGWeightIdMap_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGWeightGroupRcd'),
            tag = cms.string('EcalTPGWeightGroup_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGLutGroupRcd'),
            tag = cms.string('EcalTPGLutGroup_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainEBGroupRcd'),
            tag = cms.string('EcalTPGFineGrainEBGroup_craft')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGPhysicsConstRcd'),
            tag = cms.string('EcalTPGPhysicsConst_craft')
        ),
	cms.PSet(
            record = cms.string('EcalTPGCrystalStatusRcd'),
            tag = cms.string('EcalTPGCrystalStatus_craft')
        ),
	cms.PSet(
            record = cms.string('EcalTPGTowerStatusRcd'),
            tag = cms.string('EcalTPGTowerStatus_craft')
        ))
)

process.dbCopy = cms.EDFilter("EcalTPGDBCopy",
    timetype = cms.string('runnumber'),
    toCopy = cms.VPSet(cms.PSet(
        record = cms.string('EcalTPGPedestalsRcd'),
        container = cms.string('EcalTPGPedestals')
    ), 
        cms.PSet(
            record = cms.string('EcalTPGLinearizationConstRcd'),
            container = cms.string('EcalTPGLinearizationConst')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGSlidingWindowRcd'),
            container = cms.string('EcalTPGSlidingWindow')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainEBIdMapRcd'),
            container = cms.string('EcalTPGFineGrainEBIdMap')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainStripEERcd'),
            container = cms.string('EcalTPGFineGrainStripEE')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainTowerEERcd'),
            container = cms.string('EcalTPGFineGrainTowerEE')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGLutIdMapRcd'),
            container = cms.string('EcalTPGLutIdMap')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGWeightIdMapRcd'),
            container = cms.string('EcalTPGWeightIdMap')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGWeightGroupRcd'),
            container = cms.string('EcalTPGWeightGroup')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGLutGroupRcd'),
            container = cms.string('EcalTPGLutGroup')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGFineGrainEBGroupRcd'),
            container = cms.string('EcalTPGFineGrainEBGroup')
        ), 
        cms.PSet(
            record = cms.string('EcalTPGPhysicsConstRcd'),
            container = cms.string('EcalTPGPhysicsConst')
        ),
	cms.PSet(
            record = cms.string('EcalTPGCrystalStatusRcd'),
            container = cms.string('EcalTPGCrystalStatus')
        ),
	cms.PSet(
            record = cms.string('EcalTPGTowerStatusRcd'),
            container = cms.string('EcalTPGTowerStatus')
        ))
)

process.p = cms.Path(process.dbCopy)



