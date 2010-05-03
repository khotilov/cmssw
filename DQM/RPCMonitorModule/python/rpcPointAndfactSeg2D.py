# The following comments couldn't be translated into the new config version:

#keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSegmentEff")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.load("DQMServices.Core.DQM_cfg")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_31X_V1::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/data/CRAFT09/Cosmics/RECO/v1/000/111/125/8A62785E-158D-DE11-884B-0030486730C6.root')
)

process.rpcPointProducer = cms.EDProducer('RPCPointProducer',
  incldt = cms.untracked.bool(True),
  inclcsc = cms.untracked.bool(True),

  debug = cms.untracked.bool(False),

  rangestrips = cms.untracked.double(4.),
  rangestripsRB4 = cms.untracked.double(4.),
  MinCosAng = cms.untracked.double(0.85),
  MaxD = cms.untracked.double(80.0),
  MaxDrb4 = cms.untracked.double(150.0),
  ExtrapolatedRegion = cms.untracked.double(0.6), #in stripl/2 in Y and stripw*nstrips/2 in X

##  cscSegments = cms.untracked.string('hltCscSegments'),
##  dt4DSegments = cms.untracked.string('hltDt4DSegments'),
  cscSegments = cms.untracked.string('cscSegments'),
  dt4DSegments = cms.untracked.string('dt4DSegments')
)

process.museg = cms.EDFilter("MuonSegmentEff",

    incldt = cms.untracked.bool(True),
    incldtMB4 = cms.untracked.bool(True),
    inclcsc = cms.untracked.bool(True),

    debug = cms.untracked.bool(False),
    inves = cms.untracked.bool(True),
    
    DuplicationCorrection = cms.untracked.int32(1),

    manualalignment = cms.untracked.bool(True),
    AliFileName = cms.untracked.string('/afs/cern.ch/user/c/carrillo/endcap/CMSSW_3_0_0_pre10/src/DQM/RPCMonitorModule/data/Alignment69912.dat'),

    rangestrips = cms.untracked.double(4.),

##  cscSegments = cms.untracked.string('hltCscSegments'),
##  dt4DSegments = cms.untracked.string('hltDt4DSegments'),
    cscSegments = cms.untracked.InputTag('cscSegments'),
    dt4DSegments = cms.untracked.InputTag('dt4DSegments'),


## is missing the hltrpcRecHits
    rpcRecHits = cms.InputTag("rpcRecHits"),
    rpcDTPoints = cms.InputTag("rpcPointProducer","RPCDTExtrapolatedPoints"),
    rpcCSCPoints = cms.InputTag("rpcPointProducer","RPCCSCExtrapolatedPoints"),

    EffSaveRootFile = cms.untracked.bool(True),
    EffRootFileName = cms.untracked.string('/tmp/carrillo/eff111125.FEEFDE56-288D-DE11-9A6E-0030486780B8.root'),
    EffSaveRootFileEventsInterval = cms.untracked.int32(100)
)

process.p = cms.Path(process.rpcPointProducer*process.museg)

process.DQM.collectorHost = ''
process.DQM.collectorPort = 9090
process.DQM.debug = False


