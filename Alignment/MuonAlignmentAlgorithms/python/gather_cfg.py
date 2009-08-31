import os
import FWCore.ParameterSet.Config as cms

inputfiles = os.environ["ALIGNMENT_INPUTFILES"].split(" ")
iteration = int(os.environ["ALIGNMENT_ITERATION"])
jobnumber = int(os.environ["ALIGNMENT_JOBNUMBER"])
mapplots = (os.environ["ALIGNMENT_MAPPLOTS"] == "True")
globaltag = os.environ["ALIGNMENT_GLOBALTAG"]
inputdb = os.environ["ALIGNMENT_INPUTDB"]
trackerconnect = os.environ["ALIGNMENT_TRACKERCONNECT"]
trackeralignment = os.environ["ALIGNMENT_TRACKERALIGNMENT"]
trackerAPE = os.environ["ALIGNMENT_TRACKERAPE"]
iscosmics = (os.environ["ALIGNMENT_ISCOSMICS"] == "True")
station123params = os.environ["ALIGNMENT_STATION123PARAMS"]
station4params = os.environ["ALIGNMENT_STATION4PARAMS"]
cscparams = os.environ["ALIGNMENT_CSCPARAMS"]
minTrackPt = float(os.environ["ALIGNMENT_MINTRACKPT"])
maxTrackPt = float(os.environ["ALIGNMENT_MAXTRACKPT"])
minTrackerHits = int(os.environ["ALIGNMENT_MINTRACKERHITS"])
maxTrackerRedChi2 = float(os.environ["ALIGNMENT_MAXTRACKERREDCHI2"])
allowTIDTEC = (os.environ["ALIGNMENT_ALLOWTIDTEC"] == "True")
twoBin = (os.environ["ALIGNMENT_TWOBIN"] == "True")
minAlignmentHits = int(os.environ["ALIGNMENT_MINALIGNMENTHITS"])

process = cms.Process("GATHER")
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(*inputfiles))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.MessageLogger = cms.Service("MessageLogger",
                                    destinations = cms.untracked.vstring("cout"),
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

process.load("Alignment.MuonAlignmentAlgorithms.MuonAlignmentFromReference_cff")
process.looper.algoConfig.writeTemporaryFile = "alignment%03d.tmp" % jobnumber
process.looper.algoConfig.doAlignment = False

process.looper.algoConfig.combineME11 = False
process.looper.ParameterBuilder.Selector.alignParams = cms.vstring("MuonDTChambers,%s,stations123" % station123params, "MuonDTChambers,%s,station4" % station4params, "MuonCSCChambers,%s" % cscparams)
process.looper.algoConfig.minTrackPt = minTrackPt
process.looper.algoConfig.maxTrackPt = maxTrackPt
process.looper.algoConfig.minTrackerHits = minTrackerHits
process.looper.algoConfig.maxTrackerRedChi2 = maxTrackerRedChi2
process.looper.algoConfig.allowTIDTEC = allowTIDTEC
process.looper.algoConfig.twoBin = twoBin
process.looper.algoConfig.minAlignmentHits = minAlignmentHits

if mapplots:
    process.load("Alignment.CommonAlignmentMonitor.AlignmentMonitorMuonSystemMap1D_cfi")
    process.looper.monitorConfig = cms.PSet(monitors = cms.untracked.vstring("AlignmentMonitorMuonSystemMap1D"),
                                            AlignmentMonitorMuonSystemMap1D = process.AlignmentMonitorMuonSystemMap1D)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(globaltag)
process.looper.applyDbAlignment = True
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

if iscosmics:
    process.Path = cms.Path(process.offlineBeamSpot * process.MuonAlignmentFromReferenceGlobalCosmicRefit)
    process.looper.tjTkAssociationMapTag = cms.InputTag("MuonAlignmentFromReferenceGlobalCosmicRefit:Refitted")
else:
    process.Path = cms.Path(process.offlineBeamSpot * process.MuonAlignmentFromReferenceGlobalMuonRefit)
    process.looper.tjTkAssociationMapTag = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted")

process.MuonAlignmentFromReferenceInputDB.connect = cms.string("sqlite_file:%s" % inputdb)
process.MuonAlignmentFromReferenceInputDB.toGet = cms.VPSet(cms.PSet(record = cms.string("DTAlignmentRcd"), tag = cms.string("DTAlignmentRcd")),
                                                            cms.PSet(record = cms.string("CSCAlignmentRcd"), tag = cms.string("CSCAlignmentRcd")))

if trackerconnect != "":
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.TrackerAlignmentInputDB = cms.ESSource("PoolDBESSource",
                                                   CondDBSetup,
                                                   connect = cms.string(trackerconnect),
                                                   toGet = cms.VPSet())
    if trackeralignment != "":
        process.TrackerAlignmentInputDB.toGet.append(cms.PSet(record = cms.string("TrackerAlignmentRcd"), tag = cms.string(trackeralignment)))
    if trackerAPE != "":
        process.TrackerAlignmentInputDB.toGet.append(cms.PSet(cms.PSet(record = cms.string("TrackerAlignmentErrorRcd"), tag = cms.string(trackerAPE))))

    process.es_prefer_TrackerAlignmentInputDB = cms.ESPrefer("PoolDBESSource", "TrackerAlignmentInputDB")

process.looper.saveToDB = False
del process.PoolDBOutputService

process.TFileService = cms.Service("TFileService", fileName = cms.string("plotting%03d.root" % jobnumber))
