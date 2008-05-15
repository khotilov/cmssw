import FWCore.ParameterSet.Config as cms

# Full configuration for Tracker Geometry Comparison Tool
TrackerGeometryCompare = cms.EDFilter("TrackerGeometryCompare",
    writeToDB = cms.untracked.bool(False),
    outputFile = cms.untracked.string('output.root'),
    setCommonTrackerSystem = cms.untracked.string('NONE'), ##must be "NONE" if you don't want to use this option

    detIdFlag = cms.untracked.bool(False),
    detIdFlagFile = cms.untracked.string('blah.txt'),
    weightById = cms.untracked.bool(False),
    #	untracked vstring levels = {"PixelEndcap","PixelHalfBarrel","TID","HalfBarrel","Endcap","DetUnit"}
    levels = cms.untracked.vstring('Det'),
    weightBy = cms.untracked.string('DetUnit'),
    weightByIdFile = cms.untracked.string('blah2.txt'),
    treeName = cms.untracked.string('alignTree'),
    inputROOTFile1 = cms.untracked.string('IDEAL'),
    inputROOTFile2 = cms.untracked.string('idealtracker2.root')
)


