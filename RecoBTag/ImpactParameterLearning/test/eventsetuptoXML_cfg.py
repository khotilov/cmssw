# The following comments couldn't be translated into the new config version:

#! /bin/env cmsRun

#include "CondCore/DBCommon/data/CondDBCommon.cfi"
#replace CondDBCommon.connect = "sqlite_file:btagnew.db"
#replace CondDBCommon.catalog = "file:mycatalog.xml"
#        es_source = PoolDBESSource {
#                                  using CondDBCommon
#                                 VPSet toGet = {
#                                   {string record = "BTagTrackProbability2DRcd"
#                                     string tag = "probBTagPDF2D_tag"    },
#                                   {string record = "BTagTrackProbability3DRcd"
#                                     string tag = "probBTagPDF3D_tag"    }
#                                    }
#                                   }

import FWCore.ParameterSet.Config as cms

process = cms.Process("analyzer")
#
# one of the three
#

# 1 - frontier
#
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cfi")
#process.GlobalTag.globaltag = "IDEAL_V5::All"

# 2 - fake
#
#process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")


# 3 - file    ---  edit the file position
#
process.load("RecoBTag.TrackProbability.trackProbabilityFakeCond_cfi")
process.trackProbabilityFakeCond.connect = "sqlite_fip:RecoBTag/ImpactParameterLearning/test/btagnew_new_JANDREA.db"

# Message Logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(0)
)
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(1)
)

process.ipCalib = cms.EDFilter("ImpactParameterCalibration",
    writeToDB = cms.bool(False),
    writeToBinary = cms.bool(False),
    nBins = cms.int32(10000),
    resetHistograms = cms.bool(False),
    maxSignificance = cms.double(50.0),
    writeToRootXML = cms.bool(True),
    tagInfoSrc = cms.InputTag("impactParameterTagInfos"),
    inputCategories = cms.string('EventSetup'),
    primaryVertexSrc = cms.InputTag("offlinePrimaryVerticesFromCTFTracks")
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    authenticationMethod = cms.untracked.uint32(1),
    loadBlobStreamer = cms.untracked.bool(True),
    catalog = cms.untracked.string('file:mycatalog_new.xml'),
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0),
        authenticationPath = cms.untracked.string('.')
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string('sqlite_file:btagnew_new.db'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('BTagTrackProbability2DRcd'),
        tag = cms.string('probBTagPDF2D_tag')
    ), 
        cms.PSet(
            record = cms.string('BTagTrackProbability3DRcd'),
            tag = cms.string('probBTagPDF3D_tag')
        ))
)

process.p = cms.Path(process.ipCalib)

