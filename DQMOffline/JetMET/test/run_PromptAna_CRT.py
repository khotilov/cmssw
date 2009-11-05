
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#
# --- [do harvesting (default=True)? or read in histogram files]
harvesting = (os.environ.get('HARVESTING',True))
print 'harvesting (default=True) = '+str(harvesting)
#
# --- [reference histogram (default=jetMETMonitoring_test.root)]
reference_histogram = (os.environ.get('REFERENCE_HISTOGRAM','jetMETMonitoring_test.root'))
print 'reference_histogram = '+str(reference_histogram)
#
# --- [input file(s) for certification (default=reco_DQM_test.root)]
input_root_files = os.environ.get('INPUTFILES','file:reco_DQM_test.root').split(",")
print 'input_root_files = '+str(input_root_files)
print

#-----------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.DQM.collectorHost = ''

#-----------------------------
# DQM Environment & Specify inputs
#-----------------------------
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1)
)

#
#--- When read in RECO file including EDM from ME
process.source = cms.Source("PoolSource",
    processingMode = cms.untracked.string('RunsAndLumis'),
#   fileNames = cms.untracked.vstring('file:reco_DQM_cruzet100945.root')
    fileNames = cms.untracked.vstring(input_root_files)
)

#-----

process.load('Configuration/StandardSequences/EDMtoMEAtRunEnd_cff')
process.dqmSaver.referenceHandling = cms.untracked.string('all')

#-----------------------------
# Specify root file including reference histograms
#-----------------------------
process.DQMStore.referenceFileName = reference_histogram

#-----------------------------
# Locate a directory in DQMStore
#-----------------------------
process.dqmInfoJetMET = cms.EDFilter("DQMEventInfo",
                subSystemFolder = cms.untracked.string('JetMET')
                )

#-----------------------------
# JetMET Certification Module 
#-----------------------------
process.load("DQMOffline.JetMET.dataCertificationJetMET_cff")
process.dataCertificationJetMET = cms.EDAnalyzer('DataCertificationJetMET',
#
#--- Always define reference root file by process.DQMStore.referenceFileName
#    Use process.DQMStore.referenceFileName above. This should be empty.
                              refFileName    = cms.untracked.string(""),
#
#--- 0: harvest EDM files, 1: read in DQM root file
                              TestType       = cms.untracked.int32(0),
#
#--- When read in RECO file including EDM from ME
#                             fileName       = cms.untracked.string("jetMETMonitoring_cruzet98154.root"),
                              fileName       = cms.untracked.string(""),
#
#--- Do note save here. Save output by dqmSaver
                              OutputFile     = cms.untracked.bool(False),
                              OutputFileName = cms.untracked.string(""),
#
                              Verbose        = cms.untracked.int32(0)
)

#-----------------------------
# 
#-----------------------------
process.load("DQMOffline.Trigger.JetMETHLTOfflineClient_cfi")
from DQMOffline.Trigger.JetMETHLTOfflineClient_cfi import *

#-----------------------------
# 
#-----------------------------
#process.p = cms.Path(process.dqmInfoJetMET*process.dataCertificationJetMET)

process.p = cms.Path(process.EDMtoME
                     * process.dqmInfoJetMET
                     * process.jetMETHLTOfflineClient
                     * process.dataCertificationJetMETSequence
                     * process.dqmSaver)

