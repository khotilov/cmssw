#
# $Id: dump_MC_42X.py,v 1.2 2011/05/21 12:13:47 meridian Exp $
#
#  configuration to dump ntuples in MC
#   the only diff should be for jetmet corrections
#
import FWCore.ParameterSet.Config as cms
import sys
import os
import imp

# contains all details of what to run and defines the process and the path
#from  JetMETCorrections.GammaJet.dumper_process_cff import process

# read in process from file
filename  = 'common_dump_config.py'
handle = open(filename, 'r')
cfo = imp.load_source("pycfg", filename, handle)
process = cfo.process
handle.close()

process.p = cms.Path(process.analysisSequence)

## DO NOT CHANGE THE PATH HERE! New modules should be added ONLY in the common configuration 
#  only paramaters should be changes for data and MC
process.source.fileNames = cms.untracked.vstring(
    'file:/cmsrm/pc24_2/meridian/data/GluGluToHToGG_M-115_7TeV-powheg-pythia6_AODSIM.root'
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Global tag#
# Remember to change it and adapt to your needs #
process.GlobalTag.globaltag = cms.string('START42_V12::All')

##  apply only L2 and L3 jet corrections in MC
process.metMuonJESCorAK5.corrector = "ak5CaloL2L3"

process.myanalysis.JetCorrectionService_akt5 = cms.string('ak5CaloL1L2L3')
process.myanalysis.JetCorrectionService_akt7 = cms.string('ak7CaloL1L2L3')
process.myanalysis.JetCorrectionService_jptak5 = cms.string('ak5JPTL1L2L3')
process.myanalysis.JetCorrectionService_jptak7 = cms.string('ak7JPTL1L2L3')
process.myanalysis.JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3')
process.myanalysis.JetCorrectionService_pfakt7 = cms.string('ak7PFL1FastL2L3')
#process.myanalysis.TriggerTag = cms.untracked.InputTag("TriggerResults::REDIGI311X")

#print process.dumpPython()
