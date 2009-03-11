# Python script generated by stageEnergy.py
# For Testbeam data processing

# Jamie Ballin, Imperial College London
# December 2008

import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")
process.load("RecoParticleFlow.PFAnalyses.pflowProcessTestbeam_cff")
from RecoParticleFlow.PFAnalyses.RunDict import *
import FWCore.ParameterSet.VarParsing as VarParsing

# setup 'standard'  options
options = VarParsing.VarParsing()
options.register ('beamEnergy',
                  100, # default value
                  options.multiplicity.singleton, # singleton or list
                  options.varType.int,          # string, int, or float
                  "Beam energy to simulate")

options.register ('fileSuffix',
                  '', # default value
                  options.multiplicity.singleton, # singleton or list
                  options.varType.string,          # string, int, or float
                  "Label to append to output file names")

# setup any defaults you want
options.beamEnergy = 100
options.fileSuffix = ''

# get and parse the command line arguments
options.parseArguments()

suffix = ''
if not options.fileSuffix ==  '':
    suffix = '_' + options.fileSuffix 

outputTree = "outputtree_" + str(options.beamEnergy) + "GeV" + suffix + ".root"
outputFile = "reprocessed_" + str(options.beamEnergy) + "GeV" + suffix + ".root"
logFile = "log_" + str(options.beamEnergy) + "GeV" + suffix + ".txt"

print ("pflow_tb_general.py with options:")
print "Beam energy: " + str(options.beamEnergy)
print "File suffix: " + str(options.fileSuffix)
print "Output file: " + outputFile
print "Output tree: " + outputTree

specifiedE = energies[options.beamEnergy]
result = map(lambda x : 'rfio:///castor/cern.ch/cms/store/h2tb2006/reco/v6/h2.000' + str(x) + '.combined.OutServ_0.0-cmsswreco.root', specifiedE)


# Files to process
runs = cms.untracked.vstring(result)

# Output tree of cleaned particles
process.TFileService.fileName = cms.string(outputTree)

# New Event file
process.finishup.fileName = cms.untracked.string(outputFile)

# LogFile

process.MessageLogger = cms.Service("MessageLogger",
    destinations=cms.untracked.vstring(logFile)
)

process.source = cms.Source("PoolSource",
        fileNames = runs
)

process.p1 = cms.Path(process.pflowProcessTestbeam)
process.outpath = cms.EndPath(process.finishup)
