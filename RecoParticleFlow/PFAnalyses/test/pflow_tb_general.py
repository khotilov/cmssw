# Python script generated by stageEnergy.py
# For Testbeam data processing

# Jamie Ballin, Imperial College London
# December 2008

import FWCore.ParameterSet.Config as cms
import commands
process = cms.Process("SKIM")
process.load("RecoParticleFlow.PFAnalyses.pflowProcessTestbeam_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
from RecoParticleFlow.PFAnalyses.RunDict import *

from RecoParticleFlow.PFAnalyses.pflowOptions_cfi import *

outputTree = "PFlowTB_Tree_" + fileLabel
outputFile = "PFlowTB_Events_" + fileLabel
logFile = "log_" + logLabel


if options.notracks <> 0:
    process.faketracks.justCreateEmptyCollections = cms.bool(True)
    print "Running in notrack mode"

specifiedE = energies[options.beamEnergy]
if options.endcapMode <> 0:
    specifiedE = endcap[options.beamEnergy]


process.particleFiltration.debug = cms.int32(1)

process.particleFlowRecHitHCAL.thresh_Barrel = cms.double(0.0)
process.particleFlowRecHitHCAL.thresh_Endcap = cms.double(0.0)

#For calibration purposes
process.particleFlow.pf_nsigma_HCAL = cms.double(1.0)
process.particleFlow.pf_calibMode = cms.uint32(1)
process.particleFlowBlock.pf_chi2_ECAL_HCAL = cms.double(100.0)
#process.particleFlowBlock.debug = cms.untracked.bool(True)
#process.particleFlow.debug = cms.untracked.bool(True)
#Uncomment this lot if you want a file of noise!
#process.particleFiltration.noiseMode=cms.bool(True)
#process.extraction.applyCleaningCuts=cms.bool(False)
#process.extraction.saveJustPions=cms.bool(False)
## Also change masterConeDeltaR in pflowCalibratable_cfi to something like 1

#Need to override clustering to exclude HF components
from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *
process.pfClusteringHCAL = cms.Sequence(particleFlowRecHitHCAL * particleFlowClusterHCAL)

result = []
if options.copyToTmp <> 0:
    if options.endcapMode <> 0:
          map(lambda a : commands.getoutput('rfcp /castor/cern.ch/cms/store/data/h2tb2007/testbeam_HCalEcalCombined/DIGI-RECO/default_v1/tb07_reco_edm_run_000' + str(a) + '.root' + ' /tmp'), specifiedE)
          result = map(lambda x : 'file:/tmp/tb07_reco_edm_run_000' + str(x) + '.root', specifiedE)
    else:
        map(lambda a : commands.getoutput('rfcp ' + '/castor/cern.ch/cms/store/h2tb2006/reco/v6/h2.000' + str(a) + '.combined.OutServ_0.0-cmsswreco.root' + ' /tmp'), specifiedE)
        result = map(lambda x : 'file:/tmp/h2.000' + str(x) + '.combined.OutServ_0.0-cmsswreco.root', specifiedE)
else:
    if options.endcapMode <> 0: 
        result = map(lambda x : 'rfio:///castor/cern.ch/cms/store/data/h2tb2007/testbeam_HCalEcalCombined/DIGI-RECO/default_v1/tb07_reco_edm_run_000' + str(x) + '.root', specifiedE)
    else:   
        result = map(lambda x : 'rfio:///castor/cern.ch/cms/store/h2tb2006/reco/v6/h2.000' + str(x) + '.combined.OutServ_0.0-cmsswreco.root', specifiedE)
    
if options.endcapMode <> 0:
    process.particleFlowRecHitECAL.ecalRecHitsEB = cms.InputTag("pflowCalibEcalRechits", "ecalEBRechitsCalib")
    process.particleFlowRecHitECAL.ecalRecHitsEE = cms.InputTag("pflowCalibEcalRechits", "ecalEERechitsCalib")
    process.extraction.RawRecHitsEcalEB = cms.InputTag("pflowCalibEcalRechits", "ecalEBRechitsCalib")
    process.extraction.RawRecHitsEcalEE = cms.InputTag("pflowCalibEcalRechits", "ecalEERechitsCalib")
    process.faketracks.endcapMode = cms.bool(True)
    process.particleFiltration.isEndcap2007 = cms.bool(True)
    process.particleFlowRecHitHCAL.isEndcap2007 = cms.bool(True)
    process.extraction.isEndcap2007 = cms.bool(True)
    
process.extraction.clustersFromCandidates=cms.bool(False)
process.extraction.rechitsFromCandidates=cms.bool(False)

if options.kevents <> 0:
    process.maxEvents = cms.untracked.PSet(
        input=cms.untracked.int32(options.kevents * 1000)
    )

# Files to process
runs = cms.untracked.vstring(result)
print "Input files :"
print result

print "Log file :"
print logFile

# Output tree of cleaned particles
process.TFileService.fileName = cms.string(outputTree)

# New Event file
process.finishup.fileName = cms.untracked.string(outputFile)


# LogFile
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.destinations=cms.untracked.vstring('PFlowTB_' + logLabel, 'cout')
#process.MessageLogger.destinations=cms.untracked.vstring('cout')


process.source = cms.Source("PoolSource",
        fileNames=runs,
        inputCommands=cms.untracked.vstring('keep *', 'drop EBDataFramesSorted_*_*_*', 'drop EEDataFramesSorted_*_*_*')

        
)
#process.p1 = cms.Path(process.pflowCleaning)
process.p1 = cms.Path(process.pflowProcessTestbeam)

if options.endcapMode <> 0:
    process.p1 = cms.Path(process.pflowProcessEndcapTestbeam)

if options.outputCollections:
    process.outpath = cms.EndPath(process.finishup)
