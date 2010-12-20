import FWCore.ParameterSet.Config as cms

print "The pile up is taken from 8 TeV MinBias files."
# Take pileup events from files
PileUpSimulatorBlock = cms.PSet(
    PileUpSimulator = cms.PSet(
        # The file with the last minimum bias events read in the previous run
        # to be put in the local running directory (if desired)
        inputFile = cms.untracked.string('PileUpInputFile.txt'),
        # Special files of minimum bias events (generated with 
        # cmsRun FastSimulation/PileUpProducer/test/producePileUpEvents8TeV_cfg.py)
        fileNames = cms.untracked.vstring(
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_001.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_002.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_003.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_004.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_005.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_006.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_007.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_008.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_009.root', 
            'file:/afs/cern.ch/cms/data/CMSSW/FastSimulation/PileUpProducer/data/MinBias8TeV_010.root'),
        averageNumber = cms.double(0.0)
    )
)

