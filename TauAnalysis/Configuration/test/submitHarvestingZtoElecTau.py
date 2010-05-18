#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.Configuration.makeReplacementsHarvesting import makeReplacementsHarvesting

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputFilePath = "/castor/cern.ch/user/j/jkolb/elecTauAnalysis/spring10/hists"

inputFilePath = "rfio:" + outputFilePath

#--------------------------------------------------------------------------------
#
# Add histograms, numbers in FilterStatisticsTables and run + event numbers
# stored as DQM MonitorElements in different ROOT files
#
# NOTE: The jobs get submitted to the '1nh' queue,
#       which allows for an execution time of the cmsRun jobs of up to 1 hour
#       (the queues are {'1nh' (1 hour), '1nd' (24 hours) and '1nw' (1 week execution time limit);
#        see https://twiki.cern.ch/twiki/bin/view/CMS/CMSUKCMSSWBatch for details about the CERN batch system)           
#
#--------------------------------------------------------------------------------

# 7TeV samples

# harvest data
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "Data_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "data")

# harvest Z --> tau tau 
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "Ztautau_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest Z --> e e
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "Zee_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest QCD_BCtoE 
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt20to30_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt30to80_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt80to170_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest QCD_EMenriched
for i in range(2):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt20to30_7TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")
for i in range(4):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt30to80_7TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")
for i in range(2):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt80to170_7TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")

# harvest W/Z + jets
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "WplusJets_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "ZtautauPlusJets_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "ZeePlusJets_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest TT + jets
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "TTplusJets_7TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# 10 TeV samples

# harvest Z --> tau tau 
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "Ztautau_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest Z --> e e
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "Zee_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest QCD_BCtoE 
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt20to30_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt30to80_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "QCD_BCtoE_Pt80to170_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest QCD_EMenriched
for i in range(2):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt20to30_10TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")
for i in range(3):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt30to80_10TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")
for i in range(2):
    submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", 
				  sample = "QCD_EMenriched_Pt80to170_10TeV" + "_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "harvesting", queue = "8nh", outputFilePath = outputFilePath, type = "mc")

# harvest W/Z + jets
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "WplusJets_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "ZtautauPlusJets_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "ZeePlusJets_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

# harvest TT + jets
submitToBatch(configFile = "harvestZtoElecTauPlots_cfg.py", channel = "ZtoElecTau", sample = "TTplusJets_10TeV",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "harvesting", queue = "1nh", outputFilePath = outputFilePath, type = "mc")

