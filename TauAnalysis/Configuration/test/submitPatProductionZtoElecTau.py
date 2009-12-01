#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.Configuration.makeReplacementsPatProduction import makeReplacementsPatProduction

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputFilePath = "/castor/cern.ch/user/j/jkolb/elecTauAnalysis/patTuples/"

# small cmsRun job for testing purposes...
#submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_Ztautau_part01",
#              replFunction = makeReplacementsPatProduction, replacements = "maxEvents = 100",
#              job = "PatProduction", queue = "1nh", outputFilePath = outputFilePath)

#--------------------------------------------------------------------------------
#
# Monte Carlo samples from Summer'08 production
# reprocessed with CMSSW_2_2_3, skimmed by Jeff Kolb
#
# NOTE: The jobs get submitted to the '1nd' queue,
#       which allows for an execution time of the cmsRun jobs of up to 24 hours
#       (the queues are {'1nh' (1 hour), '1nd' (24 hours) and '1nw' (1 week execution time limit);
#        see https://twiki.cern.ch/twiki/bin/view/CMS/CMSUKCMSSWBatch for details about the CERN batch system)           
#
#--------------------------------------------------------------------------------

# Z --> tau tau jobs
for i in range(10):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_Ztautau_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)

# Z --> e e jobs
for i in range(24):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_Zee_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)

# Photon + jets jobs
submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_PhotonJets_Pt15to20",
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "8nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_PhotonJets_Pt20to25",
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "8nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_PhotonJets_Pt25to30",
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "8nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_PhotonJets_Pt30to35",
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "8nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", sample = "ZtoElecTau_PhotonJets_PtGt35",
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "8nh", outputFilePath = outputFilePath)

# QCD_BCtoE jobs
for i in range(28):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
		sample = "ZtoElecTau_QCD_BCtoE_Pt20to30_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)
for i in range(46):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
		sample = "ZtoElecTau_QCD_BCtoE_Pt30to80_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)
for i in range(26):
		submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
		sample = "ZtoElecTau_QCD_BCtoE_Pt80to170_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)    

# QCD_EMenriched jobs
for i in range(53):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
		sample = "ZtoElecTau_QCD_EMenriched_Pt20to30_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)
for i in range(220):    
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
	sample = "ZtoElecTau_QCD_EMenriched_Pt30to80_part%(i)02d" % {"i" : (i + 1)},
	replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
	job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)
for i in range(59):    
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau",
		sample = "ZtoElecTau_QCD_EMenriched_Pt80to170_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)    

# W/Z + jets jobs
for i in range(34):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_WplusJets_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)

for i in range(16):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_ZeePlusJets_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_ZtautauPlusJets_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)

# TT + jets jobs
for i in range(32):
	submitToBatch(configFile = "producePatTuple_cfg.py", channel = "ZtoElecTau", 
		sample = "ZtoElecTau_TTplusJets_part%(i)02d" % {"i" : (i + 1)},
		replFunction = makeReplacementsPatProduction, replacements = "maxEvents = -1",
		job = "PatProduction", queue = "1nd", outputFilePath = outputFilePath)

