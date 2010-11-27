import FWCore.ParameterSet.Config as cms

import copy

from TauAnalysis.Configuration.prepareConfigFile2 import getNewConfigFileName, prepareConfigFile2, PLOT_FILES_PREFIX
from TauAnalysis.Configuration.submitToGrid2 import submitToGrid

def _number_of_jobs(sample_info, preferred = 10000, max_jobs = 300):
    """
    Try to run on <preferred> events per job, unless it would
    create to many jobs.
    """
    # Check if explicitly specified
    if 'number_of_jobs' in sample_info:
        return sample_info['number_of_jobs']
    # Otherwise try and run 10k events per job, with a max of 300 jobs
    to_process = sample_info['events_processed']*sample_info['skim_eff']
    desired_njobs = to_process/preferred
    if desired_njobs < 0:
        return max_jobs
    output = (desired_njobs > max_jobs) and max_jobs or desired_njobs
    return int(output)

def submitAnalysisToGrid(configFile = None, channel = None, samples = None,
                         outputFilePath = None, jobId = None,
                         samplesToAnalyze = None, samplesToSkip = None,
                         disableFactorization = False,
                         disableSysUncertainties = False,
                         create = True, submit = True,
                         cfgdir = 'crab', inputFileMap = None,
                         outputFileMap = None,
                         enableEventDumps = False,
                         enableFakeRates = False,
                         processName = None,
                         saveFinalEvents = False):
    """
    Submit analysis job (event selection, filling of histogram)
    via crab
    """

    # check that configFile, channel, samples, outputFilePath and jobId
    # parameters are defined and non-empty
    if configFile is None:
        raise ValueError("Undefined configFile Parameter !!")
    if channel is None:
        raise ValueError("Undefined channel Parameter !!")
    if samples is None:
        raise ValueError("Undefined samples Parameter !!")
    if outputFilePath is None:
        raise ValueError("Undefined outputFilePath Parameter !!")
    if jobId is None:
        raise ValueError("Undefined jobId Parameter !!")

    # Loop over the samples to be analyzed
    for sample in samples['SAMPLES_TO_ANALYZE']:
        # Skip submitting crab job in case
        #  o list of samples for which crab jobs are to be submitted has been
        #    explicitely specified
        #  o sample has explicitely been requested to be skipped
        if samplesToAnalyze:
            if sample not in samplesToAnalyze:
                print "Skipping", sample
                continue
        if samplesToSkip:
            if sample in samplesToSkip:
                print "Skipping", sample
                continue
        print "Submitting ", sample

        sample_info = samples['RECO_SAMPLES'][sample]
        
        # Make job info
        jobInfo = {
            'channel' : channel,
            'sample' : sample,
            'id' : jobId
        }

        newConfigFile = getNewConfigFileName(configFile, cfgdir, sample, jobId)

        # Check if we want to use a special file for the produced cfg file
        # File map is a function that takes a sample name and returns a list of
        # files corresponding to that file.  If files is None, no change will be
        # made.
        input_files = None
        if inputFileMap is not None:
            input_files = inputFileMap(sample)
            if input_files is None:
                print "Warning: No special input files specified for sample%s, using default." % sample
        output_file = None        
        if outputFileMap is not None:
            output_file = outputFileMap(sample)
            
        prepareConfigFile2(
          configFile = configFile, jobInfo = jobInfo, newConfigFile = newConfigFile,
          sample_infos = samples,
          disableFactorization = disableFactorization,
          disableSysUncertainties = disableSysUncertainties,
          input_files = input_files, output_file = output_file,
          enableEventDumps = enableEventDumps, enableFakeRates = enableFakeRates,
          processName = processName,
          saveFinalEvents = saveFinalEvents)

        # Always include the plot files
        output_files = ["%s_%s_%s_%s.root" % (
            PLOT_FILES_PREFIX, jobInfo['channel'],
            jobInfo['sample'], jobInfo['id'])]

        # Add our final event skim as well
        if saveFinalEvents:
            output_files.append("final_events_%s_%s_%s.root" % (
             jobInfo['channel'], jobInfo['sample'], jobInfo['id']))

        # Build crab options
        crabOptions = {
            'number_of_jobs' : _number_of_jobs(sample_info),
            'datasetpath' : sample_info['datasetpath'],
            'dbs_url' : sample_info['dbs_url'],
            'user_remote_dir' : outputFilePath,
            'output_file' : ", ".join(output_files),
            # Default MC info
            'split_type' : (sample_info['type'] == 'Data') and 'lumis' or 'events',
            'lumi_mask' : sample_info['lumi_mask'],
            'runselection' : sample_info['runselection'],
            'SE_white_list' : sample_info['SE_white_list'],
            'SE_black_list' : sample_info['SE_black_list']
        }

        ##submitToGrid(newConfigFile, jobInfo, jobOptions, crabOptions,
        ##             create=create, submit=submit, cfgdir=cfgdir)
        submitToGrid(newConfigFile, jobInfo, crabOptions, create=False, submit=False, cfgdir=cfgdir) # CV: only for testing
