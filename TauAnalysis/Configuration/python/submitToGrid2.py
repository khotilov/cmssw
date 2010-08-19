import os
import subprocess
import copy
import string
from TauAnalysis.Configuration.cfgOptionMethods import copyCfgFileAndApplyOptions

_CRAB_TEMPLATE = string.Template('''
[CRAB]
jobtype = cmssw
use_server = 1
scheduler = glidein

[CMSSW]
datasetpath = $datasetpath
dbs_url = $dbs_url
pset = $pset
total_number_of_events = $total_number_of_events
events_per_job = $events_per_job
output_file = $output_file

[USER]
ui_working_dir = $ui_working_dir
return_data = $return_data
copy_data = $copy_data
storage_element = $storage_element
storage_path = $storage_path
user_remote_dir = $user_remote_dir
publish_data = $publish_data

[GRID]
# todo
''')

_CRAB_DEFAULTS = {
    'events_per_job' : 20000,
    'total_number_of_events' : -1,
    'return_data' : 0,
    'copy_data' : 1,
    'storage_element' : 'srm-cms.cern.ch',
    'storage_path' : '/srm/managerv2?SFN=/castor/cern.ch',
    'publish_data' : 0,
}

def submitToGrid(configFile, jobInfo, jobOptions, crabOptions, submit=True):
    workingDirectory = os.getcwd()
    submissionDirectory = os.path.join(workingDirectory, 'crab')
    configFilePath = os.path.join(workingDirectory, configFile)
    if not os.path.exists(configFilePath):
        raise ValueError("Can't find config file %s in current direcotry!" %
                         configFile)
    # Strip off _cfg.py and add sample info
    configFileName = configFile.replace(
        '_cfg.py', '_%s_%s_%s' % (jobInfo['channel'], jobInfo['sample'], jobInfo['id']))
    # New name of config file
    newConfigFile =  configFileName + '@Grid_cfg.py' 
    newConfigFilePath = os.path.join(submissionDirectory, newConfigFile)

    # Copy the config file and add our specialization options
    copyCfgFileAndApplyOptions(configFilePath, newConfigFilePath,
                               jobInfo, jobOptions)
    # Update the default crab options with our options
    fullCrabOptions = copy.copy(_CRAB_DEFAULTS)
    # Point crab to our PSET
    fullCrabOptions['pset'] = newConfigFilePath
    ui_working_dir = os.path.join(
        submissionDirectory, 'crabdir_%s' % configFileName)
    fullCrabOptions['ui_working_dir'] = ui_working_dir
    fullCrabOptions.update(crabOptions)

    # Create the crab file
    crabFileName = "crab_" + configFileName + ".cfg"
    crabFilePath = os.path.join( submissionDirectory, crabFileName)
    crabFile = open(crabFilePath, 'w')
    crabFile.write(_CRAB_TEMPLATE.substitute(fullCrabOptions))
    crabFile.close()

    crabCreateCommand = "crab -create -cfg " + crabFilePath
    print crabCreateCommand
    subprocess.call(crabCreateCommand, shell = True)
    if submit:
        crabSubmitCommand = "crab -submit -c "+ ui_working_dir
        subprocess.call(crabSubmitCommand, shell = True)
        crabStatusCommand = "crab -status -c "+ ui_working_dir
        subprocess.call(crabStatusCommand, shell = True)

