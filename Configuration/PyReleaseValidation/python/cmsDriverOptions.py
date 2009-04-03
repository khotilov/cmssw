#! /usr/bin/env python

# A Pyrelval Wrapper

import optparse
import sys
import os
import Configuration.PyReleaseValidation
from Configuration.PyReleaseValidation.ConfigBuilder import ConfigBuilder, defaultOptions
import traceback
# Prepare a parser to read the options
usage=\
"""%prog <TYPE> [options].
Example:

%prog reco -s RAW2DIGI,RECO --conditions FrontierConditions_GlobalTag,STARTUP_V4::All --eventcontent RECOSIM
"""
parser = optparse.OptionParser(usage)

expertSettings = optparse.OptionGroup(parser, '===============\n  Expert Options', 'Caution: please use only if you know what you are doing.')
famosSettings = optparse.OptionGroup(parser, '===============\n  FastSimulation options', '')
parser.add_option_group(expertSettings)

parser.add_option("-s", "--step",
                   help="The desired step. The possible values are: "+\
                        "GEN,SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO,POSTRECO,DQM,ALCA,VALIDATION,HARVESTING or ALL.",
                   default="ALL",
                   dest="step")

parser.add_option("--conditions",
                  help="What conditions to use. Default are frontier conditions 'STARTUP_V4::All' (FrontierConditions_GlobalTag,STARTUP_V4::All)",
                  default=defaultOptions.conditions,
                  dest="conditions")

parser.add_option("--eventcontent",
                  help="What event content to write out. Default=FEVTDEBUG, or FEVT (for cosmics)",
                  default='RECOSIM',
                  dest="eventcontent")

parser.add_option("--filein",
                   help="The infile name.",
                   default="",#to be changed in the default form later
                   dest="filein")

parser.add_option("--fileout",
                   help="The outfile name. If absent a default value is assigned",
                   default="", #to be changed in the default form later
                   dest="fileout")

parser.add_option("-n", "--number",
                  help="The number of events. The default is 1.",
                  default="1",
                  dest="number")

parser.add_option("--mc",
                  help="Specify that simulation is to be processed (default = guess based on options",
                  action="store_true",
                  default=False,
                  dest="isMC")
parser.add_option("--data",
                  help="Specify that data is to be processed (default = guess based on options",
                  action="store_true",
                  default=False,
                  dest="isData")

# expert settings
expertSettings.add_option("--beamspot",
                          help="What beam spot to use (from Configuration/StandardSequences). Default depends on scenario",
                          default=None,
                          dest="beamspot")

expertSettings.add_option("--customise",
                          help="Specify the file where the code to modify the process object is stored.",
                          default="",
                          dest="customisation_file")

expertSettings.add_option("--datatier",
                          help="What data tier to use.",
                          default='',
                          dest="datatier")

expertSettings.add_option( "--dirin",
                          help="The infile directory.",
                          default="",
                          dest="dirin")                    

expertSettings.add_option( "--dirout",
                          help="The outfile directory.",
                          default="",
                          dest="dirout")                

expertSettings.add_option("--filtername",
                          help="What filter name to specify in output module",
                          default="",
                          dest="filtername")

expertSettings.add_option("--geometry",
                          help="What geometry to use (from Configuration/StandardSequences). Default=Ideal",
                          default=defaultOptions.geometry,
                          dest="geometry")

expertSettings.add_option("--magField",
                          help="What magnetic field to use (from Configuration/StandardSequences).",
                          default=defaultOptions.magField,
                          dest="magField")

expertSettings.add_option("--no_output",
                          help="Do not write anything to disk. This is for "+\
                          "benchmarking purposes.",
                          action="store_true",
                          default=False,
                          dest="no_output_flag")

expertSettings.add_option("--oneoutput",
                          help="use only one output module",
                          action="store_true",
                          default="False",
                          dest="oneoutput")

expertSettings.add_option("--prefix",
                          help="Specify a prefix to the cmsRun command.",
                          default="",
                          dest="prefix")  

expertSettings.add_option("--relval",
                          help="Set total number of events and events per job.", #this does not get used but get parsed in the command by DatOps
                          default="",
                          dest="relval")

expertSettings.add_option("--dump_python",
                  help="Dump the config file in python "+\
                  "and do a full expansion of imports.",
                  action="store_true",
                  default=False,                  
                  dest="dump_python")

expertSettings.add_option("--dump_DSetName",
                          help="Dump the primary datasetname.",
                          action="store_true",
                          default=False,
                          dest="dump_dsetname_flag")

expertSettings.add_option("--pileup",
                  help="What pileup config to use. Default=NoPileUp.",
                  default=defaultOptions.pileup,
                  dest="pileup")

                                                    
expertSettings.add_option("--python_filename",
                          help="Change the name of the created config file ",
                          default='',
                          dest="python_filename")

expertSettings.add_option("--secondfilein",
                          help="The secondary infile name."+\
                                "for the two-file solution. Default is no file",
                          default="",#to be changed in the default form later
                          dest="secondfilein")

expertSettings.add_option("--writeraw",
                          help="In addition to the nominal output, write a file with just raw",
                          action="store_true",
                          default=False,
                          dest="writeraw")

expertSettings.add_option("--processName",
                          help="set process name explicitly",
                          default = None,
                          dest="name" 
                          )

expertSettings.add_option("--scenario",
                          help="Select scenario overriding standard settings (available:"+str(defaultOptions.scenarioOptions)+")",
                          default='pp',
                          dest="scenario")

expertSettings.add_option("--harvesting",
                          help="What harvesting to use (from Configuration/StandardSequences). Default=AtRunEnd",
                          default=defaultOptions.harvesting,
                          dest="harvesting")

parser.add_option("--no_exec",
                  help="Do not exec cmsRun. Just prepare the python config file.",
                  action="store_true",
                  default=False,
                  dest="no_exec_flag")   
                  
(options,args) = parser.parse_args() # by default the arg is sys.argv[1:]

# A simple check on the consistency of the arguments
if len(sys.argv)==1:
    raise "Event Type: ", "No event type specified!"

options.evt_type=sys.argv[1]

# memorize the command line arguments 
options.arguments = reduce(lambda x, y: x+' '+y, sys.argv[1:])

# now adjust the given parameters before passing it to the ConfigBuilder



# Build the IO files if necessary.
# The default form of the files is:
# <type>_<energy>_<step>.root
prec_step = {"ALL":"",
             "GEN":"",
             "SIM":"GEN",
             "DIGI":"SIM",
             "RECO":"DIGI",
             "ALCA":"RECO",
             "ANA":"RECO",
             "DIGI2RAW":"DIGI",
             "RAW2DIGI":"DIGI2RAW",
             "HARVESTING":"RECO"}

trimmedEvtType=options.evt_type.split('/')[-1]

trimmedStep=''
isFirst=0
step_list=options.step.split(',')
for s in step_list:
    stepSP=s.split(':')
    step=stepSP[0]
    if ( isFirst==0 ):
        trimmedStep=step
        isFirst=1
    else:
        trimmedStep=trimmedStep+','+step

        
first_step=trimmedStep.split(',')[0]             
if options.filein=="" and not first_step in ("ALL","GEN","SIM_CHAIN"):
    if options.dirin=="":
        options.dirin="file:"
    options.filein=trimmedEvtType+"_"+prec_step[first_step]+".root"


# Prepare the canonical file name for output / config file etc
#   (EventType_STEP1_STEP2_..._PU)
standardFileName = ""
standardFileName = trimmedEvtType+"_"+trimmedStep
standardFileName = standardFileName.replace(",","_").replace(".","_")
if options.pileup != "NoPileUp":
    standardFileName += "_PU"


# if no output file name given, set it to default
if options.fileout=="" and not first_step in ("HARVESTING"):
    options.fileout = standardFileName+".root"


# Prepare the name of the config file
# (in addition list conditions in name)
python_config_filename = standardFileName
conditionsSP = options.conditions.split(',')

if len(conditionsSP) > 1:
    # for conditions like STARTUP_V1, IDEAL_V1 we want only the STARTUP or IDEAL part
    conditionsType = conditionsSP[1].split("_")[0]
    python_config_filename += "_"+str(conditionsType)

python_config_filename+=".py"


#if desired, just add _rawonly to the end of the output file name
fileraw=''
if options.writeraw:
    fileraw=options.dirout
    wrSP=options.fileout.split('.')
    wrSPLen=len(wrSP)
    counter=0
    for w in wrSP:
        counter=counter+1
        if ( counter < wrSPLen ):
            if ( counter == 1):
                fileraw=fileraw+w
            else:    
                fileraw=fileraw+'.'+w
        else:
            fileraw=fileraw+'_rawonly.'+w



secondfilestr=''
if options.secondfilein!='':
    secondfilestr=options.dirin+options.secondfilein


# replace step aliases by right list
if options.step=='ALL':
        options.step='GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,RECO,POSTRECO,VALIDATION,DQM'
elif options.step=='DATA_CHAIN':
        options.step='RAW2DIGI,RECO,POSTRECO,DQM'
options.step = options.step.replace("SIM_CHAIN","GEN,SIM,DIGI,L1,DIGI2RAW")

# add on the end of job sequence...
# if not fastsim or harvesting...

addEndJob = True
if ("FASTSIM" in options.step and not "VALIDATION" in options.step) or "HARVESTING" in options.step: 
    addEndJob = False
if addEndJob:    
    options.step=options.step+',ENDJOB'
print options.step


# Setting name of process
# if not set explicitly it needs some thinking
if not options.name:
    if 'HLT' in trimmedStep:    
        options.name = 'HLT'
    elif 'RECO' in trimmedStep:
        options.name = 'RECO'
    else:
        options.name = trimmedStep.split(',')[-1]

# check to be sure that people run the harvesting as a separte step
isHarvesting = False
isOther = False

s_list=options.step.split(',')
if "HARVESTING" in options.step and len(s_list) > 1:
    print "The Harvesting step must be run alone"
    sys.exit(1)

# sanity check options specifying data or mc
if options.isData and options.isMC:
    print "You may specify only --data or --mc, not both"
    sys.exit(1)

# if not specified by user try to guess
if not options.isData and not options.isMC:
    if 'SIM' in trimmedStep:
        options.isMC=True
    if 'DIGI' in trimmedStep:
        options.isMC=True
    if (not (options.eventcontent == None)) and 'SIM' in options.eventcontent:
        options.isMC=True
    if 'SIM' in options.datatier:
        options.isMC=True
    if options.isMC:
        print 'We have determined that this is simulation (if not, rerun cmsDriver.py with --data)'
    else:
        print 'We have determined that this is real data (if not, rerun cmsDriver.py with --mc)'
    
options.outfile_name = options.dirout+options.fileout

