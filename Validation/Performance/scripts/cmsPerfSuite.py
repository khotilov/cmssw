#!/usr/bin/env python
import os, time, sys, re, glob
import optparse as opt
import cmsRelRegress as crr
from cmsPerfCommons import Candles, MIN_REQ_TS_EVENTS, CandFname, getVerFromLog

global ERRORS 
ERRORS = 0
_CASTOR_DIR = "/castor/cern.ch/cms/store/relval/performance/"
_dryrun   = False
_debug    = False
_unittest = False
_verbose  = True
logh = sys.stdout

try:
    #Get some environment variables to use
    cmssw_base        = os.environ["CMSSW_BASE"]
    cmssw_release_base= os.environ["CMSSW_RELEASE_BASE"]
    cmssw_version     = os.environ["CMSSW_VERSION"]
    host              = os.environ["HOST"]
    user              = os.environ["USER"]
except KeyError:
    print 'Error: An environment variable either CMSSW_{BASE, RELEASE_BASE or VERSION} HOST or USER is not available.'
    print '       Please run eval `scramv1 runtime -csh` to set your environment variables'
    sys.exit()

#Scripts used by the suite:
Scripts         =["cmsDriver.py","cmsRelvalreport.py","cmsRelvalreportInput.py","cmsScimark2"]
AuxiliaryScripts=["cmsScimarkLaunch.csh","cmsScimarkParser.py","cmsScimarkStop.pl"]


#Options handling
def optionParse():
    global _dryrun, _debug, _unittest, _verbose    
    parser = opt.OptionParser(usage='''./cmsPerfSuite.py [options]
       
Examples:
./cmsPerfSuite.py
(this will run with the default options)
OR
./cmsPerfSuite.py -o "/castor/cern.ch/user/y/yourusername/yourdirectory/"
(this will archive the results in a tarball on /castor/cern.ch/user/y/yourusername/yourdirectory/)
OR
./cmsPerfSuite.py -t 5 -i 2 -v 1
(this will run the suite with 5 events for TimeSize tests, 2 for IgProf tests, 1 for Valgrind tests)
OR
./cmsPerfSuite.py -t 200 --candle QCD_80_120 --cmsdriver="--conditions FakeConditions"
(this will run the performance tests only on candle QCD_80_120, running 200 TimeSize evts, default IgProf and Valgrind evts. It will also add the option "--conditions FakeConditions" to all cmsDriver.py commands executed by the suite)
OR
./cmsPerfSuite.py -t 200 --candle QCD_80_120 --cmsdriver="--conditions=FakeConditions --eventcontent=FEVTDEBUGHLT" --step=GEN-SIM,DIGI
(this will run the performance tests only on candle QCD_80_120, running 200 TimeSize evts, default IgProf and Valgrind evts. It will also add the option "--conditions=FakeConditions" and the option "--eventcontent=FEVTDEBUGHLT" to all cmsDriver.py commands executed by the suite. In addition it will run only 2 cmsDriver.py "steps": "GEN,SIM" and "DIGI". Note the syntax GEN-SIM for combined cmsDriver.py steps)

Legal entries for individual candles (--candle option):
%s
''' % ("\n".join(Candles)))

    parser.set_defaults(TimeSizeEvents   = 100        ,
                        IgProfEvents     = 5          ,
                        ValgrindEvents   = 1          ,
                        cmsScimark       = 10         ,
                        cmsScimarkLarge  = 10         ,  
                        cmsdriverOptions = ""         ,
                        stepOptions      = ""         ,
                        candleOptions    = ""         ,
                        profilers        = ""         ,
                        outputdir        = ""         ,
                        logfile          = None       ,
                        runonspare       = True       ,
                        bypasshlt        = False      ,
                        quicktest        = False      ,
                        unittest         = False      ,
                        dryrun           = False      ,
                        verbose          = True       ,
                        previousrel      = ""         ,
                        castordir        = _CASTOR_DIR,
                        cores            = 4          , #Number of cpu cores on the machine
                        cpu              = "1"        ) #Cpu core on which the suite is run:

    parser.add_option('-q', '--quiet'      , action="store_false", dest='verbose'   ,
        help = 'Output less information'                  )
    parser.add_option('-b', '--bypass-hlt' , action="store_true" , dest='bypasshlt' ,
        help = 'Bypass HLT root file as input to RAW2DIGI')
    parser.add_option('-n', '--notrunspare', action="store_false", dest='runonspare',
        help = 'Do not run cmsScimark on spare cores')        
    parser.add_option('-t', '--timesize'  , type='int'   , dest='TimeSizeEvents'  , metavar='<#EVENTS>'   ,
        help = 'specify the number of events for the TimeSize tests'                   )
    parser.add_option('-i', '--igprof'    , type='int'   , dest='IgProfEvents'    , metavar='<#EVENTS>'   ,
        help = 'specify the number of events for the IgProf tests'                     )
    parser.add_option('-v', '--valgrind'  , type='int'   , dest='ValgrindEvents'  , metavar='<#EVENTS>'   ,
        help = 'specify the number of events for the Valgrind tests'                   )
    parser.add_option('--cmsScimark'      , type='int'   , dest='cmsScimark'      , metavar=''            ,
        help = 'specify the number of times the cmsScimark benchmark is run before and after the performance suite on cpu1'         )
    parser.add_option('--cmsScimarkLarge' , type='int'   , dest='cmsScimarkLarge' , metavar=''            ,
        help = 'specify the number of times the cmsScimarkLarge benchmark is run before and after the performance suite on cpu1'    )
    parser.add_option('--cores'           , type='int', dest='cores'              , metavar='<CORES>'     ,
        help = 'specify the number of cores of the machine (can be used with 0 to stop cmsScimark from running on the other cores)' )        
    parser.add_option('-c', '--cmsdriver' , type='string', dest='cmsdriverOptions', metavar='<OPTION_STR>',
        help = 'specify special options to use with the cmsDriver.py commands (designed for integration build use'                  )        
    parser.add_option('-a', '--archive'   , type='string', dest='castordir'       , metavar='<DIR>'       ,
        help = 'specify the wanted CASTOR directory where to store the results tarball'                                             )
    parser.add_option('-L', '--logfile'   , type='string', dest='logfile'         , metavar='<FILE>'      ,
        help = 'file to store log output of the script'                                                                             )                
    parser.add_option('-o', '--output'    , type='string', dest='outputdir'       , metavar='<DIR>'       ,
        help = 'specify the directory where to store the output of the script'                                                      )        
    parser.add_option('-r', '--prevrel'   , type='string', dest='previousrel'     , metavar='<DIR>'       ,
        help = 'Top level dir of previous release for regression analysis'                                                          )        
    parser.add_option('--step'            , type='string', dest='stepOptions'     , metavar='<STEPS>'     ,
        help = 'specify the processing steps intended (instead of the default ones)'                                                )
    parser.add_option('--candle'          , type='string', dest='candleOptions'   , metavar='<CANDLES>'   ,
        help = 'specify the candle(s) to run (instead of all 7 default candles)'                                                    )
    parser.add_option('--cpu'             , type='string', dest='cpu'             , metavar='<CPU>'       ,
        help = 'specify the core on which to run the performance suite'                                                             )


    #####################
    #    
    # Developer options
    #

    devel  = opt.OptionGroup(parser, "Developer Options",
                                     "Caution: use these options at your own risk."
                                     "It is believed that some of them bite.\n")

    devel.add_option('-p', '--profile'  , type="str" , dest='profilers', metavar="<PROFILERS>" ,
        help = 'Profile codes to use for cmsRelvalInput' )
    devel.add_option('-f', '--false-run', action="store_true", dest='dryrun'   ,
        help = 'Dry run'                                                                                           )            
    devel.add_option('-d', '--debug'    , action='store_true', dest='debug'    ,
        help = 'Debug'                                                                                             )
    devel.add_option('--quicktest'      , action="store_true", dest='quicktest',
        help = 'Quick overwrite all the defaults to small numbers so that we can run a quick test of our chosing.' )  
    devel.add_option('--test'           , action="store_true", dest='unittest' ,
        help = 'Perform a simple test, overrides other options. Overrides verbosity and sets it to false.'         )            

    parser.add_option_group(devel)
    (options, args) = parser.parse_args()


    _debug           = options.debug
    _unittest        = options.unittest 
    _verbose         = options.verbose
    _dryrun          = options.dryrun    
    castordir        = options.castordir
    TimeSizeEvents   = options.TimeSizeEvents
    IgProfEvents     = options.IgProfEvents
    ValgrindEvents   = options.ValgrindEvents
    cmsScimark       = options.cmsScimark
    cmsScimarkLarge  = options.cmsScimarkLarge
    cmsdriverOptions = options.cmsdriverOptions
    stepOptions      = options.stepOptions
    quicktest        = options.quicktest
    candleoption     = options.candleOptions
    runonspare       = options.runonspare
    profilers        = options.profilers.strip()
    cpu              = options.cpu.strip()
    bypasshlt        = options.bypasshlt
    cores            = options.cores
    logfile          = options.logfile
    prevrel          = options.previousrel
    outputdir        = options.outputdir

    if not logfile == None:
        logfile = os.path.abspath(logfile)
        logdir = os.path.dirname(logfile)
        if not os.path.exists(logdir):
            parser.error("Directory to output logfile does not exist")
            sys.exit()
        logfile = os.path.abspath(logfile)

    if "GEN,SIM" in stepOptions:
        print "WARNING: Please use GEN-SIM with a hypen not a \",\"!"

    isnumreg = re.compile("^-?[0-9]*$")
    found    = isnumreg.search(profilers)
    if not found :
        parser.error("profile codes option contains non-numbers")
        sys.exit()

    numetcomreg = re.compile("^[0-9,]*")

    if outputdir == "":
        outputdir = os.getcwd()
    else:
        outputdir = os.path.abspath(outputdir)

    if not os.path.isdir(outputdir):
        parser.error("%s is not a valid output directory" % outputdir)
        sys.exit()

    if not numetcomreg.search(cpu):
        parser.error("cpu option needs to be a comma separted list of ints or a single int")
        sys.exit()

    cpustr = cpu
    cpu = []
    if "," in cpustr:
        cpu = map(lambda x: int(x),cpustr.split(","))
    else:
        cpu = [ int(cpustr)  ]

    if not prevrel == "":
        prevrel = os.path.abspath(prevrel)
        if not os.path.exists(prevrel):
            print "ERROR: Previous release dir %s could not be found" % prevrel
            sys.exit()

    if quicktest:
        TimeSizeEvents = 1
        IgProfEvents = 1
        ValgrindEvents = 0
        cmsScimark = 1
        cmsScimarkLarge = 1

    if _unittest:
        _verbose = False
        if candleoption == "":
            candleoption = "MinBias"
        if stepOptions == "":
            stepOptions = "GEN-SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO"
        cmsScimark      = 0
        cmsScimarkLarge = 0
        ValgrindEvents  = 0
        IgProfEvents    = 0
        TimeSizeEvents  = 1

    if not cmsdriverOptions == "":
        cmsdriverOptions = "--cmsdriver=" + cmsdriverOptions        
        #Wrapping the options with "" for the cmsSimPyRelVal.pl until .py developed
        cmsdriverOptions= '"%s"' % (cmsdriverOptions)
        
    if not stepOptions == "":
        #Wrapping the options with "" for the cmsSimPyRelVal.pl until .py developed
        stepOptions='"--usersteps=%s"' % (stepOptions)
    
    isAllCandles = candleoption == ""
    candles = {}
    if isAllCandles:
        candles = Candles
    else:
        candles = candleoption.split(",")

    return (castordir       ,
            TimeSizeEvents  ,
            IgProfEvents    ,
            ValgrindEvents  ,
            cmsScimark      ,
            cmsScimarkLarge ,
            cmsdriverOptions,
            stepOptions     ,
            quicktest       ,
            profilers       ,
            cpu             ,
            cores           ,
            prevrel         ,
            isAllCandles    ,
            candles         ,
            bypasshlt       ,
            runonspare      ,
            outputdir       ,
            logfile         )

def usage():
    return __doc__

def runCmdSet(cmd):
    exitstat = None
    if len(cmd) <= 1:
        exitstat = runcmd(cmd)
        if _verbose:
            printFlush(cmd)
    else:
        for subcmd in cmd:
            if _verbose:
                printFlush(subcmd)
        exitstat = runcmd(" && ".join(cmd))
    if _verbose:
        printFlush(getDate())
    return exitstat

def printFlush(command):
    if _verbose:
        logh.write(command + "\n")
        logh.flush()

def runcmd(command):
    process  = os.popen(command)
    cmdout   = process.read()
    exitstat = process.close()
    if _verbose:
        logh.write(cmdout + "\n")
        logh.flush()
    return exitstat

def getDate():
    return time.ctime()

def printDate():
    logh.write(getDate() + "\n")

def getPrereqRoot(rootdir,rootfile):
    logh.write("WARNING: %s file required to run QCD profiling does not exist. Now running cmsDriver.py to get Required Minbias root file\n"   % (rootdir + "/" +rootfile))

    if not os.path.exists(rootdir):
        os.system("mkdir -p %s" % rootdir)
    if not _debug:
        cmd = "cd %s ; cmsDriver.py MinBias_cfi -s GEN,SIM -n %s >& ../minbias_for_pileup_generate.log" % (rootdir,str(10))
        log.write(cmd)
        os.system(cmd)
    if not os.path.exists(rootdir + "/" + rootfile):
        logh.write("ERROR: We can not run QCD profiling please create root file %s to run QCD profiling.\n" % (rootdir + "/" + rootfile))


def checkQcdConditions(candles,TimeSizeEvents,rootdir,rootfile):
    if TimeSizeEvents < MIN_REQ_TS_EVENTS :
        logh.write("WARNING: TimeSizeEvents is less than %s but QCD needs at least that to run. PILE-UP will be ignored\n" % MIN_REQ_TS_EVENTS)
        
        
    rootfilepath = rootdir + "/" + rootfile
    if not os.path.exists(rootfilepath):
        getPrereqRoot(rootdir,rootfile)
        if not os.path.exists(rootfilepath) and not _debug:
            logh.write("ERROR: Could not create or find a rootfile %s with enough TimeSize events for QCD exiting...\n" % rootfilepath)
            sys.exit()
    else:
        logh.write("%s Root file for QCD exists. Good!!!\n" % (rootdir + "/" + rootfile))
    return candles

def mkCandleDir(pfdir,candle,profiler):
    adir = os.path.join(pfdir,"%s_%s" % (candle,profiler))
    runcmd( "mkdir -p %s" % adir )
    if _verbose:
        printDate()
    #runCmdSet(cmd)
    return adir

def cpIgProfGenSim(dir,candle):
    cmds = ("cd %s" % dir,
            "cp -pR ../%s_IgProf/%s_GEN,SIM.root ."  % (candle,candle))
    runCmdSet(cmds)

def displayErrors(file):
    global ERRORS
    try:
        for line in open(file,"r"):
            if "cerr" in line:
                logh.write("ERROR: %s\n" % line)
                ERRORS += 1
    except OSError, detail:
        logh.write("WARNING: %s\n" % detail)
        ERRORS += 1        
    except IOError, detail:
        logh.write("WARNING: %s\n" % detail)
        ERRORS += 1
    

def valFilterReport(dir,cmsver):
    cmds = ("cd %s" % dir,
            "grep -v \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp" % (cmssw_version),
            "mv tmp SimulationCandles_%s.txt"                         % (cmssw_version))
    runCmdSet(cmds)

def benchmarks(cpu,pfdir,name,bencher,large=False):
    cmd = Commands[cpu][3]
    redirect = ""
    if large:
        redirect = " -large >& "    
    else:
        redirect = " >& "

    for i in range(bencher):
        command= cmd + redirect + os.path.join(pfdir,os.path.basename(name))        
        printFlush(command + " [%s/%s]" % (i+1,bencher))
        runcmd(command)
        logh.flush()

def runCmsReport(cpu,dir,cmsver,candle):
    cmd  = Commands[cpu][1]
    cmds = ("cd %s"                 % (dir),
            "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (cmd,cmsver,candle))
    exitstat = None
    if not _debug:
        exitstat = runCmdSet(cmds)
        
    if _unittest and (not exitstat == None):
        logh.write("ERROR: CMS Report returned a non-zero exit status \n")
        sys.exit()

def testCmsDriver(cpu,dir,cmsver,candle):
    cmsdrvreg = re.compile("^cmsDriver.py")
    cmd  = Commands[cpu][0]
    noExit = True
    stepreg = re.compile("--step=([^ ]*)")
    previousCmdOnline = ""
    for line in open(os.path.join(dir,"SimulationCandles_%s.txt" % (cmsver))):
        if (not line.lstrip().startswith("#")) and not (line.isspace() or len(line) == 0): 
            cmdonline  = line.split("@@@",1)[0]
            if cmsdrvreg.search(cmdonline) and not previousCmdOnline == cmdonline:
                stepbeingrun = "Unknown"
                matches = stepreg.search(cmdonline)
                if not matches == None:
                    stepbeingrun = matches.groups()[0]
                if "PILEUP" in cmdonline:
                    stepbeingrun += "_PILEUP"
                logh.write(cmdonline + "\n")
                cmds = ("cd %s"      % (dir),
                        "%s  >& ../cmsdriver_unit_test_%s_%s.log"    % (cmdonline,candle,stepbeingrun))
                if _dryrun:
                    logh.write(cmds + "\n")
                else:
                    out = runCmdSet(cmds)                    
                    if not out == None:
                        sig     = out >> 16    # Get the top 16 bits
                        xstatus = out & 0xffff # Mask out all bits except the first 16 
                        logh.write("FATAL ERROR: CMS Driver returned a non-zero exit status (which is %s) when running %s for candle %s. Signal interrupt was %s\n" % (xstatus,stepbeingrun,candle,sig))
                        sys.exit()
            previousCmdOnline = cmdonline
    

def runCmsInput(cpu,dir,numevents,candle,cmsdrvopts,stepopt,profiles,bypasshlt):

    bypass = ""
    if bypasshlt:
        bypass = "--bypass-hlt"
    cmd = Commands[cpu][2]
    cmds = ("cd %s"                    % (dir),
            "%s %s \"%s\" %s %s %s %s" % (cmd,
                                          numevents,
                                          candle,
                                          profiles,
                                          cmsdrvopts,
                                          stepopt,
                                          bypass))
    exitstat = runCmdSet(cmds)
    if _unittest and (not exitstat == None):
        logh.write("ERROR: CMS Report Input returned a non-zero exit status \n" )

def simpleGenReport(cpus,perfdir,NumEvents,candles,cmsdriverOptions,stepOptions,cmssw_version,Name,profilers,bypasshlt):
    valgrind = Name == "Valgrind"

    profCodes = {"TimeSize" : "0123",
                 "IgProf"   : "4567",
                 "Valgrind" : "89",
                 None       : "-1"} 

    profiles = profCodes[Name]
    if not profilers == "":
        profiles = profilers        

    for cpu in cpus:
        pfdir = perfdir
        if len(cpus) > 1:
            pfdir = os.path.join(perfdir,"cpu_%s" % cpu)
        for candle in candles:
            adir = mkCandleDir(pfdir,candle,Name)

            if valgrind:
                if candle == "SingleMuMinusPt10" : 
                    logh.write("Valgrind tests **GEN,SIM ONLY** on %s candle\n" % candle    )
                else:
                    logh.write("Valgrind tests **SKIPPING GEN,SIM** on %s candle\n" % candle)
                    cpIgProfGenSim(adir,candle)                

            if _unittest:
                # Run cmsDriver.py
                runCmsInput(cpu,adir,NumEvents,candle,cmsdriverOptions,stepOptions,profiles,bypasshlt)
                testCmsDriver(cpu,adir,cmssw_version,candle)
            else:
                runCmsInput(cpu,adir,NumEvents,candle,cmsdriverOptions,stepOptions,profiles,bypasshlt)            
                if valgrind:
                    valFilterReport(adir,cmssw_version)             
                runCmsReport(cpu,adir,cmssw_version,candle)
                proflogs = []
                if   Name == "TimeSize":
                    proflogs = [ "TimingReport" ]
                elif Name == "Valgrind":
                    pass
                elif Name == "IgProf":
                    pass

                for proflog in proflogs:
                    globpath = os.path.join(adir,"%s_*_%s.log" % (CandFname[candle],proflog))
                    logh.write("Looking for logs that match %s\n" % globpath)
                    logs     = glob.glob(globpath)
                    for log in logs:
                        logh.write("Found log %s\n" % log)
                        displayErrors(log)

def runPerfSuite(castordir        = _CASTOR_DIR,
                 perfsuitedir     = os.getcwd(),
                 TimeSizeEvents   = 100        ,
                 IgProfEvents     = 5          ,
                 ValgrindEvents   = 1          ,
                 cmsScimark       = 10         ,
                 cmsScimarkLarge  = 10         ,
                 cmsdriverOptions = ""         ,
                 stepOptions      = ""         ,
                 quicktest        = False      ,
                 profilers        = ""         ,
                 cpus             = [1]        ,
                 cores            = 4          ,
                 prevrel          = ""         ,
                 isAllCandles     = False      ,
                 candles          = Candles    ,
                 bypasshlt        = False      ,
                 runonspare       = True       ,
                 logfile          = os.path.join(os.getcwd(),"cmsPerfSuite.log")):

    global Commands, logh
    #Print a time stamp at the beginning:

    if not logfile == None:
        try:
            logh = open(logfile,"a")
        except (OSError, IOError), detail:
            logh.write(detail + "\n")

    try:        

        if not cmsdriverOptions == "":
            logh.write("Running cmsDriver.py with the special user defined options: %s\n" % cmsdriverOptions)
        
        if not stepOptions == "":
            logh.write("Running user defined steps only: %s\n" % stepOptions)

        if not len(candles) == len(Candles):
            logh.write("Running only %s candle, instead of the whole suite\n" % str(candles))
        
        logh.write("This machine ( %s ) is assumed to have %s cores, and the suite will be run on cpu %s\n" %(host,cores,cpus))
        path=os.path.abspath(".")
        logh.write("Performance Suite started running at %s on %s in directory %s, run by user %s\n" % (getDate(),host,path,user))
        showtags=os.popen4("showtags -r")[1].read()
        logh.write(showtags + "\n")

        #For the log:
        if _verbose:
            logh.write("The performance suite results tarball will be stored in CASTOR at %s\n" % castordir)
            logh.write("%s TimeSize events\n" % TimeSizeEvents)
            logh.write("%s IgProf events\n"   % IgProfEvents)
            logh.write("%s Valgrind events\n" % ValgrindEvents)
            logh.write("%s cmsScimark benchmarks before starting the tests\n"      % cmsScimark)
            logh.write("%s cmsScimarkLarge benchmarks before starting the tests\n" % cmsScimarkLarge)

        #Actual script actions!
        #Will have to fix the issue with the matplotlib pie-charts:
        #Used to source /afs/cern.ch/user/d/dpiparo/w0/perfreport2.1installation/share/perfreport/init_matplotlib.sh
        #Need an alternative in the release

        #Command Handling:


        if len(cpus) > 1:
            for cpu in cpus:
                cpupath = os.path.join(perfsuitedir,"cpu_%s" % cpu)
                if not os.path.exists(cpupath):
                    os.mkdir(cpupath)


        Commands = {}
        AllScripts = Scripts + AuxiliaryScripts

        for cpu in cpus:
            Commands[cpu] = []

        for script in AllScripts:
            which="which " + script

            #Logging the actual version of cmsDriver.py, cmsRelvalreport.py, cmsSimPyRelVal.pl
            whichstdout=os.popen4(which)[1].read()
            logh.write(whichstdout + "\n")
            if script in Scripts:
                for cpu in cpus:
                    command="taskset -c %s %s" % (cpu,script)
                    Commands[cpu].append(command)

        #First submit the cmsScimark benchmarks on the unused cores:
        scimark = ""
        scimarklarge = ""
        if not _unittest:
            for core in range(cores):
                if (not core in cpus) and runonspare:
                    logh.write("Submitting cmsScimarkLaunch.csh to run on core cpu "+str(core) + "\n")
                    command="taskset -c %s \"cd %s ; cmsScimarkLaunch.csh %s\" &" % (str(core), perfsuitedir, str(core))
                    logh.write(command + "\n")

                    #cmsScimarkLaunch.csh is an infinite loop to spawn cmsScimark2 on the other
                    #cpus so it makes no sense to try reading its stdout/err 
                    os.popen4(command)
        logh.flush()

        #dont do benchmarking if in debug mode... saves time
        benching = not _debug
        if benching and not _unittest:
            #Submit the cmsScimark benchmarks on the cpu where the suite will be run:        
            scimark      = open(os.path.join(perfsuitedir,"cmsScimark2.log")      ,"w")        
            scimarklarge = open(os.path.join(perfsuitedir,"cmsScimark2_large.log"),"w")
            if cmsScimark > 0:
                logh.write("Starting with %s cmsScimark on cpu%s\n"       % (cmsScimark,cpu))
                benchmarks(cpu,perfsuitedir,scimark.name,cmsScimark)

            if cmsScimarkLarge > 0:
                logh.write("Following with %s cmsScimarkLarge on cpu%s\n" % (cmsScimarkLarge,cpu))
                benchmarks(cpu,perfsuitedir,scimarklarge.name,cmsScimarkLarge)

        if not profilers == "":
            # which profile sets should we go into if custom profiles have been selected
            runTime     = reduce(lambda x,y: x or y, map(lambda x: x in profilers, ["0", "1", "2", "3"]))
            runIgProf   = reduce(lambda x,y: x or y, map(lambda x: x in profilers, ["4", "5", "6", "7"]))
            runValgrind = reduce(lambda x,y: x or y, map(lambda x: x in profilers, ["8", "9"]))
            if not runTime:
                TimeSizeEvents = 0
            if not runIgProf:
                IgProfEvents   = 0
            if not runValgrind:
                ValgrindEvents = 0

        qcdWillRun = (not isAllCandles) and "QCD_80_120" in candles 
        if qcdWillRun:
            candles = checkQcdConditions(candles,
                                         TimeSizeEvents,
                                         os.path.join(perfsuitedir,"%s_%s" % ("MinBias","TimeSize")),
                                         "%s_cfi_GEN_SIM.root" % "MinBias")

        #TimeSize tests:
        if TimeSizeEvents > 0:
            logh.write("Launching the TimeSize tests (TimingReport, TimeReport, SimpleMemoryCheck, EdmSize) with %s events each\n" % TimeSizeEvents)
            printDate()
            logh.flush()
            simpleGenReport(cpus,perfsuitedir,TimeSizeEvents,candles,cmsdriverOptions,stepOptions,cmssw_version,"TimeSize",profilers,bypasshlt)

        #IgProf tests:
        if IgProfEvents > 0:
            logh.write("Launching the IgProf tests (IgProfPerf, IgProfMemTotal, IgProfMemLive, IgProfMemAnalyse) with %s events each\n" % IgProfEvents)
            printDate()
            logh.flush()
            IgCandles = candles
            #By default run IgProf only on QCD_80_120 candle
            if isAllCandles:
                IgCandles = [ "QCD_80_120" ]
            simpleGenReport(cpus,perfsuitedir,IgProfEvents,IgCandles,cmsdriverOptions,stepOptions,cmssw_version,"IgProf",profilers,bypasshlt)

        #Valgrind tests:
        if ValgrindEvents > 0:
            logh.write("Launching the Valgrind tests (callgrind_FCE, memcheck) with %s events each\n" % ValgrindEvents)
            printDate()
            logh.flush()
            valCandles = candles

            if isAllCandles:
                cmds=[]
                #By default run Valgrind only on QCD_80_120, skipping SIM step since it would take forever (and do SIM step on SingleMu)
                valCandles = [ "QCD_80_120" ]

            #Besides always run, only once the GEN,SIM step on SingleMu:
            valCandles.append("SingleMuMinusPt10")
            #In the user-defined candles a different behavior: do Valgrind for all specified candles (usually it will only be 1)
            #usercandles=candleoption.split(",")
            simpleGenReport(cpus,perfsuitedir,ValgrindEvents,valCandles,cmsdriverOptions,stepOptions,cmssw_version,"Valgrind",profilers,bypasshlt)

        if benching and not _unittest:
            #Ending the performance suite with the cmsScimark benchmarks again:
            if cmsScimark > 0:
                logh.write("Ending with %s cmsScimark on cpu%s\n"         % (cmsScimark,cpu))
                benchmarks(cpu,perfsuitedir,scimark.name,cmsScimark)

            if cmsScimarkLarge > 0:
                logh.write("Following with %s cmsScimarkLarge on cpu%s\n" % (cmsScimarkLarge,cpu))
                benchmarks(cpu,perfsuitedir,scimarklarge.name,cmsScimarkLarge)

        #Stopping all cmsScimark jobs and analysing automatically the logfiles
        logh.write("Stopping all cmsScimark jobs\n")
        stopcmd = "cd %s ; %s" % (perfsuitedir,AuxiliaryScripts[2])
        printFlush(stopcmd)
        printFlush(os.popen4(stopcmd)[1].read())

        if not prevrel == "":
            crr.regressReports(prevrel,os.path.abspath(perfsuitedir),oldRelName = getVerFromLog(prevrel),newRelName=cmssw_version)

        #Create a tarball of the work directory
        TarFile = "%s_%s_%s.tar" % (cmssw_version, host, user)
        AbsTarFile = os.path.join(perfsuitedir,TarFile)
        tarcmd  = "tar -cvf %s %s; gzip %s" % (AbsTarFile,os.path.join(perfsuitedir,"*"),AbsTarFile)
        printFlush(tarcmd)
        printFlush(os.popen4(tarcmd)[1].read())

        #Archive it on CASTOR
        castorcmd="rfcp %s.gz %s.gz" % (AbsTarFile,os.path.join(castordir,TarFile))

        printFlush(castorcmd)
        printFlush(os.popen4(castorcmd)[1].read())

        #End of script actions!

        #Print a time stamp at the end:
        date=time.ctime(time.time())
        logh.write("Performance Suite finished running at %s on %s in directory %s\n" % (date,host,path))
        if ERRORS == 0:
            logh.write("There were no errors detected in any of the log files!\n")
        else:
            logh.write("ERROR: There were %s errors detected in the log files, please revise!\n" % ERRORS)
    except:
        logh.flush()
        if not logh.isatty():
            logh.close()
        raise

if __name__ == "__main__":
    #Let's check the command line arguments
    (castordir       ,
     TimeSizeEvents  ,
     IgProfEvents    ,
     ValgrindEvents  ,
     cmsScimark      ,
     cmsScimarkLarge ,
     cmsdriverOptions,
     stepOptions     ,
     quicktest       ,
     profilers       ,
     cpus            ,
     cores           ,
     prevrel         ,
     isAllCandles    ,
     candles         ,
     bypasshlt       ,
     runonspare      ,
     outputdir       ,
     logfile         ) = optionParse()
   
    runPerfSuite(castordir        = castordir       ,
                 perfsuitedir     = outputdir       ,
                 TimeSizeEvents   = TimeSizeEvents  ,
                 IgProfEvents     = IgProfEvents    ,
                 ValgrindEvents   = ValgrindEvents  ,
                 cmsScimark       = cmsScimark      ,
                 cmsScimarkLarge  = cmsScimarkLarge ,
                 cmsdriverOptions = cmsdriverOptions,
                 stepOptions      = stepOptions     ,
                 quicktest        = quicktest       ,
                 profilers        = profilers       ,
                 cpus             = cpus            ,
                 cores            = cores           ,
                 prevrel          = prevrel         ,
                 isAllCandles     = isAllCandles    ,
                 candles          = candles         ,
                 bypasshlt        = bypasshlt       ,
                 runonspare       = runonspare      ,
                 logfile          = logfile         )     

