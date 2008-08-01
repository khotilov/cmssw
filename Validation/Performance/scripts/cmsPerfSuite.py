#!/usr/bin/python
'''
Usage: ./cmsPerfSuite.py [options]
       
Options:
  -o ..., --output=...   specify the wanted CASTOR directory where to store the results tarball
  -t ..., --timesize=... specify the number of events for the TimeSize tests
  -i ..., --igprof=...   specify the number of events for the IgProf tests
  -v ..., --valgrind=... specify the number of events for the Valgrind tests
  --cmsScimark=...       specify the number of times the cmsScimark benchmark is run before and after the performance suite on cpu1
  --cmsScimarkLarge=...  specify the number of times the cmsScimarkLarge benchmark is run before and after the performance suite on cpu1
  --cmsdriver=...        specify special options to use with the cmsDriver.py commands (designed for integration build use)
  --step=...             specify the processing steps intended (instead of the default ones)
  --candle=...           specify the candle(s) to run (instead of all 7 default candles)
  --cpu=...              specify the core on which to run the performance suite
  --cores=...            specify the number of cores of the machine (can be used with 0 to stop cmsScimark from running on the other cores)
  -h, --help           show this help
  -d                   show debugging information

Examples:
./cmsPerfSuite.py
(this will run with the default options)
OR
./cmsPerfSuite.py -o "/castor/cern.ch/user/y/yourusername/yourdirectory/"
(this will archive the results in a tarball on /castor/cern.ch/user/y/yourusername/yourdirectory/)
OR
./cmsPerfSuite.py -t 5 -i 2 -v 1
(this will run the suite with 5 events for TimeSize tests, 2 for IgProf tests, 0 for Valgrind tests)
OR
./cmsPerfSuite.py -t 200 --candle QCD_80_120 --cmsdriver="--conditions FakeConditions"
(this will run the performance tests only on candle QCD_80_120, running 200 TimeSize evts, default IgProf and Valgrind evts. It will also add the option "--conditions FakeConditions" to all cmsDriver.py commands executed by the suite)
OR
./cmsPerfSuite.py -t 200 --candle QCD_80_120 --cmsdriver="--conditions=FakeConditions --eventcontent=FEVTDEBUGHLT" --step=GEN-SIM,DIGI
(this will run the performance tests only on candle QCD_80_120, running 200 TimeSize evts, default IgProf and Valgrind evts. It will also add the option "--conditions=FakeConditions" and the option "--eventcontent=FEVTDEBUGHLT" to all cmsDriver.py commands executed by the suite. In addition it will run only 2 cmsDriver.py "steps": "GEN,SIM" and "DIGI". Note the syntax GEN-SIM for combined cmsDriver.py steps)

Legal entries for individual candles (--candle option):
HiggsZZ4LM200
MinBias
SingleElectronE1000
SingleMuMinusPt10
SinglePiMinusE1000
TTbar
QCD_80_120
'''
import os
import time
import getopt
import sys

MIN_REQ_TS_EVENTS = 8

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
#Some defaults:


#Options handling


def usage():
    print __doc__

def runCmdSet(cmd):
    for subcmd in cmd:
        printFlush(subcmd)
    runcmd(";".join(cmd))
    printFlush(getDate())

def printFlush(command):
    print command
    sys.stdout.flush()

def runcmd(command):
    cmdout=os.popen4(command)[1].read()
    print cmdout

def benchmarks(cmd,redirect,name,bencher):
    numofbenchs = int(bencher)
    for i in range(numofbenchs):
        command= cmd + redirect + name
        printFlush(command+" [%s/%s]"%(i+1,numofbenchs))
        runcmd(command)
        sys.stdout.flush()

def getDate():
    return time.ctime()

def printDate():
    print getDate()

def getPrereqRoot(rootdir,rootfile):
    os.system("cd %s ; cmsDriver.py MinBias_cfi -s GEN,SIM -n 10" % (rootdir))
    print "ERROR: %s file required to run QCD profiling does not exist. We can not run QCD profiling please create root file" % (mbrootfile)
    print "       to run QCD profiling."
    print "       Running cmsDriver.py to get Required MinbiasEvents"


def checkQcdConditions(isAllCandles,usercandles,AllCandles,TimeSizeEvents,rootdir,rootfile):
    if TimeSizeEvents < MIN_REQ_TS_EVENTS :
        print "WARNING: TimeSizeEvents is less than 8 but QCD needs at least that to run. Setting TimeSizeEvents to 8"
        if isAllCandles:
            AllCandles.popitem("QCD_80_120")
        else :
            usercandles.pop("QCD_80_120")
        
    rootfilepath = rootdir + "/" + rootfile
    if not os.path.exists(rootfilepath):
        getPrereqRoot(rootdir,rootfile)
        if not os.path.exists(rootfilepath):
            print "ERROR: Could not create root with enough TimeSize events for QCD exiting..."
            sys.exit()
    return (usercandles,AllCandles)

def main(argv):
    #Some default values:
    castordir = "/castor/cern.ch/cms/store/relval/performance/"
    TimeSizeEvents   = "100"
    IgProfEvents     = "5"
    ValgrindEvents   = "1"
    cmsScimark       = "10"
    cmsScimarkLarge  = "10"
    cmsdriverOptions = ""
    stepOptions      = ""
    candleoption     = ""
    #Number of cpu cores on the machine
    cores=4
    #Cpu core on which the suite is run:
    cpu=1
    #Let's check the command line arguments
    try:
        opts, args = getopt.getopt(argv, "o:t:i:v:hd", ["output=","timesize=","igprof=","valgrind=","cmsScimark=","cmsScimarkLarge=","cmsdriver=","step=","candle=","cpu=","cores=","help"])
    except getopt.GetoptError:
        print "This argument option is not accepted"
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt == '-d':
            global _debug
            _debug = 1
        elif opt in ("-o", "--output"):
            castordir= arg
        elif opt in ("-t","--timesize"):
            TimeSizeEvents = arg
        elif opt in ("-i", "--igprof"):
            IgProfEvents = arg
        elif opt in ("-v", "--valgrind"):
            ValgrindEvents = arg
        elif opt == "--cmsScimark":
            cmsScimark = arg
        elif opt == "--cmsScimarkLarge":
            cmsScimarkLarge = arg
        elif opt in ("-c","--cmsdriver"):
            cmsdriverOptions= arg
        elif opt == "--step":
            stepOptions=arg
        elif opt == "--candle":
            candleoption=arg
        elif opt == "--cpu":
            cpu=arg
        elif opt == "--cores":
            cores=arg
    #Case with no arguments (using defaults)
    if opts == []:
        print "No arguments given, so DEFAULT test will be run:"
        
    #Print a time stamp at the beginning:

    path=os.path.abspath(".")
    print "Performance Suite started running at %s on %s in directory %s, run by user %s" % (getDate(),host,path,user)
    showtags=os.popen4("showtags -r")[1].read()
    print showtags
    
    #For the log:
    print "The performance suite results tarball will be stored in CASTOR at %s" % castordir
    print "%s TimeSize events" % TimeSizeEvents
    print "%s IgProf events"   % IgProfEvents
    print "%s Valgrind events" % ValgrindEvents
    print "%s cmsScimark benchmarks before starting the tests"      % cmsScimark
    print "%s cmsScimarkLarge benchmarks before starting the tests" % cmsScimarkLarge
    if cmsdriverOptions != "":
        print "Running cmsDriver.py with the special user defined options: %s" % cmsdriverOptions
        
        #Wrapping the options with "" for the cmsSimPyRelVal.pl until .py developed
        cmsdriverOptions= '"%s"' % (cmsdriverOptions)
    if stepOptions !="":
        print "Running user defined steps only: %s" % stepOptions
        
        #Wrapping the options with "" for the cmsSimPyRelVal.pl until .py developed
        stepOptions='"--usersteps=%s"' % (stepOptions)
    if candleoption !="":
        print "Running only %s candle, instead of the whole suite" % candleoption
    print "This machine ( %s ) is assumed to have %s cores, and the suite will be run on cpu %s" %(host,cores,cpu)
    
    #Actual script actions!
    #Will have to fix the issue with the matplotlib pie-charts:
    #Used to source /afs/cern.ch/user/d/dpiparo/w0/perfreport2.1installation/share/perfreport/init_matplotlib.sh
    #Need an alternative in the release

    #Command Handling:
    Commands=[]
    AuxiliaryCommands=[]
    AllScripts=Scripts+AuxiliaryScripts
    for script in AllScripts:
        which="which "+script
        
        #Logging the actual version of cmsDriver.py, cmsRelvalreport.py, cmsSimPyRelVal.pl
        whichstdout=os.popen4(which)[1].read()
        print whichstdout
        if script in Scripts:
            command="taskset -c "+str(cpu)+" "+script
            Commands.append(command)
        elif script == "cmsScimarkLaunch.csh":
            for core in range(int(cores)):
                if core != int(cpu):
                    command="taskset -c "+str(core)+" "+script+" "+str(core)
                    AuxiliaryCommands.append(command)
        else:
            command=script
            AuxiliaryCommands.append(command)
            
    #print Commands
    #print AuxiliaryCommands
    sys.stdout.flush()
    
    #First submit the cmsScimark benchmarks on the unused cores:
    for core in range(int(cores)):
        if core != int(cpu):
            print "Submitting cmsScimarkLaunch.csh to run on core cpu"+str(core)
            command="taskset -c "+str(core)+" cmsScimarkLaunch.csh "+str(core)+"&"
            print command
            
            #cmsScimarkLaunch.csh is an infinite loop to spawn cmsScimark2 on the other
            #cpus so it makes no sense to try reading its stdout/err 
            os.popen4(command)
            
    #Submit the cmsScimark benchmarks on the cpu where the suite will be run:
    scimark      = open("cmsScimark2.log"      ,"w")
    scimarklarge = open("cmsScimark2_large.log","w")

    #Test mode... development use only... we need to add an option for this
    # but getopt is rubbish... we should move this to optparse.
    testmode = False
    #dont do benchmarking if in test mode... saves time
    benching = not testmode
    if benching:
        print "Starting with %s cmsScimark on cpu%s"%(cmsScimark,cpu)
        benchmarks(Commands[3]," >& ",scimark.name,cmsScimark)
    
        print "Following with %s cmsScimarkLarge on cpu%s"%(cmsScimarkLarge,cpu)
        benchmarks(Commands[3]," -large >& ",scimarklarge.name,cmsScimarkLarge)
        
    #Here the real performance suite starts
    #List of Candles
    Candles={"HiggsZZ4LM200"      : "HZZLLLL",
             "MinBias"            : "MINBIAS",
             "SingleElectronE1000": "E -e 1000",
             "SingleMuMinusPt10"  : "MU- -e pt10",
             "SinglePiMinusE1000" : "PI- -e 1000",
             "TTbar"              : "TTBAR",
             "QCD_80_120"         : "QCD -e 80_120"
             }
    AllCandles=Candles.keys()
    
    #Sort the candles to make sure MinBias is executed before QCD_80_120, otherwise DIGI PILEUP would not find its MinBias root files
    AllCandles.sort()

    isAllCandles = candleoption == ""

    if not isAllCandles:
        usercandles=candleoption.split(",")

    qcdWillRun = isAllCandles or ((not isAllCandles) and "QCD_80_120" in usercandles )

    #TimeSize tests:
    if int(TimeSizeEvents)>0:

        print "Launching the TimeSize tests (TimingReport, TimeReport, SimpleMemoryCheck, EdmSize) with %s events each" % TimeSizeEvents
        printDate()

        if qcdWillRun:
            (usercandles,
             AllCandles)     = checkQcdConditions(isAllCandles,
                                                  usercandles,
                                                  AllCandles,
                                                  TimeSizeEvents,
                                                  "./MinBias_TimeSize",
                                                  "MinBias_GEN,SIM.root")

        if isAllCandles:
            cmds=[]
            sys.stdout.flush()
            for candle in AllCandles:
                cmd = ("mkdir %s_TimeSize"                                                % (candle),
                       "cd %s_TimeSize"                                                   % (candle),
                       "%s %s \"%s\" 0123 %s %s"                                          % (Commands[2],
                                                                                             TimeSizeEvents,
                                                                                             Candles[candle],
                                                                                             cmsdriverOptions,
                                                                                             stepOptions),
                       "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
                runCmdSet(cmd)

        else:
            for candle in usercandles:
                cmd = ("mkdir %s_TimeSize"                                                % (candle),
                       "cd %s_TimeSize"                                                   % (candle),
                       "%s %s \"%s\" 0123 %s %s"                                          % (Commands[2],
                                                                                             TimeSizeEvents,
                                                                                             Candles[candle],
                                                                                             cmsdriverOptions,
                                                                                             stepOptions),
                       "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
                runCmdSet(cmd)

    #IgProf tests:
    if int(IgProfEvents)>0:
        print "Launching the IgProf tests (IgProfPerf, IgProfMemTotal, IgProfMemLive, IgProfMemAnalyse) with %s events each" % IgProfEvents
        printDate()

        if qcdWillRun:
            (usercandles,
             AllCandles) = checkQcdConditions(isAllCandles,
                                                  usercandles,
                                                  AllCandles,
                                                  TimeSizeEvents,
                                                  "./MinBias_TimeSize",
                                                  "MinBias_GEN,SIM.root")           
        
        if isAllCandles:
            cmds=[]
            sys.stdout.flush()
            
            #By default run IgProf only on QCD_80_120 candle
            candle = "QCD_80_120"
            cmd = ("mkdir %s_IgProf"                                                  % (candle),
                   "cd %s_IgProf"                                                     % (candle),
                   "%s %s \"%s\" 4567 %s %s"                                          % (Commands[2],
                                                                                         IgProfEvents,
                                                                                         Candles[candle],
                                                                                         cmsdriverOptions,
                                                                                         stepOptions),
                   "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))

            runCmdSet(cmd)
        else:
            
            #In the user-defined candles a different behavior: do IgProf for all specified candles (usually it will only be 1)
            usercandles=candleoption.split(",")
            for candle in usercandles:
                cmd = ("mkdir %s_IgProf"                                                  % (candle),
                       "cd %s_IgProf"                                                     % (candle),
                       "%s %s \"%s\" 4567 %s %s"                                          % (Commands[2],
                                                                                             IgProfEvents,
                                                                                             Candles[candle],
                                                                                             cmsdriverOptions,
                                                                                             stepOptions),
                       "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
                runtCmdSet(cmd)
    #Valgrind tests:
    if int(ValgrindEvents)>0:
        print "Launching the Valgrind tests (callgrind_FCE, memcheck) with %s events each" % ValgrindEvents
        printDate()

        if qcdWillRun:
            (usercandles,
             AllCandles) = checkQcdConditions(isAllCandles,
                                              usercandles,
                                              AllCandles,
                                              TimeSizeEvents,
                                              "./MinBias_Valgrind",
                                              "MinBias_GEN,SIM.root")                                                
        
        if isAllCandles:
            cmds=[]
            sys.stdout.flush()
            
            #By default run Valgrind only on QCD_80_120, skipping SIM step since it would take forever (and do SIM step on SingleMu)
            candle = "QCD_80_120"
            print "Valgrind tests **SKIPPING GEN,SIM** on %s candle" % candle
            cmd = ("mkdir %s_Valgrind"                                                % (candle),
                   "cd %s_Valgrind"                                                   % (candle),
                   "cp -pR ../%s_IgProf/%s_GEN,SIM.root ."                            % (candle,candle),
                   "%s %s \"%s\" 89 %s %s"                                            % (Commands[2],
                                                                                         ValgrindEvents,
                                                                                         Candles[candle],
                                                                                         cmsdriverOptions,
                                                                                         stepOptions),
                   "grep -v \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp"          % (cmssw_version),
                   "mv tmp SimulationCandles_%s.txt"                                  % (cmssw_version),
                   "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
            runCmdSet(cmd)
            
            #By default run Valgring GEN,SIM profiling only on SingleMu (fastest) candle
            candle = "SingleMuMinusPt10"
            print "Valgrind tests **GEN,SIM ONLY** on %s candle" % candle
            cmd = ("mkdir %s_Valgrind"                                                % (candle),
                   "cd %s_Valgrind"                                                   % (candle),
                   "%s %s \"%s\" 89 %s %s"                                            % (Commands[2],
                                                                                         ValgrindEvents,
                                                                                         Candles[candle],
                                                                                         cmsdriverOptions,
                                                                                         stepOptions),
                   "grep \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp"             % (cmssw_version),
                   "mv tmp SimulationCandles_%s.txt"                                  % (cmssw_version),
                   "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
            
            runCmdSet(cmd)
        else:
            
            #In the user-defined candles a different behavior: do Valgrind for all specified candles (usually it will only be 1)
            usercandles=candleoption.split(",")
            for candle in usercandles:
                print "Valgrind tests **SKIPPING GEN,SIM** on %s candle" % candle
                cmd = ("mkdir %s_Valgrind"                                                % (candle),
                       "cd %s_Valgrind"                                                   % (candle),
                       "cp -pR ../%s_IgProf/%s_GEN,SIM.root ."                            % (candle,candle) ,
                       "%s %s \"%s\" 89 %s %s"                                            % (Commands[2],
                                                                                             ValgrindEvents,
                                                                                             Candles[candle],
                                                                                             cmsdriverOptions,
                                                                                             stepOptions),
                       "grep -v \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp"          % (cmssw_version),
                       "mv tmp SimulationCandles_%s.txt"                                  % (cmssw_version),
                       "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle))
 
                runCmdSet(cmd)
                
            #Besides always run, only once the GEN,SIM step on SingleMu:
            candle = "SingleMuMinusPt10"
            print "Valgrind tests **GEN,SIM ONLY** on %s candle" % candle
            cmd = ("mkdir %s_Valgrind"                                                % (candle),
                   "cd %s_Valgrind"                                                   % (candle),
                   "%s %s \"%s\" 89 %s %s"                                            % (candle,
                                                                                         Commands[2],
                                                                                         ValgrindEvents,
                                                                                         Candles[candle],
                                                                                         cmsdriverOptions),
                   "grep \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp"             % (stepOptions),
                   "mv tmp SimulationCandles_%s.txt"                                  % (cmssw_version),
                   "%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log" % (Commands[1],cmssw_version,candle)
                   )
#cmd = ("mkdir %s_Valgrind;cd %s_Valgrind;%s %s \"%s\" 89 %s %s;grep \"step=GEN,SIM\" SimulationCandles_%s.txt > tmp;mv tmp SimulationCandles_%s.txt;%s -i SimulationCandles_%s.txt -t perfreport_tmp -R -P >& %s.log% (candle,candle,Commands[2],ValgrindEvents,Candles[candle],cmsdriverOptions,stepOptions,cmssw_version,cmssw_version,Commands[1],cmssw_version,candle))

            runCmdSet(cmd)
        

    if benching:
    #Ending the performance suite with the cmsScimark benchmarks again: 
        print "Ending with %s cmsScimark on cpu%s"%(cmsScimark,cpu)
        benchmarks(" >& ",scimark.name,cmsScimark)
    
        print "Following with %s cmsScimarkLarge on cpu%s"%(cmsScimarkLarge,cpu)
        benchmarks(" -large >& ",scimarklarge.name,cmsScimarkLarge)
    
    #Stopping all cmsScimark jobs and analysing automatically the logfiles
    print "Stopping all cmsScimark jobs"
    printFlush(AuxiliaryScripts[2])

    printFlush(os.popen4(AuxiliaryScripts[2])[1].read())

    #Create a tarball of the work directory
    TarFile = cmssw_version + "_"     +     host    + "_"     + user + ".tar"
    tarcmd  = "tar -cvf "   + TarFile + " *; gzip " + TarFile
    printFlush(tarcmd)
    printFlush(os.popen4(tarcmd)[1].read())
    
    #Archive it on CASTOR
    castorcmd="rfcp "+TarFile+".gz "+castordir+TarFile+".gz"
    printFlush(castorcmd)
    printFlush(os.popen4(castorcmd)[1].read())

    #End of script actions!

    #Print a time stamp at the end:
    date=time.ctime(time.time())
    print "Performance Suite finished running at %s on %s in directory %s" % (date,host,path)
    
if __name__ == "__main__":
    main(sys.argv[1:])

