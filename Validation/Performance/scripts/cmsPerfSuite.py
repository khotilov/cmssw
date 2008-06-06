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
  --candle=...           specify the candle(s) to run (instead of all 7 default candles)
  --cpu=...        spedify the core on which to run the performance suite
  --cores=...            specify the number of cores of the machine
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
'''
import os
#Get some environment variables to use
cmssw_base=os.environ["CMSSW_BASE"]
cmssw_release_base=os.environ["CMSSW_RELEASE_BASE"]
cmssw_version=os.environ["CMSSW_VERSION"]
host=os.environ["HOST"]
user=os.environ["USER"]

#Scripts used by the suite:
Scripts=["cmsDriver.py","cmsRelvalreport.py","cmsSimPyRelVal.pl","cmsScimark2"]
AuxiliaryScripts=["cmsScimarkLaunch.csh","cmsScimarkParser.py","cmsScimarkStop.pl"]
#Some defaults:


#Options handling
import getopt
import sys

def usage():
    print __doc__

def main(argv):
    #Some default values:
    castordir = "/castor/cern.ch/cms/store/relval/performance/"
    TimeSizeEvents = "100"
    IgProfEvents = "5"
    ValgrindEvents = "1"
    cmsScimark = "10"
    cmsScimarkLarge = "10"
    cmsdriverOptions = ""
    candleoption=""
    #Number of cpu cores on the machine
    cores=4
    #Cpu core on which the suite is run:
    cpu=1
    #Let's check the command line arguments
    try:
        opts, args = getopt.getopt(argv, "o:t:i:v:hd", ["output=","timesize=","igprof=","valgrind=","cmsScimark=","cmsScimarkLarge=","cmsdriver=","candle=","cpu=","cores=","help"])
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
    import time
    date=time.ctime()
    path=os.path.abspath(".")
    print "Performance Suite started running at %s on %s in directory %s, run by user %s" % (date,host,path,user)
    showtags=os.popen4("showtags -r")[1].read()
    print showtags
    #For the log:
    print "The performance suite results tarball will be stored in CASTOR at %s" % castordir
    print "%s TimeSize events" % TimeSizeEvents
    print "%s IgProf events" % IgProfEvents
    print "%s Valgrind events" % ValgrindEvents
    print "%s cmsScimark benchmarks before starting the tests" % cmsScimark
    print "%s cmsScimarkLarge benchmarks before starting the tests" % cmsScimarkLarge
    if cmsdriverOptions != "":
        print "Running cmsDriver.py with the special user defined options: %s" % cmsdriverOptions
        #Wrapping the options with "" for the cmsSimPyRelVal.pl until .py developed
        cmsdriverOptions='"'+cmsdriverOptions+'"'
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
                if core != cpu:
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
        if core != cpu:
            print "Submitting cmsScimarkLaunch.csh to run on core cpu"+str(core)
            command="taskset -c "+str(core)+" cmsScimarkLaunch.csh "+str(core)+"&"
            print command
            #cmsScimarkLaunch.csh is an infinite loop to spawn cmsScimark2 on the other
            #cpus so it makes no sense to try reading its stdout/err 
            os.popen4(command)
    #Submit the cmsScimark benchmarks on the cpu where the suite will be run:
    scimark=open("cmsScimark2.log","w")
    scimarklarge=open("cmsScimark2_large.log","w")
    print "Starting with %s cmsScimark on cpu%s"%(cmsScimark,cpu)
    for i in range(int(cmsScimark)):
        command= Commands[3]+" >& "+scimark.name
        print command+" [%s/%s]"%(i+1,int(cmsScimark))
        cmsScimarkstdout=os.popen4(command)[1].read()
        print cmsScimarkstdout
    print "Following with %s cmsScimarkLarge on cpu%s"%(cmsScimarkLarge,cpu)
    for i in range(int(cmsScimarkLarge)):
        command= Commands[3]+" -large >& "+scimarklarge.name
        print command+" [%s/%s]"%(i+1,int(cmsScimarkLarge))
        cmsScimarkLargestdout=os.popen2(command)[1].read()
        print cmsScimarkLargestdout
    #Here the real performance suite starts
    #List of Candles
    Candles={"HiggsZZ4LM190":"HZZLLLL",
             "MinBias":"MINBIAS",
             "SingleElectronE1000":"E -e 1000",
             "SingleMuMinusPt10":"MU- -e pt10",
             "SinglePiMinusE1000":"PI- -e 1000",
             "TTbar":"TTBAR",
             "QCD_80_120":"QCD -e 80_120"
             }
    AllCandles=Candles.keys()
    #Sort the candles to make sure MinBias is executed before QCD_80_120, otherwise DIGI PILEUP would not find its MinBias root files
    AllCandles.sort()
    #TimeSize tests:
    if int(TimeSizeEvents)>0:
        print "Launching the TimeSize tests (TimingReport, TimeReport, SimpleMemoryCheck, EdmSize) with %s events each" % TimeSizeEvents
        if candleoption == "":
            cmds=[]
            sys.stdout.flush()
            for candle in AllCandles:
                cmd = 'mkdir '+candle+'_TimeSize;cd '+candle+'_TimeSize;'+Commands[2]+' '+TimeSizeEvents+' "'+Candles[candle]+'" 0123 '+cmsdriverOptions+';'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
                for subcmd in cmd.split(";"):
                    print subcmd
                    date=time.ctime()
                    print date
                    sys.stdout.flush()
                TimeSizecmdstdout=os.popen4(cmd)[1].read()
                print TimeSizecmdstdout
                date=time.ctime()
                print date
                sys.stdout.flush()
        else:
            usercandles=candleoption.split(",")
            for candle in usercandles:
                cmd = 'mkdir '+candle+'_TimeSize;cd '+candle+'_TimeSize;'+Commands[2]+' '+TimeSizeEvents+' "'+Candles[candle]+'" 0123 '+cmsdriverOptions+';'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
                for subcmd in cmd.split(";"):
                    print subcmd
                    date=time.ctime()
                    print date
                    sys.stdout.flush()
                TimeSizecmdstdout=os.popen4(cmd)[1].read()
                print TimeSizecmdstdout
                date=time.ctime()
                print date
                sys.stdout.flush()
    #IgProf tests:
    if int(IgProfEvents)>0:
        print "Launching the IgProf tests (IgProfPerf, IgProfMemTotal, IgProfMemLive, IgProfMemAnalyse) with %s events each" % IgProfEvents
        if candleoption == "":
            cmds=[]
            sys.stdout.flush()
            #By default run IgProf only on QCD_80_120 candle
            candle = "QCD_80_120"
            cmd = 'mkdir '+candle+'_IgProf;cd '+candle+'_IgProf;'+Commands[2]+' '+IgProfEvents+' "'+Candles[candle]+'" 4567 '+cmsdriverOptions+';'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
            for subcmd in cmd.split(";"):
                print subcmd
                date=time.ctime()
                print date
                sys.stdout.flush()
            IgProfstdout=os.popen4(cmd)[1].read()
            print IgProfstdout
            date=time.ctime()
            print date
            sys.stdout.flush()
        else:
            #In the user-defined candles a different behavior: do IgProf for all specified candles (usually it will only be 1)
            usercandles=candleoption.split(",")
            for candle in usercandles:
                cmd = 'mkdir '+candle+'_IgProf;cd '+candle+'_IgProf;'+Commands[2]+' '+IgProfEvents+' "'+Candles[candle]+'" 4567 '+cmsdriverOptions+';'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
                for subcmd in cmd.split(";"):
                    print subcmd
                    date=time.ctime()
                    print date
                    sys.stdout.flush()
                IgProfstdout=os.popen4(cmd)[1].read()
                print IgProfstdout
                date=time.ctime()
                print date
                sys.stdout.flush()
    #Valgrind tests:
    if int(ValgrindEvents)>0:
        print "Launching the Valgrind tests (callgrind_FCE, memcheck) with %s events each" % ValgrindEvents
        if candleoption == "":
            cmds=[]
            sys.stdout.flush()
            #By default run Valgrind only on QCD_80_120, skipping SIM step since it would take forever (and do SIM step on SingleMu)
            candle = "QCD_80_120"
            print "Valgrind tests **SKIPPING GEN,SIM** on %s candle" % candle
            cmd = 'mkdir '+candle+'_Valgrind;cd '+candle+'_Valgrind;cp -pR ../'+candle+'_IgProf/'+candle+'_GEN,SIM.root .;'+Commands[2]+' '+ValgrindEvents+' "'+Candles[candle]+'" 89 '+cmsdriverOptions+';grep -v "step=GEN,SIM" SimulationCandles_'+cmssw_version+'.txt > tmp;mv tmp SimulationCandles_'+cmssw_version+'.txt;'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
            for subcmd in cmd.split(";"):
                print subcmd
                date=time.ctime()
                print date
                sys.stdout.flush()
            Valgrindstdout=os.popen4(cmd)[1].read()
            print Valgrindstdout
            date=time.ctime()
            print date
            sys.stdout.flush()
            #By default run Valgring GEN,SIM profiling only on SingleMu (fastest) candle
            candle = "SingleMuMinusPt10"
            print "Valgrind tests **GEN,SIM ONLY** on %s candle" % candle
            cmd = 'mkdir '+candle+'_Valgrind;cd '+candle+'_Valgrind;'+Commands[2]+' '+ValgrindEvents+' "'+Candles[candle]+'" 89 '+cmsdriverOptions+';grep "step=GEN,SIM" SimulationCandles_'+cmssw_version+'.txt > tmp;mv tmp SimulationCandles_'+cmssw_version+'.txt;'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
            for subcmd in cmd.split(";"):
                print subcmd
                date=time.ctime()
                print date
                sys.stdout.flush()
            Valgrindstdout=os.popen4(cmd)[1].read()
            print Valgrindstdout
            date=time.ctime()
            print date
            sys.stdout.flush()
        else:
            #In the user-defined candles a different behavior: do Valgrind for all specified candles (usually it will only be 1)
            usercandles=candleoption.split(",")
            for candle in usercandles:
                print "Valgrind tests **SKIPPING GEN,SIM** on %s candle" % candle
                cmd = 'mkdir '+candle+'_Valgrind;cd '+candle+'_Valgrind;cp -pR ../'+candle+'_IgProf/'+candle+'_GEN,SIM.root .;'+Commands[2]+' '+ValgrindEvents+' "'+Candles[candle]+'" 89 '+cmsdriverOptions+';grep -v "step=GEN,SIM" SimulationCandles_'+cmssw_version+'.txt > tmp;mv tmp SimulationCandles_'+cmssw_version+'.txt;'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
                for subcmd in cmd.split(";"):
                    print subcmd
                    date=time.ctime()
                    print date
                    sys.stdout.flush()
                Valgrindstdout=os.popen4(cmd)[1].read()
                print Valgrindstdout
                date=time.ctime()
                print date
                sys.stdout.flush()
            #Besides always run, only once the GEN,SIM step on SingleMu:
            candle = "SingleMuMinusPt10"
            print "Valgrind tests **GEN,SIM ONLY** on %s candle" % candle
            cmd = 'mkdir '+candle+'_Valgrind;cd '+candle+'_Valgrind;'+Commands[2]+' '+ValgrindEvents+' "'+Candles[candle]+'" 89 '+cmsdriverOptions+';grep "step=GEN,SIM" SimulationCandles_'+cmssw_version+'.txt > tmp;mv tmp SimulationCandles_'+cmssw_version+'.txt;'+Commands[1]+' -i SimulationCandles_'+cmssw_version+'.txt -t perfreport_tmp -R -P >& '+candle+'.log'
            for subcmd in cmd.split(";"):
                print subcmd
                date=time.ctime()
                print date
                sys.stdout.flush()
            Valgrindstdout=os.popen4(cmd)[1].read()
            print Valgrindstdout
            date=time.ctime()
            print date
            sys.stdout.flush()
        

    #Ending the performance suite with the cmsScimark benchmarks again:
    print "Ending with %s cmsScimark on cpu%s"%(cmsScimark,cpu)
    for i in range(int(cmsScimark)):
        command= Commands[3]+" >& "+scimark.name
        print command+" [%s/%s]"%(i+1,int(cmsScimark))
        sys.stdout.flush()
        cmsScimarkstdout=os.popen4(command)[1].read()
        print cmsScimarkstdout
        sys.stdout.flush()
    print "Following with %s cmsScimarkLarge on cpu%s"%(cmsScimarkLarge,cpu)
    for i in range(int(cmsScimarkLarge)):
        command= Commands[3]+" -large >& "+scimarklarge.name
        print command+" [%s/%s]"%(i+1,int(cmsScimarkLarge))
        sys.stdout.flush()
        cmsScimarkLargestdout=os.popen4(command)[1].read()
        print cmsScimarkstdoutLarge
        sys.stdout.flush()
    #Stopping all cmsScimark jobs and analysing automatically the logfiles
    print "Stopping all cmsScimark jobs"
    print AuxiliaryScripts[2]
    sys.stdout.flush()
    cmsScimarkStopstdout=os.popen4(AuxiliaryScripts[2])[1].read()
    print cmsScimarkStopstdout
    sys.stdout.flush()

    #Create a tarball of the work directory
    TarFile=cmssw_version+"_"+host+"_"+user+".tar"
    tarcmd="tar -cvf "+TarFile+" *; gzip "+TarFile
    print tarcmd
    sys.stdout.flush()
    tarstdout=os.popen4(tarcmd)[1].read()
    print tarstdout
    sys.stdout.flush()
    #Archive it on CASTOR
    castorcmd="rfcp "+TarFile+".gz "+castordir+TarFile+".gz"
    print castorcmd
    sys.stdout.flush()
    castorstdout=os.popen4(castorcmd)[1].read()
    print castorstdout
    sys.stdout.flush()
    #End of script actions!

    #Print a time stamp at the end:
    date=time.ctime(time.time())
    print "Performance Suite finished running at %s on %s in directory %s" % (date,host,path)
    
if __name__ == "__main__":
    main(sys.argv[1:])

