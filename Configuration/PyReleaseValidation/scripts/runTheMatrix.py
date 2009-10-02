#!/usr/bin/env python

import os, sys, re, time

import random
from threading import Thread
        
class WorkFlowRunner(Thread):
    def __init__(self, wf):
        Thread.__init__(self)
        self.wf = wf

        self.status=-1
        self.report=''
        self.nfail=0
        self.npass=0

        return

    def doCmd(self, cmd, dryRun=False):

        msg = "\n# in: " +os.getcwd()
        if dryRun: msg += " dryRun for '"
        else:      msg += " going to execute "
        msg += cmd.replace(';','\n')
        print msg

        cmdLog = open(self.wf.numId+'_'+self.wf.nameId+'/cmdLog','a')
        cmdLog.write(msg+'\n')
        cmdLog.close()
        
        ret = 0
        if not dryRun:
            ret = os.system(cmd)
            if ret != 0:
                print "ERROR executing ",cmd,'ret=', ret

        return ret
    
    def run(self):

        startDir = os.getcwd()

        wfDir = self.wf.numId+'_'+self.wf.nameId
        if not os.path.exists(wfDir):
            os.makedirs(wfDir)
        # os.chdir(wfDir)

        # clean up first:
        cmd = 'rm -f step*.py *.root;'
        self.doCmd(cmd)

        preamble = ''
        if os.path.exists( os.path.join(os.environ["CMS_PATH"],'cmsset_default.sh') ) :
            preamble = 'source $CMS_PATH/cmsset_default.sh; '
        else:
            preamble = 'source $CMS_PATH/sw/cmsset_default.sh; '
        preamble += 'eval `scram run -sh`; '
        preamble += 'cd '+wfDir+'; '
        
        startime='date %s' %time.asctime()

        # set defaults for the statuses
        stat1 = 'PASSED'
        stat2 = 'PASSED' 
        stat3 = 'PASSED'
        stat4 = 'PASSED'
        if not self.wf.cmdStep2: stat2 = 'NOSTEP'
        if not self.wf.cmdStep3: stat3 = 'NOSTEP'
        if not self.wf.cmdStep4: stat4 = 'NOSTEP'
        
        # run the first workflow:
        cmd = preamble
        cmd += self.wf.cmdStep1 + ' --fileout file:raw.root '
        cmd += ' > %s 2>&1; ' % ('step1_'+self.wf.nameId+'.log ',)
        retStep1 = self.doCmd(cmd)
        print " ... ret: " , retStep1

        # prepare and run the next workflows -- if the previous step was OK :
        # set some defaults
        retStep2 = 0
        retStep3 = 0
        retStep4 = 0
        if self.wf.cmdStep2 and retStep1 == 0:
            fullcmd = preamble
            fullcmd += self.wf.cmdStep2 + ' --filein file:raw.root --fileout file:reco.root '
            fullcmd += ' > %s 2>&1; ' % ('step2_'+self.wf.nameId+'.log ',)
            # print fullcmd
            retStep2 = self.doCmd(fullcmd)
#            if random.randint(0,100) < 20 : retStep2 = -42

            if self.wf.cmdStep3 and retStep2 == 0:
                fullcmd = preamble
                fullcmd += self.wf.cmdStep3 + ' --filein file:reco.root --fileout file:step3.root '
                fullcmd += ' > %s 2>&1; ' % ('step3_'+self.wf.nameId+'.log ',)
                # print fullcmd
                retStep3 = self.doCmd(fullcmd)
#                if random.randint(0,100) < 40 : retStep3 = -42
                if self.wf.cmdStep4 and retStep3 == 0:
                    fullcmd = preamble
                    fullcmd += self.wf.cmdStep4 + ' --filein file:step3.root '
                    fullcmd += ' > %s 2>&1; ' % ('step4_'+self.wf.nameId+'.log ',)
                    # print fullcmd
                    retStep4 = self.doCmd(fullcmd)
#                    if random.randint(0,100) < 40 : retStep4 = -42

        os.chdir(startDir)

        endtime='date %s' %time.asctime()
        tottime='%s-%s'%(endtime,startime)

        self.nfail = [0,0,0,0]
        self.npass = [1,1,1,1]
        if 'NOSTEP' in stat2: # don't say reco/alca is passed if we don't have to run them
            self.npass = [1,0,0,0]
        else: # we have a reco step, check for alca:
            if 'NOSTEP' in stat3 :
                self.npass = [1,1,0,0]
                if 'NOSTEP' in stat4 :
                    self.npass = [1,1,1,0]
        if retStep1 != 0 :
            stat1 = 'FAILED'
            stat2 = 'NOTRUN'
            stat3 = 'NOTRUN'
            stat4 = 'NOTRUN'
            self.npass = [0,0,0,0]
            self.nfail = [1,0,0,0]

        if retStep2 != 0 :
            stat2 = 'FAILED'
            stat3 = 'NOTRUN'
            stat4 = 'NOTRUN'
            self.npass = [1,0,0,0]
            self.nfail = [0,1,0,0]

        if retStep3 != 0 :
            stat3 = 'FAILED'
            stat4 = 'NOTRUN'
            self.npass = [1,1,0,0]
            self.nfail = [0,0,1,0]

        if retStep4 != 0 :
            stat4 = 'FAILED'
            self.npass = [1,1,1,0]
            self.nfail = [0,0,0,1]

        logStat = 'Step1-'+stat1+' Step2-'+stat2+' Step3-'+stat3+' '+' Step4-'+stat4+' ' 
        self.report+='%s_%s %s - time %s; exit: %s %s %s %s \n' % (self.wf.numId, self.wf.nameId, logStat, tottime, retStep1,retStep2,retStep3, retStep4)

        return


# ================================================================================

class WorkFlow(object):

    def __init__(self, num, nameID, cmd1, cmd2=None, cmd3=None, cmd4=None, real=None):
        self.numId  = num.strip()
        self.nameId = nameID
        self.cmdStep1 = cmd1
        if self.cmdStep1: self.cmdStep1 = self.cmdStep1.replace('--no_exec', '') # make sure the commands execute
        self.cmdStep2 = cmd2
        if self.cmdStep2: self.cmdStep2 = self.cmdStep2.replace('--no_exec', '') # make sure the commands execute
        self.cmdStep3 = cmd3
        if self.cmdStep3: self.cmdStep3 = self.cmdStep3.replace('--no_exec', '') # make sure the commands execute
        self.cmdStep4 = cmd4
        if self.cmdStep4: self.cmdStep4 = self.cmdStep4.replace('--no_exec', '') # make sure the commands execute

        # run on real data requested:
        self.real = real
        return

# ================================================================================

class MatrixReader(object):

    def __init__(self):

        self.reset()

        return

    def reset(self):

        self.step1WorkFlows = {}
        self.step2WorkFlows = {}
        self.step3WorkFlows = {}
        self.step4WorkFlows = {}

        self.workFlows = []
        self.nameList  = {}
        
        return

    def readMatrix(self, fileNameIn, prefix='', offset=0):
        
        print "processing ", fileNameIn
        lines = []
        try:
            inFile = open(fileNameIn, 'r')
            lines = inFile.readlines()
            inFile.close()
        except Exception, e:
            print "ERROR reading in file ", fileNameIn, str(e)
            return
        
        realRe = re.compile('\s*([1-9][0-9]*\.*[0-9]*)\s*\+\+\s*(.*?)\s*\+\+\s*(.*?)\s*\+\+\s*(.*?)\s*@@@\s*(.*)\s*')
        step1Re = re.compile('\s*([1-9][0-9]*\.*[0-9]*)\s*\+\+\s*(.*?)\s*\+\+\s*(.*?)\s*@@@\s*(.*)\s*')
        step2Re = re.compile('\s*STEP2\s*\+\+\s*(\S*)\s*@@@\s*(.*)\s*')
        step3Re = re.compile('\s*STEP3\s*\+\+\s*(\S*)\s*@@@\s*(.*)\s*')
        step4Re = re.compile('\s*STEP4\s*\+\+\s*(\S*)\s*@@@\s*(.*)\s*')
        for lineIn in lines:
            line = lineIn.strip()

            realMatch = realRe.match(line)
            if realMatch :
                continue
##-not now                num  = realMatch.group(1).strip()
##-not now                name = realMatch.group(2).strip().replace('<','').replace('>','').replace(':','')
##-not now                next = realMatch.group(3).strip().replace('+','').replace(',', ' ')
##-not now                real = realMatch.group(4).strip()
##-not now                dummy, inFile, run, runNr = real.split()
##-not now                cmd  = ""
##-not now
##-not now                # for run on real data. use:
##-not now                cmdReal = 'cmsDriver.py step2 '
##-not now                cmdReal += '-s RAW2DIGI,RECO:reconstructionCosmics_woDeDx,ALCA:MuAlCalIsolatedMu+RpcCalHLT '
##-not now                cmdReal += '--relval 25000,100 --datatier RECO --eventcontent RECO '
##-not now                cmdReal += '--conditions FrontierConditions_GlobalTag,CRAFT_30X::All '
##-not now                cmdReal += '--scenario cosmics --data '
##-not now                cmdReal += '--fileIn file:'
##-not now                # ... and the file/run from the line ...
##-not now
##-not now                step2 = "None"
##-not now                step3 = "None"
##-not now                try:
##-not now                    step2, step3 = next.split()
##-not now                    step2 = step2.strip()
##-not now                    step3 = step3.strip()
##-not now                except ValueError:
##-not now                    if len(next) > 0:
##-not now                        step2 = next.strip()
##-not now                    pass
##-not now                self.step1WorkFlows[float(num)] = (num, name, step2, step3, cmd, real)
##-not now                continue
                
            step1Match = step1Re.match(line)
            if step1Match :
                num  = step1Match.group(1).strip()
                name = step1Match.group(2).strip().replace('<','').replace('>','').replace(':','')
                next = step1Match.group(3).strip().replace('+','').replace(',', ' ')
                cmd  = step1Match.group(4).strip()
                step2 = "None"
                step3 = "None"
                step4 = "None"

                steps = next.split()
                if len(steps) > 0:
                    step2 = steps[0].strip()
                if len(steps) > 1:
                    step3 = steps[1].strip()
                if len(steps) > 2:
                    step4 = steps[2].strip()
                
                self.step1WorkFlows[float(num)+offset] = (str(float(num)+offset), name, step2, step3, step4, cmd, None)
                continue
            
            step2Match = step2Re.match(line)
            if step2Match :
                name = step2Match.group(1).strip()
                cmd  = step2Match.group(2).strip()
                self.step2WorkFlows[name] = (cmd.replace('--no_exec','') ) # make sure the command is really run
                continue

            step3Match = step3Re.match(line)
            if step3Match :
                name = step3Match.group(1).strip()
                cmd  = step3Match.group(2).strip()
                self.step3WorkFlows[name] = ( cmd.replace('--no_exec','') ) # make sure the command is really run
                continue

            step4Match = step4Re.match(line)
            if step4Match :
                name = step4Match.group(1).strip()
                cmd  = step4Match.group(2).strip()
                self.step4WorkFlows[name] = ( cmd.replace('--no_exec','') ) # make sure the command is really run
                continue

        return

    def showRaw(self):

        print "found ", len(self.step1WorkFlows.keys()), ' workflows for step1:'
        ids = self.step1WorkFlows.keys()
        ids.sort()
        for key in ids:
            val = self.step1WorkFlows[key]
            print key, ':', val
        
        print "found ", len(self.step2WorkFlows.keys()), ' workflows for step2:'
        for key, val in self.step2WorkFlows.items():
            print key, ':', val
        
        print "found ", len(self.step3WorkFlows.keys()), ' workflows for step3:'
        for key, val in self.step3WorkFlows.items():
            print key, ':', val
        
        print "found ", len(self.step4WorkFlows.keys()), ' workflows for step4:'
        for key, val in self.step4WorkFlows.items():
            print key, ':', val
        
        return

    def showWorkFlows(self, selected=None):

        print "found ", len(self.workFlows), ' workflows:'
        if selected:
            print "   of which the following", len(selected), 'were selected.'
            print selected
            
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        maxLen = 100
        for wf in self.workFlows:
            if selected and float(wf.numId) not in selected: continue
            n1+=1
            print "%-6s %-35s [1]: %s ..." % (wf.numId, wf.nameId, wf.cmdStep1[:maxLen])
            if wf.cmdStep2:
                n2+=1
                print "       %35s [2]: %s ..." % ( ' ', wf.cmdStep2[:maxLen])
                if wf.cmdStep3:
                    n3+=1
                    print "       %35s [3]: %s ..." % ( ' ', wf.cmdStep3[:maxLen])
                    if wf.cmdStep4:
                        n4+=1
                        print "       %35s [4]: %s ..." % ( ' ', wf.cmdStep4[:maxLen])

        print n1, 'workflows for step1,'
        print n2, 'workflows for step1 + step2,'
        print n3, 'workflows for step1 + step2 + step3'
        print n4, 'workflows for step1 + step2 + step3 + step4'

        return
    
    def createWorkFlows(self, prefix=''):

        ids = self.step1WorkFlows.keys()
        ids.sort()
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        for key in ids:
            val = self.step1WorkFlows[key]
            num, name, step2, step3, step4, cmd, real = val
            nameId = num+'_'+name
            if step2.lower() != 'none':
                name += '+'+step2
                if step3.lower() != 'none':
                    name += '+'+step3
                    if step4.lower() != 'none':
                        name += '+'+step4
            if nameId in self.nameList.keys():
                print "==> duplicate name found for ", nameId
                print '    keeping  : ', self.nameList[nameId]
                print '    ignoring : ', val
            else:
                self.nameList[nameId] = val

            cmd2 = None
            cmd3 = None
            cmd4 = None
            
            n1 += 1

            if step2.lower() != 'none':
                n2 += 1
                cmd2 = self.step2WorkFlows[step2]
                if step3.lower() != 'none':
                    n3 += 1
                    cmd3 = self.step3WorkFlows[step3]
                    if step4.lower() != 'none':
                        n4 += 1
                        cmd4 = self.step4WorkFlows[step4]
                    #print '\tstep3 : ', self.step3WorkFlows[step3]
            self.workFlows.append( WorkFlow(num, name, cmd, cmd2, cmd3, cmd4) )

        return

    def show(self, selected=None):
        # self.showRaw()
        self.showWorkFlows(selected)
        print '\n','-'*80,'\n'


    def updateDB(self):

        import pickle
        pickle.dump(self.workFlows, open('theMatrix.pkl', 'w') )

        return

# ================================================================================

class MatrixRunner(object):

    def __init__(self, wfIn=None, nThrMax=8):

        self.workFlows = wfIn

        self.threadList = []
        self.maxThreads = int(nThrMax) # make sure we get a number ...


    def activeThreads(self):

        nActive = 0
        for t in self.threadList:
            if t.isAlive() : nActive += 1

        return nActive

        
    def runTests(self, testList=None):

        startDir = os.getcwd()

    	# make sure we have a way to set the environment in the threads ...
    	if not os.environ.has_key('CMS_PATH'):
    	    cmsPath = '/afs/cern.ch/cms'
    	    print "setting default for CMS_PATH to", cmsPath
    	    os.environ['CMS_PATH'] = cmsPath

    	report=''    	
    	print 'Running in %s thread(s)' % self.maxThreads
                
        for wf in self.workFlows:

            if testList and float(wf.numId) not in [float(x) for x in testList]: continue

            # don't know yet how to treat real data WF ...
            if wf.real : continue

            item = wf.nameId
            if os.path.islink(item) : continue # ignore symlinks
            
    	    # make sure we don't run more than the allowed number of threads:
    	    while self.activeThreads() >= self.maxThreads:
    	        time.sleep(10)
                continue
    	    
    	    print '\nPreparing to run %s %s' % (wf.numId, item)
          
##            if testList: # if we only run a selection, run only 5 events instead of 10
##                wf.cmdStep1 = wf.cmdStep1.replace('-n 10', '-n 5')
                
    	    current = WorkFlowRunner(wf)
    	    self.threadList.append(current)
    	    current.start()
            time.sleep(random.randint(1,5)) # try to avoid race cond by sleeping random amount of time [1,5] sec 

    	# wait until all threads are finished
        while self.activeThreads() > 0:
    	    time.sleep(5)
    	    
    	# all threads are done now, check status ...
    	nfail1 = 0
    	nfail2 = 0
        nfail3 = 0
        nfail4 = 0
    	npass  = 0
        npass1 = 0
        npass2 = 0
        npass3 = 0
        npass4 = 0
    	for pingle in self.threadList:
    	    pingle.join()
            nfail1 += pingle.nfail[0]
            nfail2 += pingle.nfail[1]
            nfail3 += pingle.nfail[2]
            nfail4 += pingle.nfail[3]
    	    npass1 += pingle.npass[0]
    	    npass2 += pingle.npass[1]
    	    npass3 += pingle.npass[2]
    	    npass4 += pingle.npass[3]
    	    npass  += npass1+npass2+npass3+npass4
    	    report += pingle.report
    	    # print pingle.report
    	    
    	report+='\n %s %s %s %s tests passed, %s %s %s %s failed\n' %(npass1, npass2, npass3, npass4, nfail1, nfail2, nfail3, nfail4)
    	print report
    	
    	runall_report_name='runall-report-step123-.log'
    	runall_report=open(runall_report_name,'w')
    	runall_report.write(report)
    	runall_report.close()

        os.chdir(startDir)
    	
    	return

        
# ================================================================================

def runSelected(testList, nThreads=4, show=False) :

    stdList = ['5.2', # SingleMu10 FastSim
               '7',   # Cosmics+RECOCOS+ALCACOS
               '8',   # BeamHalo+RECOCOS+ALCABH
#              '24',  # TTbar+RECO1+ALCATT1 IDEAL
               '25',  # TTbar+RECO2+ALCATT2  STARTUP
               ]
    hiStatList = [
#                 '15',  # SingleMuPt10
#                 '119', # ZTT+RECO1
                  '123.3', # TTBar FastSim
                   ]

    mrd = MatrixReader()
    files = ['cmsDriver_standard_hlt.txt', 'cmsDriver_highstats_hlt.txt']
    offset = 0
    for matrixFile in files:
        try:
            mrd.readMatrix(matrixFile, offset=offset)
        except Exception, e:
            print "ERROR reading file:", matrixFile, str(e)
        offset += 100

    try:
        mrd.createWorkFlows()
    except Exception, e:
        print "ERROR creating workflows :", str(e)

    if testList == []:
        testList = stdList+hiStatList

    ret = 0
    if show:
        mrd.show([float(x) for x in testList])
        print 'selected items:', testList
    else:
        mRunnerHi = MatrixRunner(mrd.workFlows, nThreads)
        ret = mRunnerHi.runTests(testList)

    return ret

# --------------------------------------------------------------------------------

def runAll(testList=None, nThreads=4, show=False) :

    mrd = MatrixReader()
    files = ['cmsDriver_standard_hlt.txt', 'cmsDriver_highstats_hlt.txt']
    offset = 0
    for matrixFile in files:
        try:
            mrd.readMatrix(matrixFile, offset=offset)
        except Exception, e:
            print "ERROR reading file:", matrixFile, str(e)
        offset += 100
        
    try:
        mrd.createWorkFlows()
    except Exception, e:
        print "ERROR creating workflows : "+str(e)


    ret = 0
    
    if show:
        mrd.show()
        print "nThreads = ",nThreads
    else:
        mRunnerHi = MatrixRunner(mrd.workFlows, nThreads)
        ret = mRunnerHi.runTests()

    return ret


# --------------------------------------------------------------------------------

def runOnly(only, show, nThreads=4):

    if not only: return
    
    for what in only:
        print "found request to run relvals only for ",what


# ================================================================================

if __name__ == '__main__':

    import getopt
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "j:sl:nqo:", ["nproc=",'selected','list=','showMatrix','only='])
    except getopt.GetoptError, e:
        print "unknown option", str(e)
        sys.exit(2)
        
# check command line parameter

    np=4 # default: four threads
    sel = None
    show = False
    only = None
    for opt, arg in opts :
        if opt in ('-j', "--nproc" ):
            np=int(arg)
        if opt in ('-n','-q','--showMatrix', ):
            show = True
        if opt in ('-s','--selected',) :
            sel = []
        if opt in ('-o','--only',) :
            only = []
        if opt in ('-l','--list',) :
            sel = arg.split(',')

    # print "sel",sel
    ret = 0
    if sel != None: # explicit distinguish from empty list (which is also false)
        ret = runSelected(testList=sel, nThreads=np, show=show)
    elif only != None:
        ret = runOnly(only=only, show=show, nThreads=np)
    else:
        ret = runAll(show=show, nThreads=np)

    sys.exit(ret)
