import sys
import json
import os
import copy

def performInjectionOptionTest(opt):
    if opt.show:
        print 'Not injecting to wmagent in --show mode. Need to run the worklfows.'
        sys.exit(-1)
    if opt.wmcontrol=='init':
        #init means it'll be in test mode
        opt.nThreads=0
    if opt.wmcontrol=='test':
        #means the wf were created already, and we just dryRun it.
        opt.dryRun=True
    if opt.wmcontrol=='submit' and opt.nThreads==0:
        print 'Not injecting to wmagent in -j 0 mode. Need to run the worklfows.'
        sys.exit(-1)
    if opt.wmcontrol=='force':
        print "This is an expert setting, you'd better know what you're doing"
        opt.dryRun=True


class MatrixInjector(object):

    def __init__(self,mode='init'):
        self.count=1040
        self.testMode=((mode!='submit') and (mode!='force'))
        self.version =1

        #wagemt stuff
        self.wmagent=os.getenv('WMAGENT_REQMGR')
        if not self.wmagent:
            self.wmagent = 'cmsweb.cern.ch'
            
        #couch stuff
        self.couch = 'https://'+self.wmagent+'/couchdb'
        self.couchDB = 'reqmgr_config_cache'
        self.couchCache={} # so that we do not upload like crazy, and recyle cfgs
        self.user = os.getenv('USER')
        self.group = 'ppd'
        self.label = 'RelValSet_'+os.getenv('CMSSW_VERSION').replace('-','')+'_v'+str(self.version)


        if not os.getenv('WMCORE_ROOT'):
            print '\n\twmclient is not setup properly. Will not be able to upload or submit requests.\n'
            if not self.testMode:
                print '\n\t QUIT\n'
                sys.exit(-18)
        else:
            print '\n\tFound wmclient\n'
            
        self.defaultChain={
            "RequestType" :   "TaskChain",                    #this is how we handle relvals
            "AcquisitionEra": {},                             #Acq Era
            "Requestor": self.user,                           #Person responsible
            "Group": self.group,                              #group for the request
            "CMSSWVersion": os.getenv('CMSSW_VERSION'),       #CMSSW Version (used for all tasks in chain)
            "ScramArch": os.getenv('SCRAM_ARCH'),             #Scram Arch (used for all tasks in chain)
            "ProcessingVersion": self.version,                #Processing Version (used for all tasks in chain)
            "GlobalTag": None,                                #Global Tag (overridden per task)
            "CouchURL": self.couch,                           #URL of CouchDB containing Config Cache
            "CouchDBName": self.couchDB,                      #Name of Couch Database containing config cache
            #- Will contain all configs for all Tasks
            "SiteWhitelist" : ["T2_CH_CERN", "T1_US_FNAL"],   #Site whitelist
            "TaskChain" : None,                                  #Define number of tasks in chain.
            "nowmTasklist" : [],  #a list of tasks as we put them in
            "unmergedLFNBase" : "/store/unmerged",
            "mergedLFNBase" : "/store/relval",
            "dashboardActivity" : "relval",
            "Memory" : 2000,
            "SizePerEvent" : 1234,
            "TimePerEvent" : 20
            }

        self.defaultScratch={
            "TaskName" : None,                            #Task Name
            "ConfigCacheID" : None,                   #Generator Config id
            "GlobalTag": None,
            "SplittingAlgorithm"  : "EventBased",             #Splitting Algorithm
            "SplittingArguments" : {"events_per_job" : None},  #Size of jobs in terms of splitting algorithm
            "RequestNumEvents" : None,                      #Total number of events to generate
            "Seeding" : "AutomaticSeeding",                          #Random seeding method
            "PrimaryDataset" : None,                          #Primary Dataset to be created
            "nowmIO": {}
            }
        self.defaultInput={
            "TaskName" : "DigiHLT",                                      #Task Name
            "ConfigCacheID" : None,                                      #Processing Config id
            "GlobalTag": None,
            "InputDataset" : None,                                       #Input Dataset to be processed
            "SplittingAlgorithm"  : "LumiBased",                        #Splitting Algorithm
            "SplittingArguments" : {"lumis_per_job" : 10},               #Size of jobs in terms of splitting algorithm
            "nowmIO": {}
            }
        self.defaultTask={
            "TaskName" : None,                                 #Task Name
            "InputTask" : None,                                #Input Task Name (Task Name field of a previous Task entry)
            "InputFromOutputModule" : None,                    #OutputModule name in the input task that will provide files to process
            "ConfigCacheID" : None,                            #Processing Config id
            "GlobalTag": None,
            "SplittingAlgorithm"  : "LumiBased",                        #Splitting Algorithm
            "SplittingArguments" : {"lumis_per_job" : 10},               #Size of jobs in terms of splitting algorithm
            "nowmIO": {}
            }

        self.chainDicts={}


    def prepare(self,mReader, directories, mode='init'):
        
        for (n,dir) in directories.items():
            chainDict=copy.deepcopy(self.defaultChain)
            print "inspecting",dir
            nextHasDSInput=None
            for (x,s) in mReader.workFlowSteps.items():
                #x has the format (num, prefix)
                #s has the format (num, name, commands, stepList)
                if x[0]==n:
                    #print "found",n,s[3]
                    chainDict['RequestString']='RV'+chainDict['CMSSWVersion']+s[1].split('+')[0]
                    index=0
                    splitForThisWf=None
                    for step in s[3]:
                        if 'INPUT' in step or (not isinstance(s[2][index],str)):
                            nextHasDSInput=s[2][index]
                        else:
                            if 'HARVEST' in step:
                                continue
                            if (index==0):
                                #first step and not input -> gen part
                                chainDict['nowmTasklist'].append(copy.deepcopy(self.defaultScratch))
                                chainDict['nowmTasklist'][-1]['PrimaryDataset']='RelVal'+s[1].split('+')[0]
                                if not '--relval' in s[2][index]:
                                    print 'Impossible to create task from scratch'
                                    return -12
                                else:
                                    arg=s[2][index].split()
                                    ns=map(int,arg[arg.index('--relval')+1].split(','))
                                    chainDict['nowmTasklist'][-1]['RequestNumEvents'] = ns[0]
                                    chainDict['nowmTasklist'][-1]['SplittingArguments']['events_per_job'] = ns[1]
                            elif nextHasDSInput:
                                chainDict['nowmTasklist'].append(copy.deepcopy(self.defaultInput))
                                chainDict['nowmTasklist'][-1]['InputDataset']=nextHasDSInput.dataSet
                                splitForThisWf=nextHasDSInput.split
                                chainDict['nowmTasklist'][-1]['SplittingArguments']['lumis_per_job']=splitForThisWf
                                # get the run numbers or #events
                                if len(nextHasDSInput.run):
                                    chainDict['nowmTasklist'][-1]['RunWhitelist']=nextHasDSInput.run
                                nextHasDSInput=None
                            else:
                                #not first step and no inputDS
                                chainDict['nowmTasklist'].append(copy.deepcopy(self.defaultTask))
                                if splitForThisWf:
                                    chainDict['nowmTasklist'][-1]['SplittingArguments']['lumis_per_job']=splitForThisWf
                            #print step
                            chainDict['nowmTasklist'][-1]['TaskName']=step
                            try:
                                chainDict['nowmTasklist'][-1]['nowmIO']=json.loads(open('%s/%s.io'%(dir,step)).read())
                            except:
                                print "Failed to find",'%s/%s.io'%(dir,step),".The workflows were probably not run on cfg not created"
                                return -15
                            chainDict['nowmTasklist'][-1]['ConfigCacheID']='%s/%s.py'%(dir,step)
                            chainDict['nowmTasklist'][-1]['GlobalTag']=chainDict['nowmTasklist'][-1]['nowmIO']['GT'] # copy to the proper parameter name
                            chainDict['GlobalTag']=chainDict['nowmTasklist'][-1]['nowmIO']['GT'] #set in general to the last one of the chain
                            #chainDict['AcquisitionEra']=(chainDict['CMSSWVersion']+'-'+chainDict['GlobalTag']).replace('::All','')
                            chainDict['AcquisitionEra'][step]=(chainDict['CMSSWVersion']+'-'+chainDict['nowmTasklist'][-1]['GlobalTag']).replace('::All','')
                        index+=1
                        
            #wrap up for this one
            #print 'wrapping up'
            chainDict['TaskChain']=len(chainDict['nowmTasklist'])
            #loop on the task list
            for i_second in reversed(range(len(chainDict['nowmTasklist']))):
            #for t_second in reversed(chainDict['nowmTasklist']):
                t_second=chainDict['nowmTasklist'][i_second]
                #print "t_second taskname", t_second['TaskName']
                if 'primary' in t_second['nowmIO']:
                    #print t_second['nowmIO']['primary']
                    primary=t_second['nowmIO']['primary'][0].replace('file:','')
                    #for t_input in reversed(chainDict['nowmTasklist']):
                    for i_input in reversed(range(0,i_second)):
                        t_input=chainDict['nowmTasklist'][i_input]
                        for (om,o) in t_input['nowmIO'].items():
                            if primary in o:
                                #print "found",primary,"procuced by",om,"of",t_input['TaskName']
                                t_second['InputTask'] = t_input['TaskName']
                                t_second['InputFromOutputModule'] = om
                                #print 't_second',t_second
                                break

            ## there is in fact only one acquisition era
            ##if len(set(chainDict['AcquisitionEra'].values()))==1:
            chainDict['AcquisitionEra'] = chainDict['AcquisitionEra'].values()[0]
                
            ## clean things up now
            for (i,t) in enumerate(chainDict['nowmTasklist']):
                t.pop('nowmIO')
                chainDict['Task%d'%(i+1)]=t

            chainDict.pop('nowmTasklist')
            self.chainDicts[n]=chainDict

            
        return 0

    def uploadConf(self,filePath,label,where):
        labelInCouch=self.label+'_'+label
        cacheName=filePath.split('/')[-1]
        if self.testMode:
            self.count+=1
            print '\tFake upload to couch with label',labelInCouch
            return self.count
        else:
            try:
                from modules.wma import upload_to_couch
            except:
                print '\n\tUnable to find wmcontrol modules. Please include it in your python path\n'
                print '\n\t QUIT\n'
                sys.exit(-16)
            if cacheName in self.couchCache:
                print "Not re-uploading",filePath,"to",where,"for",label
                cacheId=self.couchCache[cacheName]
            else:
                print "Loading",filePath,"to",where,"for",label
                cacheId=upload_to_couch(filePath,
                                        labelInCouch,
                                        self.user,
                                        self.group,
                                        test_mode=False,
                                        url=where
                                        )
                self.couchCache[cacheName]=cacheId
            return cacheId
    
    def upload(self):
        for (n,d) in self.chainDicts.items():
            for it in d:
                if it.startswith("Task") and it!='TaskChain':
                    #upload
                    couchID=self.uploadConf(d[it]['ConfigCacheID'],
                                            str(n)+d[it]['TaskName'],
                                            d['CouchURL']
                                            )
                    print d[it]['ConfigCacheID']," uploaded to couchDB for",str(n),"with ID",couchID
                    d[it]['ConfigCacheID']=couchID
            
    def submit(self):
        try:
            from modules.wma import makeRequest,approveRequest
            from wmcontrol import random_sleep
            print '\n\tFound wmcontrol\n'
        except:
            print '\n\tUnable to find wmcontrol modules. Please include it in your python path\n'
            if not self.testMode:
                print '\n\t QUIT\n'
                sys.exit(-17)

        import pprint
        for (n,d) in self.chainDicts.items():
            if self.testMode:
                print "Only viewing request",n
                print pprint.pprint(d)
            else:
                #submit to wmagent each dict
                print "For eyes before submitting",n
                print pprint.pprint(d)
                print "Submitting",n,"..........."
                workFlow=makeRequest(self.wmagent,d,encodeDict=True)
                approveRequest(self.wmagent,workFlow)
                print "...........",n,"submitted"
                random_sleep()
            

        
